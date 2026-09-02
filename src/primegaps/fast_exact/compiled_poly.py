"""Compiled exact bivariate polynomial operations via FLINT univariate encoding."""

from __future__ import annotations


def _numerator_denominator(value):
    """Return Python integers for Fraction, mpq, or Sage rational values."""
    numerator = value.numerator
    denominator = value.denominator
    if callable(numerator):
        numerator = numerator()
    if callable(denominator):
        denominator = denominator()
    return int(numerator), int(denominator)


class FlintEncodedPolynomialBackend:
    """Encode x^a z^b as q^(a*stride+b) in a FLINT fmpq polynomial.

    The stride must exceed every z exponent that can occur in a product, so
    ordinary univariate exponent addition exactly represents bivariate
    multiplication without carries.
    """

    def __init__(self, *, stride: int = 256, rational):
        if stride < 2:
            raise ValueError("stride must be at least two")
        from sage.all import PolynomialRing, QQ

        self.stride = int(stride)
        self.rational = rational
        self.base_ring = QQ
        self.ring = PolynomialRing(QQ, "q")

    def _coefficient(self, value):
        numerator, denominator = _numerator_denominator(value)
        return self.base_ring(numerator) / denominator

    def scalar_zero(self):
        return self.rational(0)

    def scalar_add(self, left, right):
        return left + right

    def from_dict(self, polynomial):
        encoded = {}
        for (x_power, z_power), coefficient in polynomial.items():
            if not 0 <= z_power < self.stride:
                raise ValueError(
                    f"z exponent {z_power} exceeds stride {self.stride}"
                )
            encoded[x_power * self.stride + z_power] = self._coefficient(
                coefficient
            )
        return self.ring(encoded)

    def zero(self):
        return self.ring.zero()

    def multiply(self, left, right):
        return left * right

    def add_scaled(self, target, source, scalar):
        return target + int(scalar) * source

    def contract_bilinear(self, left, right, moments):
        """Evaluate an exact moment-Hankel bilinear form inside Sage/FLINT."""
        from sage.all import matrix, vector

        left_exponents = tuple(left)
        right_exponents = tuple(right)
        left_vector = vector(
            self.base_ring,
            [self._coefficient(left[exponent]) for exponent in left_exponents],
        )
        right_vector = vector(
            self.base_ring,
            [self._coefficient(right[exponent]) for exponent in right_exponents],
        )
        hankel = matrix(
            self.base_ring,
            len(left_exponents),
            len(right_exponents),
            [
                self._coefficient(
                    moments[(left_x + right_x, left_z + right_z)]
                )
                for left_x, left_z in left_exponents
                for right_x, right_z in right_exponents
            ],
        )
        value = left_vector * hankel * right_vector
        return self.rational(
            int(value.numerator()), int(value.denominator())
        )

    def terms(self, polynomial):
        for encoded_power, coefficient in polynomial.dict().items():
            x_power, z_power = divmod(int(encoded_power), self.stride)
            yield (
                x_power,
                z_power,
                self.rational(
                    int(coefficient.numerator()),
                    int(coefficient.denominator()),
                ),
            )

    def integrate_product(
        self,
        density,
        candidate,
        *,
        kind,
        cell,
        verifier,
    ):
        if not density or not candidate:
            return self.rational(0)
        integrand = density * candidate
        if kind == "polygons":
            polygon, _sample = cell
            answer = self.rational(0)
            for x_power, z_power, coefficient in self.terms(integrand):
                answer += coefficient * verifier._polygon_monomial_moment(
                    x_power, z_power, polygon
                )
            return answer
        if kind in ("xintervals", "zintervals"):
            start, end, _sample = cell
            variable_is_x = kind == "xintervals"
            answer = self.rational(0)
            for x_power, z_power, coefficient in self.terms(integrand):
                if variable_is_x and z_power == 0:
                    degree = x_power
                elif not variable_is_x and x_power == 0:
                    degree = z_power
                else:
                    continue
                answer += coefficient * (
                    end ** (degree + 1) - start ** (degree + 1)
                ) / (degree + 1)
            return answer
        if kind == "point":
            coefficient = integrand[0]
            return self.rational(
                int(coefficient.numerator()), int(coefficient.denominator())
            )
        if kind == "empty":
            return self.rational(0)
        raise ValueError(f"unknown geometry kind {kind!r}")


class FlintModularEncodedPolynomialBackend:
    """FLINT polynomial contractions over one prime field.

    This uses the same carry-free bivariate-to-univariate encoding as the
    rational backend, but moves every large polynomial multiply and
    accumulation into ``F_p[q]``.  Only the small geometry moments are formed
    as rationals before being reduced modulo ``p``.
    """

    def __init__(self, prime: int, *, stride: int = 256):
        if stride < 2:
            raise ValueError("stride must be at least two")
        from sage.all import GF, PolynomialRing

        self.prime = int(prime)
        self.stride = int(stride)
        self.base_ring = GF(self.prime)
        self.ring = PolynomialRing(self.base_ring, "q")

    def _coefficient(self, value):
        numerator, denominator = _numerator_denominator(value)
        denominator %= self.prime
        if not denominator:
            raise ZeroDivisionError(
                f"coefficient denominator vanishes modulo {self.prime}"
            )
        return self.base_ring(numerator) / self.base_ring(denominator)

    def from_dict(self, polynomial):
        encoded = {}
        for (x_power, z_power), coefficient in polynomial.items():
            if not 0 <= z_power < self.stride:
                raise ValueError(
                    f"z exponent {z_power} exceeds stride {self.stride}"
                )
            encoded[x_power * self.stride + z_power] = self._coefficient(
                coefficient
            )
        return self.ring(encoded)

    def zero(self):
        return self.ring.zero()

    def scalar_zero(self):
        return 0

    def scalar_add(self, left, right):
        return (int(left) + int(right)) % self.prime

    def multiply(self, left, right):
        return left * right

    def add_scaled(self, target, source, scalar):
        return target + self.base_ring(int(scalar)) * source

    def terms(self, polynomial):
        for encoded_power, coefficient in polynomial.dict().items():
            x_power, z_power = divmod(int(encoded_power), self.stride)
            yield x_power, z_power, coefficient

    def integrate_product(
        self,
        density,
        candidate,
        *,
        kind,
        cell,
        verifier,
    ):
        if not density or not candidate:
            return 0
        integrand = density * candidate
        answer = self.base_ring.zero()
        if kind == "polygons":
            polygon, _sample = cell
            for x_power, z_power, coefficient in self.terms(integrand):
                moment = verifier._polygon_monomial_moment(
                    x_power, z_power, polygon
                )
                answer += coefficient * self._coefficient(moment)
            return int(answer)
        if kind in ("xintervals", "zintervals"):
            start, end, _sample = cell
            variable_is_x = kind == "xintervals"
            for x_power, z_power, coefficient in self.terms(integrand):
                if variable_is_x and z_power == 0:
                    degree = x_power
                elif not variable_is_x and x_power == 0:
                    degree = z_power
                else:
                    continue
                moment = (
                    end ** (degree + 1) - start ** (degree + 1)
                ) / (degree + 1)
                answer += coefficient * self._coefficient(moment)
            return int(answer)
        if kind == "point":
            return int(integrand[0])
        if kind == "empty":
            return 0
        raise ValueError(f"unknown geometry kind {kind!r}")


class ProductPolynomialBackend:
    """Evaluate several polynomial coefficient rings in one geometry pass.

    This is primarily for batching CRT primes. Exact density construction,
    support-cell enumeration, slice construction, and cached geometry moments
    are shared; only ring conversion, multiplication, and accumulation are
    repeated for each component backend.
    """

    def __init__(self, backends):
        self.backends = tuple(backends)
        if not self.backends:
            raise ValueError("at least one component backend is required")

    def from_dict(self, polynomial):
        return tuple(backend.from_dict(polynomial) for backend in self.backends)

    def zero(self):
        return tuple(backend.zero() for backend in self.backends)

    def scalar_zero(self):
        return tuple(backend.scalar_zero() for backend in self.backends)

    def scalar_add(self, left, right):
        return tuple(
            backend.scalar_add(left_value, right_value)
            for backend, left_value, right_value in zip(
                self.backends, left, right
            )
        )

    def multiply(self, left, right):
        return tuple(
            backend.multiply(left_value, right_value)
            for backend, left_value, right_value in zip(
                self.backends, left, right
            )
        )

    def add_scaled(self, target, source, scalar):
        return tuple(
            backend.add_scaled(target_value, source_value, scalar)
            for backend, target_value, source_value in zip(
                self.backends, target, source
            )
        )

    def integrate_product(self, density, candidate, **kwargs):
        return tuple(
            backend.integrate_product(
                density_value, candidate_value, **kwargs
            )
            for backend, density_value, candidate_value in zip(
                self.backends, density, candidate
            )
        )
