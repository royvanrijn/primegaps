from __future__ import annotations

import argparse

from .basis import symmetric_basis
from .support import contains, sample_uniform_simplex, stadlmann_240_parameters


def main() -> None:
    parser = argparse.ArgumentParser(description="Stadlmann support/basis diagnostic")
    parser.add_argument("--samples", type=int, default=200_000)
    parser.add_argument("--seed", type=int, default=260831126)
    args = parser.parse_args()

    p = stadlmann_240_parameters()
    outer = p.A[-1] + p.epsilon
    print("Published H1<=240 support")
    print(f"  epsilon={p.epsilon} delta={p.delta} A={p.A}")
    print(f"  B[1,1:4]={p.B[0][:4]} outer total={outer}")

    for degree in (21, 27):
        print(f"basis D={degree}: {len(symmetric_basis(degree, k=49))} functions")

    print("\nUniform-simplex geometry diagnostic (NOT a sieve score):")
    for k in (49, 48, 47, 46):
        points = sample_uniform_simplex(k, outer, args.samples, args.seed + k)
        acceptance = contains(points, p).mean()
        print(f"  k={k:2d}: support acceptance {acceptance:.8f}")

    print("\nNext required milestone: exact Section 5 recurrences -> M1/M2 -> reproduce k=49 eigenvalue > 1.")


if __name__ == "__main__":
    main()
