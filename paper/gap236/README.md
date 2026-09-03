# Conditional Gap 236 manuscript

This directory contains the source of the working paper
*Bounded gaps between primes: towards the bound* $H_1\leq236$.

The manuscript deliberately states a conditional assembly theorem. It must not be
presented as an unconditional improvement while the D27 boundary certificate, the
analytic soundness bridge, the complete Type-IIc proof review, and the shaped-support
sieve-to-`DHL48` formalization remain open.

Build the PDF with:

```bash
cd paper/gap236
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Generated LaTeX files and the PDF are build products and are not tracked.

