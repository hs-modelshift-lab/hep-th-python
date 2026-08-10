# hep-th-python

Consistency-check scripts accompanying manuscripts by Hergen Scheck
([ORCID 0000-0003-1806-4048](https://orcid.org/0000-0003-1806-4048)).

## `composition_check.py`

Independent verification of the algebraic identities used in *Spinor
Transformations and the Free Dirac Equation from Probability Bookkeeping*.

No result of that manuscript depends on this script. The identities it
evaluates are elementary and are stated in the text; the script is a check on
the algebra rather than on the claim.

The computation is exact rather than numerical: coefficients are Gaussian
rationals and expressions are polynomials in the half-angle symbols
`ca, sa, cb, sb, ch, sh`, reduced by `sa^2 -> 1-ca^2`, `sb^2 -> 1-cb^2`,
`sh^2 -> ch^2-1`. The checks are therefore polynomial identities. Two of the
sections additionally sample numerically over random angles.

What is checked, section by section:

1. **Generators** — unit determinant, unitarity, and the Axiom 1 transport
   `M† σ_k M = Σ_l R_kl σ_l` against the classical `SO(3)` matrix, per axis.
2. **Composition** for non-parallel rotation axes: the product of two
   bookkeeping matrices is the bookkeeping matrix of the composed classical
   rotation, with the Rodrigues relations and non-commutativity checked on
   the spinor and the classical side alike.
3. **Both sign branches** of the residual square-root freedom.
4. **The 720° property** of the composite.
5. **Numerical sample** over 200 random angle pairs (Rodrigues axis and
   angle written out explicitly).
6. **Boost sector** — unit determinant, branch inversion, the chiral
   factorisation, and the density and invariance laws.
7. **Which pair the boost mixes** — that the transformed pair is not the two
   spin projections, via the commutator pattern, and the four-dimensional
   carrier.
8. **The free Dirac equation** in momentum space, including the mass shell
   and the Clifford relations in 3+1 dimensions.
9. **Subsidiary claims** — that the azimuthal relative phase is forced to `φ`
   and its split to `∓φ/2`, and the spin-`j` binomial law.

Sections 7 and 8 reach beyond what the manuscript claims; the comments there
refer to a longer unpublished draft in which the same bookkeeping is carried
further.

### Running it

Python 3 and the standard library only — `fractions`, `cmath`, `math`,
`random`, `sys`. There are no dependencies and nothing to install.

```
python composition_check.py
```

The script prints one line per check and exits `0` if every check passes.
Current status: **63 PASS, 0 FAIL** (Python 3.13).
