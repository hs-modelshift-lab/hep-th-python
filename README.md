# hep-th-python

Consistency-check scripts by Hergen Scheck
([ORCID 0000-0003-1806-4048](https://orcid.org/0000-0003-1806-4048)).

## `composition_check.py`

Checks the algebra of a probability-bookkeeping construction of the spinor
transformations. The construction starts from two axioms about a two-outcome
analyser — that the difference of the two probabilities transforms as a
projection, `P↑(α) - P↓(α) = cos α`, and that the two probabilities sum to one
— and obtains the half-angle law, the Pauli generators, the boost sector with
its indefinite norm, and the free Dirac equation in momentum space. The script
verifies that the matrices so obtained transport, compose and reduce as the
construction says they do.

The script is self-contained and is not referenced by any publication. The
identities it evaluates are elementary; it is a check on the algebra, not on a
claim.

### Exact, not numerical

Coefficients are Gaussian rationals — a `fractions.Fraction` for the real and
one for the imaginary part — and expressions are polynomials in six half-angle
symbols,

| symbol | meaning | |
|---|---|---|
| `ca`, `sa` | `cos(α/2)`, `sin(α/2)` | polar rotation, angle `α` |
| `cb`, `sb` | `cos(φ/2)`, `sin(φ/2)` | azimuthal rotation, angle `φ` |
| `ch`, `sh` | `cosh(θ/2)`, `sinh(θ/2)` | boost, rapidity `θ` |

reduced by `sa² → 1-ca²`, `sb² → 1-cb²`, `sh² → ch²-1`. A structural check is
therefore a polynomial identity or it is nothing — no tolerance is involved.
Two sections sample numerically in addition, where the quantity checked is not
polynomial in these symbols: the Rodrigues axis and angle, and the transverse
means as functions of a free phase.

The building blocks are the two-component bookkeeping laws of the
construction, transcribed as they stand; nothing further is assumed.
Axiom 1 is tested as classical vector transport,
`M† σ_k M = Σ_l R_kl σ_l`, against the classical `SO(3)` matrix; axiom 2
appears as unitarity in the rotation sector and as `det = 1` throughout. The
residual sign of the square root is carried as a branch label and every
structural check is run on both branches.

### What is checked

1. **Generators** — unit determinant, unitarity, and the axiom-1 transport,
   per plane and per branch. The azimuthal matrix is *derived* here by the
   change of basis, not postulated; the even split `∓φ/2` is what comes out.
2. **Composition for non-parallel axes** — the product of two bookkeeping
   matrices is the bookkeeping matrix of the composed classical rotation,
   with the Rodrigues relations, non-commutativity on the spinor and the
   classical side alike, and both sign branches propagating through products.
3. **The 720° property** of the composite: the spinor returns `-1` where the
   classical rotation returns `+1`.
4. **Numerical sample** over 200 random angle pairs, with the Rodrigues axis
   and angle written out explicitly.
5. **Boost sector** — determinant, the inverse relation between the branches,
   the chiral factorisation, and the density and invariance laws.
6. **Which pair the boost mixes** — the commutator argument showing the pair
   is not the two spin projections, the diagonalisation, the four-dimensional
   carrier and its Wigner covariance.
7. **The free Dirac equation** in momentum space, with linearity in the
   momentum, the mass shell, and the Clifford relations in 3+1 dimensions.
8. **Subsidiary claims** — that the azimuthal relative phase is forced to `φ`
   and its split to `∓φ/2`, and the spin-`j` binomial law.

### Running it

Python 3 and the standard library only — `fractions`, `cmath`, `math`,
`random`, `sys`. There are no dependencies and nothing to install.

```
python composition_check.py
```

One line per check; exits `0` if every check passes. Current status:
**63 PASS, 0 FAIL** (Python 3.13).
