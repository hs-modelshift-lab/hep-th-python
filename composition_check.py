# -*- coding: utf-8 -*-
"""
composition_check.py -- consistency checks for the spinor transformations of a
probability-bookkeeping construction

The construction takes two axioms about a two-outcome analyser -- that the
difference of the two probabilities transforms as a projection,
P_up(alpha) - P_dn(alpha) = cos alpha, and that the two sum to one -- and
obtains the half-angle law, the Pauli generators, the boost sector with its
indefinite norm, and the free Dirac equation in momentum space.  This script
checks the algebra of that construction: the composition of non-parallel axes,
the boost sector, and the reduction to the Dirac equation in momentum space.
The identities it evaluates are elementary; this is a check on the algebra,
not on a claim, and nothing depends on it.

EXACT, NOT NUMERICAL
--------------------
Coefficients are Gaussian rationals (a fractions.Fraction for the real and one
for the imaginary part) and expressions are polynomials in six half-angle
symbols,

    ca = cos(alpha/2)   sa = sin(alpha/2)    polar rotation, angle alpha
    cb = cos(phi/2)     sb = sin(phi/2)      azimuthal rotation, angle phi
    ch = cosh(theta/2)  sh = sinh(theta/2)   boost, rapidity theta

reduced by  sa^2 -> 1-ca^2,  sb^2 -> 1-cb^2,  sh^2 -> ch^2-1.  A structural
check is therefore a polynomial identity or it is nothing; no tolerance is
involved.  Two sections sample numerically in addition, where the quantity
checked is not polynomial in these symbols (the Rodrigues axis and angle, and
the transverse means as functions of a free phase).

THE BUILDING BLOCKS
-------------------
The two-component bookkeeping laws of the construction, transcribed as they
stand; nothing further is assumed.

  polar rotation about y, on the z doublet
      A_up(a) = A_up(0) cos(a/2) -/+ A_dn(0) sin(a/2)
      A_dn(a) = +/- A_up(0) sin(a/2) + A_dn(0) cos(a/2)

  azimuthal rotation about z
      the same half-angle law, run with the analyzer pointed into the x-y
      plane, and carried to the z basis by the change of basis.  The even
      split -+phi/2 is derived here, not postulated.

  boost along z
      A_up' =    A_up cosh(t/2) + i A_dn sinh(t/2)
      A_dn' = -i A_up sinh(t/2) +   A_dn cosh(t/2)

The residual sign of the square root is carried as a branch label
b = +1 / -1 and every structural check is run on both branches.

Axiom 1 is tested as classical vector transport,

    M^dagger sigma_k M = sum_l R_kl sigma_l ,

against the classical SO(3) matrix R.  Axiom 2 appears as unitarity in the
rotation sector and as det = 1 throughout.

Python 3, standard library only.  Exits 0 if every check passes.
"""
from fractions import Fraction
import cmath
import math
import random
import sys

# --------------------------------------------------------- Gaussian rationals
class GC(object):
    """Gaussian rational re + i*im, both exact Fractions."""
    __slots__ = ('re', 'im')

    def __init__(self, re=0, im=0):
        self.re = Fraction(re)
        self.im = Fraction(im)

    def __add__(self, o):
        return GC(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return GC(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        return GC(self.re * o.re - self.im * o.im,
                  self.re * o.im + self.im * o.re)

    def __neg__(self):
        return GC(-self.re, -self.im)

    def conj(self):
        return GC(self.re, -self.im)

    def iszero(self):
        return self.re == 0 and self.im == 0

    def tocomplex(self):
        return complex(float(self.re), float(self.im))

    def __repr__(self):
        if self.im == 0:
            return str(self.re)
        sign = '+' if self.im >= 0 else '-'
        return '(%s%s%si)' % (self.re, sign, abs(self.im))

ONE = GC(1)
I = GC(0, 1)
HALF = GC(Fraction(1, 2))
MINUS_I = GC(0, -1)
HALF_I = GC(0, Fraction(1, 2))

# ------------------------------------------------------------------ polynomials
VARNAMES = ['ca', 'sa', 'cb', 'sb', 'ch', 'sh']
NVAR = len(VARNAMES)

_m = lambda **kw: tuple(kw.get(v, 0) for v in VARNAMES)

# square-reduction rules, keyed by variable index
RULES = {
    1: [(_m(), ONE), (_m(ca=2), -ONE)],          # sa^2 = 1 - ca^2
    3: [(_m(), ONE), (_m(cb=2), -ONE)],          # sb^2 = 1 - cb^2
    5: [(_m(ch=2), ONE), (_m(), -ONE)],          # sh^2 = ch^2 - 1
}

def _reduce(mono, coef):
    """Reduce a monomial to normal form (degree in sa, sb, sh at most 1)."""
    stack = [(mono, coef)]
    out = {}
    while stack:
        m, c = stack.pop()
        for v, repl in RULES.items():
            if m[v] >= 2:
                rest = list(m)
                rest[v] -= 2
                rest = tuple(rest)
                for rm, rc in repl:
                    stack.append((tuple(rest[i] + rm[i] for i in range(NVAR)),
                                  c * rc))
                break
        else:
            acc = out.get(m)
            out[m] = c if acc is None else acc + c
    return out

class Poly(object):
    __slots__ = ('t',)

    def __init__(self, terms=None):
        self.t = {}
        if terms:
            for m, c in terms.items():
                self._add(m, c)

    def _add(self, mono, coef):
        for m, c in _reduce(mono, coef).items():
            acc = self.t.get(m)
            s = c if acc is None else acc + c
            if s.iszero():
                self.t.pop(m, None)
            else:
                self.t[m] = s

    @staticmethod
    def const(c):
        p = Poly()
        if not c.iszero():
            p.t[_m()] = c
        return p

    @staticmethod
    def var(name):
        p = Poly()
        p.t[_m(**{name: 1})] = ONE
        return p

    def __add__(self, o):
        p = Poly()
        p.t = dict(self.t)
        for m, c in o.t.items():
            p._add(m, c)
        return p

    def __sub__(self, o):
        return self + (-o)

    def __neg__(self):
        p = Poly()
        p.t = {m: -c for m, c in self.t.items()}
        return p

    def __mul__(self, o):
        p = Poly()
        for m1, c1 in self.t.items():
            for m2, c2 in o.t.items():
                p._add(tuple(m1[i] + m2[i] for i in range(NVAR)), c1 * c2)
        return p

    def conj(self):
        # the symbols are real, so only the coefficients conjugate
        p = Poly()
        p.t = {m: c.conj() for m, c in self.t.items()}
        return p

    def iszero(self):
        return not self.t

    def isreal(self):
        return all(c.im == 0 for c in self.t.values())

    def eval(self, vals):
        z = 0j
        for m, c in self.t.items():
            w = c.tocomplex()
            for i, e in enumerate(m):
                if e:
                    w *= vals[VARNAMES[i]] ** e
            z += w
        return z

    def __repr__(self):
        if not self.t:
            return '0'
        parts = []
        for m in sorted(self.t):
            mono = '*'.join('%s^%d' % (VARNAMES[i], e) if e > 1 else VARNAMES[i]
                            for i, e in enumerate(m) if e)
            parts.append('%s%s' % (self.t[m], ('*' + mono) if mono else ''))
        return ' + '.join(parts)

P0, P1, Pi = Poly.const(GC(0)), Poly.const(ONE), Poly.const(I)
CA, SA = Poly.var('ca'), Poly.var('sa')
CB, SB = Poly.var('cb'), Poly.var('sb')
CH, SH = Poly.var('ch'), Poly.var('sh')
TWO = Poly.const(GC(2))

# --------------------------------------------------------------------- matrices
def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum((A[i][x] * B[x][j] for x in range(k)), P0)
             for j in range(m)] for i in range(n)]

def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def msub(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def mscal(c, A):
    return [[Poly.const(c) * A[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]

def mscalp(p, A):
    """Scale a matrix by a polynomial."""
    return [[p * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def dagger(A):
    return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]

def det2(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]

def iszero_m(A):
    return all(A[i][j].iszero() for i in range(len(A)) for j in range(len(A[0])))

def meval(A, vals):
    return [[A[i][j].eval(vals) for j in range(len(A[0]))] for i in range(len(A))]

def zeros(n):
    return [[P0] * n for _ in range(n)]

def kron(A, B):
    n = len(A) * len(B)
    K = zeros(n)
    for i1 in range(len(A)):
        for j1 in range(len(A)):
            for i2 in range(len(B)):
                for j2 in range(len(B)):
                    K[2 * i1 + i2][2 * j1 + j2] = A[i1][j1] * B[i2][j2]
    return K

ID2 = [[P1, P0], [P0, P1]]
SIGX = [[P0, P1], [P1, P0]]
SIGY = [[P0, Poly.const(MINUS_I)], [Pi, P0]]
SIGZ = [[P1, P0], [P0, -P1]]
SIG = [SIGX, SIGY, SIGZ]
ID4 = kron(ID2, ID2)

# ------------------------------------------------------- bookkeeping generators
def M_polar(b):
    """Polar rotation about y on the z doublet, branch b."""
    s = SA if b > 0 else -SA
    return [[CA, -s], [s, CA]]

def M_inplane(b):
    """The same half-angle law with the analyzer pointed into the x-y plane
    branch b.  This is the polar construction run word for word with phi in
    the place of alpha."""
    s = SB if b > 0 else -SB
    return [[CB, -s], [s, CB]]

def basis_change(b):
    """Change of basis between the z doublet and the transverse doublet: a
    state up along z has transverse amplitudes (1, b*i)/sqrt2.  Returns
    (T, Tinv); the sqrt2 factors cancel in the conjugation."""
    bi = Pi if b > 0 else Poly.const(MINUS_I)
    T = [[P1, bi], [bi, P1]]                       # times 1/sqrt2
    Tinv = mscal(HALF, [[P1, -bi], [-bi, P1]])     # times sqrt2/2
    return T, Tinv

def M_azimuthal(b):
    """Azimuthal rotation about z on the z doublet, obtained from the in-plane
    law by the change of basis rather than postulated.  The even
    split -+phi/2 is what comes out."""
    T, Tinv = basis_change(b)
    return mmul(Tinv, mmul(M_inplane(b), T))

def M_boost(b):
    """Boost along z on the z doublet; b = +1 as written above,
    b = -1 with the placement of the imaginary unit exchanged."""
    ish = Poly.const(I if b > 0 else MINUS_I) * SH
    return [[CH, ish], [-ish, CH]]

# ------------------------------------------------------------- classical side
COSA, SINA = CA * CA - SA * SA, TWO * CA * SA
COSB, SINB = CB * CB - SB * SB, TWO * CB * SB

def R_polar_cl(sign):
    c, s = COSA, (SINA if sign > 0 else -SINA)
    return [[c, P0, s], [P0, P1, P0], [-s, P0, c]]

def R_azimuthal_cl(sign):
    c, s = COSB, (SINB if sign > 0 else -SINB)
    return [[c, -s, P0], [s, c, P0], [P0, P0, P1]]

def transport_matches(M, R):
    """Axiom 1 as vector transport: M^dagger sigma_k M == sum_l R[k][l] sigma_l."""
    Md = dagger(M)
    for k in range(3):
        lhs = mmul(Md, mmul(SIG[k], M))
        rhs = [[P0, P0], [P0, P0]]
        for l in range(3):
            rhs = madd(rhs, mscalp(R[k][l], SIG[l]))
        if not iszero_m(msub(lhs, rhs)):
            return False
    return True

def orientation(M, Rfactory):
    """Return the sign s with transport(M) == Rfactory(s), or 0 if neither."""
    for s in (+1, -1):
        if transport_matches(M, Rfactory(s)):
            return s
    return 0

def quaternion(P):
    """Decompose P = w*1 - i*(x sx + y sy + z sz); return (w, x, y, z)."""
    return (Poly.const(HALF) * (P[0][0] + P[1][1]),
            Poly.const(HALF_I) * (P[0][1] + P[1][0]),
            Poly.const(HALF) * (P[1][0] - P[0][1]),
            Poly.const(HALF_I) * (P[0][0] - P[1][1]))

def axial(R):
    h = Poly.const(HALF)
    return (h * (R[2][1] - R[1][2]), h * (R[0][2] - R[2][0]),
            h * (R[1][0] - R[0][1]))

def trace3(R):
    return R[0][0] + R[1][1] + R[2][2]

# --------------------------------------------------------------- test harness
NPASS, NFAIL = [0], [0]

def check(name, ok):
    print('  [%s] %s' % ('PASS' if ok else 'FAIL', name))
    (NPASS if ok else NFAIL)[0] += 1

def section(title):
    print()
    print('=' * 74)
    print(title)
    print('=' * 74)

# ==============================================================================
section('1. GENERATORS: unit determinant, axiom 2, and axiom-1 transport')
# Every plane carries the same half-angle law; the checks run on both branches
# of the residual square-root sign.
for b, tag in ((+1, 'branch +'), (-1, 'branch -')):
    mp, ma = M_polar(b), M_azimuthal(b)
    check('polar %s: det = 1' % tag, (det2(mp) - P1).iszero())
    check('polar %s: unitary' % tag, iszero_m(msub(mmul(dagger(mp), mp), ID2)))
    check('azimuthal %s: det = 1' % tag, (det2(ma) - P1).iszero())
    check('azimuthal %s: unitary' % tag,
          iszero_m(msub(mmul(dagger(ma), ma), ID2)))
    sp, sa_ = orientation(mp, R_polar_cl), orientation(ma, R_azimuthal_cl)
    check('polar %s: transport = classical R_y (orientation %+d)' % (tag, sp),
          sp != 0)
    check('azimuthal %s: transport = classical R_z (orientation %+d)'
          % (tag, sa_), sa_ != 0)
print('  azimuthal matrix, derived by change of basis (branch +):')
print('    [[%s, %s], [%s, %s]]'
      % (M_azimuthal(+1)[0][0], M_azimuthal(+1)[0][1],
         M_azimuthal(+1)[1][0], M_azimuthal(+1)[1][1]))
print('  -- diagonal, with the split -+phi/2; not postulated.')

# ==============================================================================
section('2. COMPOSITION for non-parallel axes (the hard test)')
SP = orientation(M_polar(+1), R_polar_cl)
SA_ = orientation(M_azimuthal(+1), R_azimuthal_cl)
P = mmul(M_azimuthal(+1), M_polar(+1))
Rc = mmul(R_azimuthal_cl(SA_), R_polar_cl(SP))

check('product is unitary with det = 1',
      iszero_m(msub(mmul(dagger(P), P), ID2)) and (det2(P) - P1).iszero())
check('axiom-1 transport of the product: P^dag sigma P = (R_z R_y) sigma',
      transport_matches(P, Rc))

w, qx, qy, qz = quaternion(P)
check('quaternion form: w, x, y, z all real',
      all(q.isreal() for q in (w, qx, qy, qz)))
check('w^2 + x^2 + y^2 + z^2 = 1',
      (w * w + qx * qx + qy * qy + qz * qz - P1).iszero())

cosg = Poly.const(HALF) * (trace3(Rc) - P1)          # cos(gamma) from the trace
check('Rodrigues I:   2w^2 - 1 = (tr R - 1)/2   [cos gamma]',
      (TWO * w * w - P1 - cosg).iszero())
check('Rodrigues II:  2w(x,y,z) = axial(R)      [sin(gamma) * axis]',
      all((TWO * w * q - a).iszero()
          for q, a in zip((qx, qy, qz), axial(Rc))))
check('Rodrigues III: x^2+y^2+z^2 = (3 - tr R)/4 [sin^2(gamma/2)]',
      (qx * qx + qy * qy + qz * qz
       - Poly.const(GC(Fraction(1, 4))) * (Poly.const(GC(3))
                                           - trace3(Rc))).iszero())
check('the composed axis is genuinely tilted (x component nonzero)',
      not qx.iszero())

Q = mmul(M_polar(+1), M_azimuthal(+1))
check('non-commutativity: M_z M_y <> M_y M_z', not iszero_m(msub(P, Q)))
check('the mirrored order composes too: transport of M_y M_z = R_y R_z',
      transport_matches(Q, mmul(R_polar_cl(SP), R_azimuthal_cl(SA_))))
check('the classical side is non-commutative in the same way',
      not iszero_m(msub(Rc, mmul(R_polar_cl(SP), R_azimuthal_cl(SA_)))))

# both sign branches propagate through products (requirement (ii))
for bz, bp in ((-1, -1), (+1, -1), (-1, +1)):
    Pb = mmul(M_azimuthal(bz), M_polar(bp))
    found = None
    for sz in (+1, -1):
        for sy in (+1, -1):
            if transport_matches(Pb, mmul(R_azimuthal_cl(sz), R_polar_cl(sy))):
                found = (sz, sy)
    check('branches (azimuthal %+d, polar %+d): composition exact, '
          'orientations (%s, %s)'
          % (bz, bp, found and '%+d' % found[0], found and '%+d' % found[1]),
          found is not None)

# ==============================================================================
section('3. THE 720 DEGREE PROPERTY of the composite')
v360 = {'ca': -1.0, 'sa': 0.0, 'cb': 1.0, 'sb': 0.0, 'ch': 1.0, 'sh': 0.0}
Pn, Rn = meval(P, v360), meval(Rc, v360)
check('at alpha = 360 deg the spinor is -1 while the classical R is +1',
      max(abs(Pn[i][j] - (-1.0 if i == j else 0))
          for i in range(2) for j in range(2)) < 1e-12
      and max(abs(Rn[i][j] - (1.0 if i == j else 0))
              for i in range(3) for j in range(3)) < 1e-12)

# ==============================================================================
section('4. NUMERICAL SAMPLE: Rodrigues axis and angle written out')
random.seed(20260714)
ok = True
for _ in range(200):
    a, p = random.uniform(-3, 3), random.uniform(-3, 3)
    vals = {'ca': math.cos(a / 2), 'sa': math.sin(a / 2),
            'cb': math.cos(p / 2), 'sb': math.sin(p / 2),
            'ch': 1.0, 'sh': 0.0}
    Pn = meval(P, vals)
    Rn = [[x.real for x in row] for row in meval(Rc, vals)]
    tr = Rn[0][0] + Rn[1][1] + Rn[2][2]
    g = math.acos(max(-1.0, min(1.0, (tr - 1) / 2)))
    if abs(math.sin(g)) > 1e-9:
        n = [(Rn[2][1] - Rn[1][2]) / (2 * math.sin(g)),
             (Rn[0][2] - Rn[2][0]) / (2 * math.sin(g)),
             (Rn[1][0] - Rn[0][1]) / (2 * math.sin(g))]
        cg, sg = math.cos(g / 2), math.sin(g / 2)
        B = [[cg - 1j * sg * n[2], -1j * sg * n[0] - sg * n[1]],
             [-1j * sg * n[0] + sg * n[1], cg + 1j * sg * n[2]]]
        d1 = max(abs(Pn[i][j] - B[i][j]) for i in range(2) for j in range(2))
        d2 = max(abs(Pn[i][j] + B[i][j]) for i in range(2) for j in range(2))
        if min(d1, d2) > 1e-10:
            ok = False
check('200 random angle pairs: the product is +-[cos(g/2) - i sin(g/2) n.sigma] '
      'of the classically composed rotation', ok)

# ==============================================================================
section('5. BOOST SECTOR: determinant, branches, chirality, the two sums')
B, Bm = M_boost(+1), M_boost(-1)
cosh_t, sinh_t = CH * CH + SH * SH, TWO * CH * SH

check('det = 1, now by cosh^2 - sinh^2 = 1', (det2(B) - P1).iszero())
check('the boost matrix is hermitian, not unitary', iszero_m(msub(dagger(B), B)))
check('the two branches boost inversely: B(-) = B(+)^(-1)',
      iszero_m(msub(mmul(B, Bm), ID2)))
check('the boost generator sits on the sigma_y axis: B = ch*1 - sh*sigma_y',
      iszero_m(msub(B, msub([[CH, P0], [P0, CH]], mscalp(SH, SIGY)))))

# chiral combinations: A_R = A_up + i A_dn scales as e^{+theta/2}
row_R = [P1 * B[0][0] + Pi * B[1][0], P1 * B[0][1] + Pi * B[1][1]]
check('A_R = A_up + i A_dn scales with e^{+theta/2} = ch + sh',
      (row_R[0] - (CH + SH)).iszero()
      and (row_R[1] - Pi * (CH + SH)).iszero())

# the unsigned sum is a density, the signed sum is invariant
BdB = mmul(dagger(B), B)
check('density law: tr(B^dag B) = 2 cosh(theta)',
      (BdB[0][0] + BdB[1][1] - TWO * cosh_t).iszero())
BT = [[B[j][i] for j in range(2)] for i in range(2)]
check('invariance law: B^T B = 1, so the plain squares are boost invariant',
      iszero_m(msub(mmul(BT, B), ID2)))

# covariance: conjugating by a rotation gives the boost about the rotated axis
Mz1 = M_azimuthal(+1)
target = msub([[CH, P0], [P0, CH]],
              madd(mscalp(COSB * SH, SIGY),
                   mscalp((SINB if SA_ < 0 else -SINB) * SH, SIGX)))
check('covariance: M_z B M_z^(-1) = ch*1 - sh*(rotated generator)',
      iszero_m(msub(mmul(Mz1, mmul(B, dagger(Mz1))), target)))

# ==============================================================================
section('6. WHICH PAIR THE BOOST MIXES: not the spin projections')
# A boost along z must commute with rotations about z and must not commute
# with polar rotations.  On the spin doublet the matrix does the opposite,
# which is the commutator argument in explicit form.
c_polar = msub(mmul(M_polar(+1), B), mmul(B, M_polar(+1)))
c_azim = msub(mmul(Mz1, B), mmul(B, Mz1))
check('read as spin: [M_z, B] <> 0, though for a z boost it would have to '
      'vanish', not iszero_m(c_azim))
check('read as spin: [M_y, B] = 0, though for a z boost it must not vanish',
      iszero_m(c_polar))

# diagonalised, the same matrix is the correct z boost on the chiral pair
E = [[P1, P1], [Poly.const(MINUS_I), Pi]]              # columns (1,-i), (1,i)
Einv = [[Poly.const(HALF), Poly.const(HALF_I)],
        [Poly.const(HALF), Poly.const(GC(0, -Fraction(1, 2)))]]
Dz = [[CH + SH, P0], [P0, CH - SH]]                    # diag(e^{t/2}, e^{-t/2})
check('B = E diag(e^{t/2}, e^{-t/2}) E^(-1), E the chiral combinations',
      iszero_m(msub(B, mmul(E, mmul(Dz, Einv)))))
check('the diagonal form has the right pattern: [M_z, D] = 0',
      iszero_m(msub(mmul(Mz1, Dz), mmul(Dz, Mz1))))
check('the diagonal form has the right pattern: [M_y, D] <> 0',
      not iszero_m(msub(mmul(M_polar(+1), Dz), mmul(Dz, M_polar(+1)))))

# the four-dimensional carrier: the mixed pair tensored with the spin doublet
TAU = [[P0, Pi], [Poly.const(MINUS_I), P0]]            # acts on (large, small)
B4 = madd(mscalp(CH, ID4), mscalp(SH, kron(TAU, SIGZ)))
R4y, R4z = kron(ID2, M_polar(+1)), kron(ID2, Mz1)
check('four-dimensional carrier: [R_z, B4] = 0',
      iszero_m(msub(mmul(R4z, B4), mmul(B4, R4z))))
check('four-dimensional carrier: [R_y, B4] <> 0',
      not iszero_m(msub(mmul(R4y, B4), mmul(B4, R4y))))
nsig_a = [[COSA, SINA], [SINA, -COSA]]                 # n.sigma, n = R_y(a) z
check('Wigner covariance: R_y B4 R_y^(-1) is the boost of equal rapidity '
      'about R_y(a) z',
      iszero_m(msub(mmul(R4y, mmul(B4, dagger(R4y))),
                    madd(mscalp(CH, ID4), mscalp(SH, kron(TAU, nsig_a))))))
# restricting the four-dimensional boost to one spin slot returns the 2x2 law
sub_up = [[B4[0][0], B4[0][2]], [B4[2][0], B4[2][2]]]
sub_dn = [[B4[1][1], B4[1][3]], [B4[3][1], B4[3][3]]]
check('restricted to spin up, B4 is the boost matrix above',
      iszero_m(msub(sub_up, B)))
check('restricted to spin down, B4 is its inverse (the other branch)',
      iszero_m(msub(sub_dn, Bm)))
check('helicity conserved: B4 mixes no spin slots',
      all(B4[i][j].iszero() for i in range(4) for j in range(4)
          if (i % 2) != (j % 2)))
print('  The partner of the boost is the large/small pair at fixed spin;')
print('  the spin projection along the boost axis is untouched.')

# ==============================================================================
section('7. THE FREE DIRAC EQUATION in momentum space')
# Ingredients, all obtained above: two branches boosting
# inversely, equal in the rest frame, and the rest frequency as the mass.
# The boost axis n = R_y(a) z = (sin a, 0, cos a) is genuinely tilted.
NSIG = [[COSA, SINA], [SINA, -COSA]]
BP = madd([[CH, P0], [P0, CH]], mscalp(SH, NSIG))      # exp(+(t/2) n.sigma)
BM = msub([[CH, P0], [P0, CH]], mscalp(SH, NSIG))      # exp(-(t/2) n.sigma)
EN, PN = cosh_t, sinh_t                                # E/m and |p|/m
EmP = msub([[EN, P0], [P0, EN]], mscalp(PN, NSIG))     # (E - sigma.p)/m
EpP = madd([[EN, P0], [P0, EN]], mscalp(PN, NSIG))     # (E + sigma.p)/m

check('the axis is normalised: (n.sigma)^2 = 1',
      iszero_m(msub(mmul(NSIG, NSIG), ID2)))
check('branch inversion on an arbitrary axis: B(-) B(+) = 1',
      iszero_m(msub(mmul(BM, BP), ID2)))
check('(E -+ sigma.p)/m = B(-+)^2, so momentum enters only through squared '
      'boost factors', iszero_m(msub(EmP, mmul(BM, BM)))
      and iszero_m(msub(EpP, mmul(BP, BP))))
check('DIRAC I:  (E - sigma.p) psi_R = m psi_L',
      iszero_m(msub(mmul(EmP, BP), BM)))
check('DIRAC II: (E + sigma.p) psi_L = m psi_R',
      iszero_m(msub(mmul(EpP, BM), BP)))
check('the mass shell is automatic: E^2 - p^2 = m^2',
      (EN * EN - PN * PN - P1).iszero())

def block4(A, Bb, C, D):
    M4 = zeros(4)
    for i in range(2):
        for j in range(2):
            M4[i][j], M4[i][j + 2] = A[i][j], Bb[i][j]
            M4[i + 2][j], M4[i + 2][j + 2] = C[i][j], D[i][j]
    return M4

ZERO2 = [[P0, P0], [P0, P0]]
DIRAC = block4(mscal(GC(-1), ID2), EmP, EpP, mscal(GC(-1), ID2))
U4 = [[BM[0][0], BM[0][1]], [BM[1][0], BM[1][1]],
      [BP[0][0], BP[0][1]], [BP[1][0], BP[1][1]]]      # u = (psi_L, psi_R)
check('(gamma.p - m) u(p) = 0 for both rest spinors, in the 4x4 Weyl form',
      iszero_m(mmul(DIRAC, U4)))
G0 = block4(ZERO2, ID2, ID2, ZERO2)
GN = block4(ZERO2, NSIG, mscal(GC(-1), NSIG), ZERO2)
ID4b = block4(ID2, ZERO2, ZERO2, ID2)
check('Clifford relations in 3+1: (g0)^2 = +1, (gn)^2 = -1, {g0, gn} = 0',
      iszero_m(msub(mmul(G0, G0), ID4b))
      and iszero_m(madd(mmul(GN, GN), ID4b))
      and iszero_m(madd(mmul(G0, GN), mmul(GN, G0))))
v0 = {'ca': 1.0, 'sa': 0.0, 'cb': 1.0, 'sb': 0.0, 'ch': 1.0, 'sh': 0.0}
check('in the rest frame psi_R = psi_L, which is A_R(0) = A_L(0)',
      max(abs(meval(msub(BP, BM), v0)[i][j])
          for i in range(2) for j in range(2)) < 1e-15)
E4, E4inv = kron(E, ID2), kron(Einv, ID2)
check('the Weyl blocks are the chiral eigenbranches of the four-dimensional '
      'boost, for n = z',
      iszero_m(msub(mmul(E4inv, mmul(B4, E4)),
                    madd(mscalp(CH, ID4), mscalp(SH, kron(SIGZ, SIGZ))))))
print('  Linearity in p is automatic, since p enters only through B^2;')
print('  E^2 = p^2 + m^2 comes out as the determinant.')

# ==============================================================================
section('8. SUBSIDIARY CLAIMS: the azimuthal phase, and spin j')
# 8a) The azimuthal plane: axiom 1 forces the relative phase to phi itself,
#     and det = 1 forces the common phase to zero, hence the split -+phi/2.
def transverse_means(al, ph, k, dl):
    """<sigma_x>, <sigma_y> for (cos(a/2) e^{i(dl - k ph/2)},
                                 sin(a/2) e^{i(dl + k ph/2)})."""
    u = math.cos(al / 2) * cmath.exp(1j * (dl - k * ph / 2))
    d = math.sin(al / 2) * cmath.exp(1j * (dl + k * ph / 2))
    return 2 * (u.conjugate() * d).real, 2 * (u.conjugate() * d).imag

random.seed(20260808)
ok_k1, ok_k2, ok_det = True, False, True
for _ in range(300):
    al, ph, dl = (random.uniform(-3, 3), random.uniform(-3, 3),
                  random.uniform(-3, 3))
    ex, ey = transverse_means(al, ph, 1.0, dl)
    if (abs(ex - math.sin(al) * math.cos(ph)) > 1e-12
            or abs(ey - math.sin(al) * math.sin(ph)) > 1e-12):
        ok_k1 = False
    if abs(transverse_means(al, ph, 2.0, dl)[0]
           - math.sin(al) * math.cos(ph)) > 1e-9:
        ok_k2 = True
    # the determinant is e^{2 i delta} for any split k, not just for k = 1
    k = random.uniform(-3, 3)
    det = (cmath.exp(1j * (dl - k * ph / 2))
           * cmath.exp(1j * (dl + k * ph / 2)))
    if abs(det - cmath.exp(2j * dl)) > 1e-12:
        ok_det = False
check('azimuthal: a relative phase of phi reproduces the transverse ledgers '
      '(axiom 1)', ok_k1)
check('azimuthal: a multiple of phi does not (k = 2 as the counter-check)',
      ok_k2)
check('azimuthal: det = e^{2 i delta} regardless of phi, so det = 1 forces '
      'the even split -+phi/2', ok_det)

# 8b) Spin j as the symmetric product of 2j aligned doublets:
#     the binomial law is |d^j_{m,j}(alpha)|^2, sums to one, mean j cos(alpha).
def wigner_d(j2, m2, mp2, beta):
    """d^j_{m,mp}(beta), with doubled quantum numbers (j2 = 2j)."""
    j, m, mp = j2 / 2.0, m2 / 2.0, mp2 / 2.0
    pre = math.sqrt(math.factorial(int(j + m)) * math.factorial(int(j - m))
                    * math.factorial(int(j + mp)) * math.factorial(int(j - mp)))
    tot = 0.0
    for s in range(0, int(2 * j) + 1):
        a1, a2, a3 = int(j + mp - s), int(j - m - s), int(m - mp + s)
        if a1 < 0 or a2 < 0 or a3 < 0:
            continue
        tot += ((-1) ** (m - mp + s)
                / (math.factorial(a1) * math.factorial(s)
                   * math.factorial(a3) * math.factorial(a2))
                * math.cos(beta / 2) ** int(2 * j + mp - m - 2 * s)
                * math.sin(beta / 2) ** int(m - mp + 2 * s))
    return pre * tot

ok_bin, ok_sum, ok_mean = True, True, True
for n in range(1, 7):                        # n = 2j, so j = 1/2 ... 3
    for beta in (0.3, 0.9, 1.7, 2.6):
        tot, mean = 0.0, 0.0
        for i in range(n + 1):               # i = number of 'up' among the 2j
            m2 = 2 * i - n
            Pm = (math.comb(n, i) * math.cos(beta / 2) ** (2 * i)
                  * math.sin(beta / 2) ** (2 * (n - i)))
            tot += Pm
            mean += (m2 / 2.0) * Pm
            if abs(wigner_d(n, m2, n, beta) ** 2 - Pm) > 1e-11:
                ok_bin = False
        if abs(tot - 1.0) > 1e-11:
            ok_sum = False
        if abs(mean - (n / 2.0) * math.cos(beta)) > 1e-11:
            ok_mean = False
check('spin j: the binomial law equals |d^j_(m,j)(alpha)|^2 for 2j = 1..6',
      ok_bin)
check('spin j: the probabilities sum to one (axiom 2)', ok_sum)
check('spin j: the mean is j cos(alpha) (axiom 1)', ok_mean)

# ==============================================================================
print()
print('=' * 74)
print('RESULT: %d PASS, %d FAIL' % (NPASS[0], NFAIL[0]))
print('=' * 74)
sys.exit(0 if NFAIL[0] == 0 else 1)
