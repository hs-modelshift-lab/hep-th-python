# -*- coding: utf-8 -*-
"""
Kompositions-Check des Wahrscheinlichkeits->Amplituden-Verfahrens (Modell V6).

Haertetest aus Feedback 3.1 / To-do 11 / Paper-1-Gliederung Abschnitt 7:
Ist das Produkt zweier Buchfuehrungs-Matrizen fuer nichtkommutierende
Rotationen (verschiedene Achsen) wieder die Buchfuehrungs-Matrix der
klassisch komponierten Gesamtrotation?

Bausteine, transkribiert aus model_v6.tex (nichts weiter angenommen):
 - Polare Rotation um y (z-Dublett, EQ:amplitude_alpha):
     A_up(a)   = A_up(0)*cos(a/2) -/+ A_dn(0)*sin(a/2)
     A_dn(a)   = +/-A_up(0)*sin(a/2) + A_dn(0)*cos(a/2)
   (Vorzeichen-Freiheitsgrad = Flavour, Z. 756)
 - Azimutale Rotation um z (y-Dublett, EQ:amplitude_alphaz_imaginaryunit):
     dieselbe Halbwinkel-Matrix auf (A_yup, A_ydn); Flavour-Zweig
     B+- = A_yup +- i*A_ydn nimmt Phase e^{+-i b/2} auf.
 - i-Marker-Bruecke (Z. 781): Spin-up-z hat y-Amplituden (1, +-i)/sqrt2
   -> Basiswechsel T zwischen z-Dublett und y-Dublett.
 - Boost in z-Richtung (EQ:amplitude_rotation4):
     A_up' = A_up*cosh(t/2) + i*A_dn*sinh(t/2)
     A_dn' = -i*A_up*sinh(t/2) + A_dn*cosh(t/2)

Axiom 1 (klassische Transformation der Erwartungswerte) wird als
Vektor-Transport geprueft:  M^dagger sigma_k M = sum_l R_kl sigma_l
mit der klassischen SO(3)-Matrix R. Axiom 2 = Unitaritaet/det = 1.

Alle Rechnungen exakt: Gauss-rationale Koeffizienten, Polynome in den
Halbwinkel-Symbolen ca,sa,cb,sb,ch,sh mit Reduktionsregeln
sa^2 -> 1-ca^2, sb^2 -> 1-cb^2, sh^2 -> ch^2-1.
"""
from fractions import Fraction
import cmath
import math
import random
import sys

# ------------------------------------------------------------ Gauss-rational
class GC(object):
    """Gauss-rationale Zahl re + i*im mit exakten Fractions."""
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

# ------------------------------------------------------------------- Polynom
# Variablen: 0:ca 1:sa 2:cb 3:sb 4:ch 5:sh  (Halbwinkel-Groessen, reell)
NVAR = 6
VARNAMES = ['ca', 'sa', 'cb', 'sb', 'ch', 'sh']
# Reduktionsregeln var^2 -> Polynom (als Liste (Monom, Koeffizient))
_m = lambda **kw: tuple(kw.get(v, 0) for v in VARNAMES)
RULES = {
    1: [(_m(), ONE), (_m(ca=2), -ONE)],          # sa^2 = 1 - ca^2
    3: [(_m(), ONE), (_m(cb=2), -ONE)],          # sb^2 = 1 - cb^2
    5: [(_m(ch=2), ONE), (_m(), -ONE)],          # sh^2 = ch^2 - 1
}

def _reduce(mono, coef):
    """Reduziert ein Monom auf Normalform (Grad in sa,sb,sh <= 1)."""
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
                    nm = tuple(rest[i] + rm[i] for i in range(NVAR))
                    stack.append((nm, c * rc))
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

    def conj(self):  # Symbole sind reell -> nur Koeffizienten konjugieren
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

P0 = Poly.const(GC(0))
P1 = Poly.const(ONE)
Pi = Poly.const(I)
CA, SA = Poly.var('ca'), Poly.var('sa')
CB, SB = Poly.var('cb'), Poly.var('sb')
CH, SH = Poly.var('ch'), Poly.var('sh')

# ------------------------------------------------------------------ Matrizen
def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum((A[i][x] * B[x][j] for x in range(k)), P0)
             for j in range(m)] for i in range(n)]

def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def msub(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def mscal(c, A):
    return [[Poly.const(c) * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def dagger(A):
    return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]

def det2(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]

def iszero_m(A):
    return all(A[i][j].iszero() for i in range(len(A)) for j in range(len(A[0])))

def meval(A, vals):
    return [[A[i][j].eval(vals) for j in range(len(A[0]))] for i in range(len(A))]

ID2 = [[P1, P0], [P0, P1]]
SIGX = [[P0, P1], [P1, P0]]
SIGY = [[P0, Poly.const(-I)], [Pi, P0]]
SIGZ = [[P1, P0], [P0, -P1]]
SIG = [SIGX, SIGY, SIGZ]

# ------------------------------------------------- Buchfuehrungs-Generatoren
def My(f):
    """Polare Rotation um y auf dem z-Dublett (EQ:amplitude_alpha), Flavour f."""
    s = SA if f > 0 else -SA
    return [[CA, -s], [s, CA]]

def Rhalf_y(f):
    """Azimutale Halbwinkel-Matrix auf dem y-Dublett
    (EQ:amplitude_alphaz_imaginaryunit), Flavour f."""
    s = SB if f > 0 else -SB
    return [[CB, -s], [s, CB]]

def bridge(f):
    """i-Marker-Bruecke z-Dublett -> y-Dublett (Z. 781): Spin-up-z hat
    y-Amplituden (1, f*i)/sqrt2. Rueckgabe (T, Tinv) unitaer, det 1
    (sqrt2-Faktoren kuerzen sich in der Konjugation)."""
    fi = Pi if f > 0 else Poly.const(-I)
    T = [[P1, fi], [fi, P1]]                       # *1/sqrt2
    Tinv = mscal(HALF, [[P1, -fi], [-fi, P1]])     # *sqrt2/2
    return T, Tinv

def Mz(f):
    """Azimutale Rotation um z auf dem z-Dublett: aus dem y-Basis-Gesetz
    per Bruecken-Konjugation abgeleitet (nicht postuliert)."""
    T, Tinv = bridge(f)
    return mmul(Tinv, mmul(Rhalf_y(f), T))

def Boost(f):
    """Boost in z-Richtung auf dem z-Dublett (EQ:amplitude_rotation4);
    f = +1 wie im Manuskript, f = -1 mit getauschter i-Zuweisung (Z. 850)."""
    ish = Poly.const(I if f > 0 else -I) * SH
    return [[CH, ish], [-ish, CH]]

# --------------------------------------------------- klassische SO(3)-Seite
COSA, SINA = CA * CA - SA * SA, Poly.const(GC(2)) * CA * SA
COSB, SINB = CB * CB - SB * SB, Poly.const(GC(2)) * CB * SB

def Ry_cl(sign):
    c, s = COSA, (SINA if sign > 0 else -SINA)
    return [[c, P0, s], [P0, P1, P0], [-s, P0, c]]

def Rz_cl(sign):
    c, s = COSB, (SINB if sign > 0 else -SINB)
    return [[c, -s, P0], [s, c, P0], [P0, P0, P1]]

def transport_matches(M, R):
    """Axiom-1-Transport: M^dagger sigma_k M == sum_l R[k][l] sigma_l ?"""
    Md = dagger(M)
    for k in range(3):
        lhs = mmul(Md, mmul(SIG[k], M))
        rhs = [[P0, P0], [P0, P0]]
        for l in range(3):
            rhs = madd(rhs, mscal(GC(1), [[Poly.const(GC(1)) * R[k][l] * SIG[l][i][j]
                                           for j in range(2)] for i in range(2)]))
        if not iszero_m(msub(lhs, rhs)):
            return False
    return True

def quaternion(P):
    """Zerlegung P = w*1 - i*(x sx + y sy + z sz); liefert (w,x,y,z)."""
    w = Poly.const(HALF) * (P[0][0] + P[1][1])
    x = Poly.const(GC(0, Fraction(1, 2))) * (P[0][1] + P[1][0])
    y = Poly.const(HALF) * (P[1][0] - P[0][1])
    z = Poly.const(GC(0, Fraction(1, 2))) * (P[0][0] - P[1][1])
    return w, x, y, z

def axial(R):
    h = Poly.const(HALF)
    return (h * (R[2][1] - R[1][2]), h * (R[0][2] - R[2][0]),
            h * (R[1][0] - R[0][1]))

def trace3(R):
    return R[0][0] + R[1][1] + R[2][2]

# ---------------------------------------------------------------- Testrahmen
NPASS = [0]
NFAIL = [0]

def check(name, ok):
    print('  [%s] %s' % ('PASS' if ok else 'FAIL', name))
    (NPASS if ok else NFAIL)[0] += 1

# =========================================================================
print('=' * 74)
print('1. GENERATOREN (Axiom 2 und Axiom-1-Transport pro Achse)')
print('=' * 74)
for f, tag in ((+1, 'Flavour +'), (-1, 'Flavour -')):
    m_y, m_z = My(f), Mz(f)
    check('M_y %s: det = 1' % tag, (det2(m_y) - P1).iszero())
    check('M_y %s: unitaer' % tag, iszero_m(msub(mmul(dagger(m_y), m_y), ID2)))
    check('M_z %s: det = 1' % tag, (det2(m_z) - P1).iszero())
    check('M_z %s: unitaer' % tag, iszero_m(msub(mmul(dagger(m_z), m_z), ID2)))
    sy = +1 if transport_matches(m_y, Ry_cl(+1)) else (
        -1 if transport_matches(m_y, Ry_cl(-1)) else 0)
    sz = +1 if transport_matches(m_z, Rz_cl(+1)) else (
        -1 if transport_matches(m_z, Rz_cl(-1)) else 0)
    check('M_y %s: Transport = klassisches R_y (Orientierung %+d)' % (tag, sy),
          sy != 0)
    check('M_z %s: Transport = klassisches R_z (Orientierung %+d)' % (tag, sz),
          sz != 0)
print('  abgeleitete M_z (Flavour +): [[%s, %s], [%s, %s]]'
      % (Mz(+1)[0][0], Mz(+1)[0][1], Mz(+1)[1][0], Mz(+1)[1][1]))

print()
print('=' * 74)
print('2. HAUPTSATZ: M_z(b)*M_y(a) = Buchfuehrungs-Matrix der komponierten')
print('   Rotation R_z(b)*R_y(a)  (nichtkommutierende Achsen, exakt)')
print('=' * 74)
# Orientierungen aus Sektion 1 uebernehmen (Flavour +)
SY = +1 if transport_matches(My(+1), Ry_cl(+1)) else -1
SZ = +1 if transport_matches(Mz(+1), Rz_cl(+1)) else -1
P = mmul(Mz(+1), My(+1))
Rc = mmul(Rz_cl(SZ), Ry_cl(SY))

check('P unitaer, det = 1',
      iszero_m(msub(mmul(dagger(P), P), ID2)) and (det2(P) - P1).iszero())
check('Axiom-1-Transport der Komposition: P^dag sigma P = (R_z R_y) sigma',
      transport_matches(P, Rc))

w, qx, qy, qz = quaternion(P)
check('Quaternion-Form: w, x, y, z alle reell',
      all(q.isreal() for q in (w, qx, qy, qz)))
check('w^2 + x^2 + y^2 + z^2 = 1',
      (w * w + qx * qx + qy * qy + qz * qz - P1).iszero())

two = Poly.const(GC(2))
half3 = Poly.const(HALF)
cosg = half3 * (trace3(Rc) - P1)                      # cos(gamma) aus Spur
check('Rodrigues I:  2w^2 - 1 = (tr R - 1)/2  [cos gamma]',
      (two * w * w - P1 - cosg).iszero())
ax = axial(Rc)
check('Rodrigues II: 2w*(x,y,z) = axial(R)  [sin(gamma) * Achse]',
      all((two * w * q - a).iszero() for q, a in zip((qx, qy, qz), ax)))
check('Rodrigues III: x^2+y^2+z^2 = (3 - tr R)/4  [sin^2(gamma/2)]',
      (qx * qx + qy * qy + qz * qz
       - Poly.const(GC(Fraction(1, 4))) * (Poly.const(GC(3)) - trace3(Rc))).iszero())
check('komponierte Achse echt verkippt (x-Komponente <> 0)', not qx.iszero())

Q = mmul(My(+1), Mz(+1))
check('Nichtkommutativitaet: M_z M_y <> M_y M_z', not iszero_m(msub(P, Q)))
check('Spiegelordnung exakt: Transport von M_y M_z = R_y R_z',
      transport_matches(Q, mmul(Ry_cl(SY), Rz_cl(SZ))))
check('klassische Seite ebenso nichtkommutativ: R_z R_y <> R_y R_z',
      not iszero_m(msub(Rc, mmul(Ry_cl(SY), Rz_cl(SZ)))))

print()
print('=' * 74)
print('3. FLAVOUR-ZWEIGE (Vorzeichen-Freiheitsgrad pro Rotationsebene)')
print('=' * 74)
for fy, fz in ((-1, -1), (+1, -1), (-1, +1)):
    Pf = mmul(Mz(fz), My(fy))
    found = None
    for sz in (+1, -1):
        for sy in (+1, -1):
            if transport_matches(Pf, mmul(Rz_cl(sz), Ry_cl(sy))):
                found = (sz, sy)
    check('M_z(f=%+d) M_y(f=%+d): Komposition exakt, Orientierungen (z:%s, y:%s)'
          % (fz, fy, found and '%+d' % found[0], found and '%+d' % found[1]),
          found is not None)

print()
print('=' * 74)
print('4. DOPPELUEBERDECKUNG (720-Grad-Eigenschaft der Komposition)')
print('=' * 74)
v360 = {'ca': -1.0, 'sa': 0.0, 'cb': 1.0, 'sb': 0.0, 'ch': 1.0, 'sh': 0.0}
Pn = meval(P, v360)
Rn = meval(Rc, v360)
check('alpha = 360 Grad: P = -1 (Spinor), klassisches R = +1',
      max(abs(Pn[i][j] - (-1.0 if i == j else 0)) for i in range(2)
          for j in range(2)) < 1e-12
      and max(abs(Rn[i][j] - (1.0 if i == j else 0)) for i in range(3)
              for j in range(3)) < 1e-12)

print()
print('=' * 74)
print('5. NUMERISCHER STICHPROBEN-CHECK (Rodrigues-Achse/Winkel explizit)')
print('=' * 74)
random.seed(20260714)
ok = True
for _ in range(200):
    a, b = random.uniform(-3, 3), random.uniform(-3, 3)
    vals = {'ca': math.cos(a / 2), 'sa': math.sin(a / 2),
            'cb': math.cos(b / 2), 'sb': math.sin(b / 2), 'ch': 1.0, 'sh': 0.0}
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
check('200 Zufallswinkel: P = +-[cos(g/2) - i sin(g/2) n.sigma] der '
      'klassisch komponierten Rotation', ok)

print()
print('=' * 74)
print('6. BOOST-SEKTOR (Ausblick ueber den Rotations-Haertetest hinaus)')
print('=' * 74)
B = Boost(+1)
Bm = Boost(-1)
check('det B = 1', (det2(B) - P1).iszero())
check('B hermitesch (nicht unitaer): B^dag = B', iszero_m(msub(dagger(B), B)))
check('Zweig-Inversion (Feedback 3.2): B(-) = B(+)^(-1)',
      iszero_m(msub(mmul(B, Bm), ID2)))
# chirale Faktorisierung (EQ:amplitude_chiral): (1,+-i)-Kombination -> e^{+-t/2}
row_R = [P1 * B[0][0] + Pi * B[1][0], P1 * B[0][1] + Pi * B[1][1]]
check('A_R = A_up + i A_dn skaliert mit e^{+t/2} = ch + sh',
      (row_R[0] - (CH + SH)).iszero() and (row_R[1] - Pi * (CH + SH)).iszero())
# ungesignete Summe (Dichte): B^dag B = cosh(t)*1 + sinh(t)*(Boost-Generator)
BdB = mmul(dagger(B), B)
cosh_t, sinh_t = CH * CH + SH * SH, two * CH * SH
gen = mscal(HALF, msub(BdB, [[cosh_t, P0], [P0, cosh_t]]))
check('Dichte-Gesetz: tr(B^dag B) = 2 cosh(t)  [EQ:boost_density]',
      (BdB[0][0] + BdB[1][1] - two * cosh_t).iszero())
# gesignete Summe: psi^T psi invariant  <=>  B^T B = 1
BT = [[B[j][i] for j in range(2)] for i in range(2)]
check('Invarianz-Gesetz: B^T B = 1 (plain squares) [EQ:boost_invariant_sum]',
      iszero_m(msub(mmul(BT, B), ID2)))
# Kommutator-Muster mit den Rotations-Generatoren
c_my = msub(mmul(My(+1), B), mmul(B, My(+1)))
c_mz = msub(mmul(Mz(+1), B), mmul(B, Mz(+1)))
print('  Befund: [M_y, B] %s, [M_z, B] %s'
      % ('= 0' if iszero_m(c_my) else '<> 0',
         '= 0' if iszero_m(c_mz) else '<> 0'))
check('Boost-Generator liegt auf der sigma_y-Achse: B = ch*1 - sh*sigma_y',
      iszero_m(msub(B, msub([[CH, P0], [P0, CH]],
                            [[P0 * SH, Poly.const(-I) * SH],
                             [Pi * SH, P0 * SH]]))))
# Kovarianz: M_z(b) B M_z(b)^(-1) = Boost gleicher Rapiditaet um gedrehte Achse
Mzi = dagger(Mz(+1))  # unitaer -> Inverse = Adjungierte
BB = mmul(Mz(+1), mmul(B, Mzi))
target = msub([[CH, P0], [P0, CH]],
              madd([[COSB * SH * SIGY[i][j] for j in range(2)] for i in range(2)],
                   [[(SINB if SZ < 0 else -SINB) * SH * SIGX[i][j]
                     for j in range(2)] for i in range(2)]))
check('Kovarianz: M_z B M_z^(-1) = ch*1 - sh*(R_z-gedrehter Generator)',
      iszero_m(msub(BB, target)))

print()
print('=' * 74)
print('7. AUFLOESUNG DES ACHSEN-WOERTERBUCHS (Boost-Slots = m-Paar, nicht Spin)')
print('=' * 74)
# 7a) Widerlegung der Spin-Dublett-Lesart: ein z-Boost muss mit azimutalen
#     Rotationen um z kommutieren und darf mit polaren nicht kommutieren.
#     Die Manuskript-Matrix auf dem Spin-Dublett tut exakt das Gegenteil.
check('Widerlegung: [M_z, B] <> 0 (fuer einen z-Boost auf dem Spin-Dublett '
      'muesste er verschwinden)', not iszero_m(c_mz))
check('Widerlegung: [M_y, B] = 0 (fuer einen z-Boost auf dem Spin-Dublett '
      'duerfte er nicht verschwinden)', iszero_m(c_my))
# 7b) Diagonalisierung: B ist der korrekte z-Boost, geschrieben in der
#     i-Marker-Kombinationsbasis (1, -+i) = chirale Amplituden A_R/A_L.
E = [[P1, P1], [Poly.const(-I), Pi]]                      # Spalten (1,-i),(1,i)
Einv = [[Poly.const(HALF), Poly.const(GC(0, Fraction(1, 2)))],
        [Poly.const(HALF), Poly.const(GC(0, -Fraction(1, 2)))]]
Dz = [[CH + SH, P0], [P0, CH - SH]]                       # diag(e^{t/2}, e^{-t/2})
check('B = E diag(e^{t/2}, e^{-t/2}) E^{-1} mit E = i-Marker-Kombinationen',
      iszero_m(msub(B, mmul(E, mmul(Dz, Einv)))))
check('diagonaler z-Boost hat das korrekte Muster: [M_z, D_z] = 0',
      iszero_m(msub(mmul(Mz(+1), Dz), mmul(Dz, Mz(+1)))))
check('diagonaler z-Boost hat das korrekte Muster: [M_y, D_z] <> 0',
      not iszero_m(msub(mmul(My(+1), Dz), mmul(Dz, My(+1)))))
# 7c) Vier-dimensionaler Traeger (m-Dublett x s-Dublett):
#     K_z ~ tau (x) sigma_z mit tau = [[0,i],[-i,0]] auf dem m-Dublett.
def kron(A, Bm):
    n = len(A) * len(Bm)
    K = [[P0] * n for _ in range(n)]
    for i1 in range(len(A)):
        for j1 in range(len(A)):
            for i2 in range(len(Bm)):
                for j2 in range(len(Bm)):
                    K[2 * i1 + i2][2 * j1 + j2] = A[i1][j1] * Bm[i2][j2]
    return K

ID4 = kron(ID2, ID2)
TAU = [[P0, Pi], [Poly.const(-I), P0]]                    # (tau)^2 = 1
K4 = kron(TAU, SIGZ)
B4 = madd([[CH * ID4[i][j] for j in range(4)] for i in range(4)],
          [[SH * K4[i][j] for j in range(4)] for i in range(4)])
R4y, R4z = kron(ID2, My(+1)), kron(ID2, Mz(+1))
check('4-dim Traeger: [R_z, B4] = 0 (Boost kommutiert mit Rotation um '
      'die Boost-Achse)', iszero_m(msub(mmul(R4z, B4), mmul(B4, R4z))))
check('4-dim Traeger: [R_y, B4] <> 0 (polare Rotation dreht die Boost-Achse)',
      not iszero_m(msub(mmul(R4y, B4), mmul(B4, R4y))))
# Wigner-Kovarianz: R_y(a) B4(t) R_y(a)^{-1} = Boost um die gedrehte Achse
nsig = [[COSA, SINA], [SINA, -COSA]]                      # n.sigma, n = R_y(a) z
B4rot = madd([[CH * ID4[i][j] for j in range(4)] for i in range(4)],
             [[SH * kron(TAU, nsig)[i][j] for j in range(4)] for i in range(4)])
check('Wigner-Kovarianz: R_y B4 R_y^(-1) = Boost gleicher Rapiditaet um '
      'R_y(a)*z', iszero_m(msub(mmul(R4y, mmul(B4, dagger(R4y))), B4rot)))
# Schatten-Theoreme: beide Manuskript-Gesetze sind exakte Restriktionen.
sub_up = [[B4[0][0], B4[0][2]], [B4[2][0], B4[2][2]]]     # (a up, b up)
sub_dn = [[B4[1][1], B4[1][3]], [B4[3][1], B4[3][3]]]     # (a dn, b dn)
check('Schatten I: B4 auf dem m-Paar bei Spin up = EQ:amplitude_rotation4(_m)',
      iszero_m(msub(sub_up, B)))
check('Schatten II: B4 auf dem m-Paar bei Spin down = B^(-1) '
      '(Zweige boosten invers, Feedback 3.2)', iszero_m(msub(sub_dn, Bm)))
check('Helizitaets-Erhaltung: B4 mischt keine Spin-Slots '
      '(alle (up,down)-Bloecke = 0)',
      all(B4[i][j].iszero() for i in range(4) for j in range(4)
          if (i % 2) != (j % 2)))
print('  Lesart: Der Misch-Partner des Boosts ist das m-Dublett (Dirac:')
print('  grosse/kleine Komponente bei festem Spin); die Spin-Projektion')
print('  entlang der Boost-Achse bleibt erhalten. EQ:amplitude_rotation4')
print('  ist der korrekte z-Boost in der i-Marker-Kombinationsbasis.')

print()
print('=' * 74)
print('8. DIE FREIE DIRAC-GLEICHUNG ALS THEOREM (Impulsraum, 3+1)')
print('=' * 74)
# Zutaten (alle oben bzw. im Manuskript hergeleitet):
#  - zwei Zweige, die invers boosten (EQ:amplitude_chiral, B- = B+^{-1}),
#  - im Ruhesystem identisch (A_R(0) = A_L(0)),
#  - Ruhefrequenz = Masse (EQ:omega_m0).
# Boost-Achse n = R_y(a)*z = (sin a, 0, cos a) -- echt verkippt, exakt.
NSIG = [[COSA, SINA], [SINA, -COSA]]                     # n.sigma
BP = madd([[CH, P0], [P0, CH]], [[SH * NSIG[i][j] for j in range(2)]
                                 for i in range(2)])     # exp(+(t/2) n.sigma)
BM = msub([[CH, P0], [P0, CH]], [[SH * NSIG[i][j] for j in range(2)]
                                 for i in range(2)])     # exp(-(t/2) n.sigma)
EN, PN = cosh_t, sinh_t                                  # E/m, |p|/m
EmP = msub([[EN, P0], [P0, EN]], [[PN * NSIG[i][j] for j in range(2)]
                                  for i in range(2)])    # (E - sigma.p)/m
EpP = madd([[EN, P0], [P0, EN]], [[PN * NSIG[i][j] for j in range(2)]
                                  for i in range(2)])    # (E + sigma.p)/m
check('Achse normiert: (n.sigma)^2 = 1', iszero_m(msub(mmul(NSIG, NSIG), ID2)))
check('Zweig-Inversion auf beliebiger Achse: B(-) B(+) = 1',
      iszero_m(msub(mmul(BM, BP), ID2)))
check('(E - sigma.p)/m = B(-)^2 und (E + sigma.p)/m = B(+)^2 '
      '(Impuls betritt die Rechnung nur ueber Boost-Quadrate)',
      iszero_m(msub(EmP, mmul(BM, BM))) and iszero_m(msub(EpP, mmul(BP, BP))))
check('DIRAC I:  (E - sigma.p) psi_R = m psi_L  [als Matrixidentitaet '
      '(E-sigma.p)B(+) = m B(-)]', iszero_m(msub(mmul(EmP, BP), BM)))
check('DIRAC II: (E + sigma.p) psi_L = m psi_R  [(E+sigma.p)B(-) = m B(+)]',
      iszero_m(msub(mmul(EpP, BM), BP)))
check('Massenschale automatisch: E^2 - p^2 = m^2  [cosh^2 - sinh^2 = 1]',
      (EN * EN - PN * PN - P1).iszero())
# 4x4-Form: (gamma.p - m) u(p) = 0 in der Weyl-Darstellung, u = (psi_R, psi_L)
def block4(A, B_, C, Dm):
    M4 = [[P0] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            M4[i][j], M4[i][j + 2] = A[i][j], B_[i][j]
            M4[i + 2][j], M4[i + 2][j + 2] = C[i][j], Dm[i][j]
    return M4

DIRAC = block4(mscal(GC(-1), ID2), EmP, EpP, mscal(GC(-1), ID2))
U4 = [[BM[0][0], BM[0][1]], [BM[1][0], BM[1][1]],
      [BP[0][0], BP[0][1]], [BP[1][0], BP[1][1]]]        # u = (psi_L, psi_R)
check('(gamma.p - m) u(p) = 0 fuer beide Ruhe-Spinoren (4x4 Weyl-Form)',
      iszero_m(mmul(DIRAC, U4)))
G0 = block4([[P0, P0], [P0, P0]], ID2, ID2, [[P0, P0], [P0, P0]])
GN = block4([[P0, P0], [P0, P0]], NSIG, mscal(GC(-1), NSIG), [[P0, P0], [P0, P0]])
ID4b = block4(ID2, [[P0, P0], [P0, P0]], [[P0, P0], [P0, P0]], ID2)
check('Clifford in 3+1: (gamma0)^2 = +1, (gamma_n)^2 = -1, '
      '{gamma0, gamma_n} = 0',
      iszero_m(msub(mmul(G0, G0), ID4b))
      and iszero_m(madd(mmul(GN, GN), ID4b))
      and iszero_m(madd(mmul(G0, GN), mmul(GN, G0))))
# Ruhesystem: B(+) = B(-) = 1 -> psi_R = psi_L (Manuskript: A_R(0) = A_L(0)),
# Dirac-Gleichung reduziert auf E = m = hbar*omega (EQ:omega_m0)
v0 = {'ca': 1.0, 'sa': 0.0, 'cb': 1.0, 'sb': 0.0, 'ch': 1.0, 'sh': 0.0}
check('Ruhesystem: psi_R = psi_L (Manuskript-Bedingung A_R(0) = A_L(0)), '
      'Gleichung reduziert auf E = m',
      max(abs(meval(msub(BP, BM), v0)[i][j]) for i in range(2)
          for j in range(2)) < 1e-15)
# Anschluss an Sektion 7: die Weyl-Bloecke B(+-) sind exakt die i-Marker-
# Eigenzweige des m(x)s-Boosts (Nachtrag 33), hier fuer n = z
E4 = kron(E, ID2)
E4inv = kron(Einv, ID2)
target8 = madd([[CH * ID4[i][j] for j in range(4)] for i in range(4)],
               [[SH * kron(SIGZ, SIGZ)[i][j] for j in range(4)]
                for i in range(4)])
check('Weyl-Bloecke = i-Marker-Eigenzweige des 4-dim Boosts: '
      'E^(-1) B4 E = diag(B(+), B(-)) fuer n = z',
      iszero_m(msub(mmul(E4inv, mmul(B4, E4)), target8)))
print('  Lesart: Aus (i) invers boostenden Zweigen, (ii) Zweig-Gleichheit')
print('  im Ruhesystem und (iii) omega = m_0 folgt die freie Dirac-')
print('  Gleichung als Identitaet -- Linearitaet in p automatisch, da p')
print('  nur ueber B^2 eintritt; E^2 = p^2 + m^2 als Determinante.')

print()
print('=' * 74)
print('9. GUTACHTER-NACHTRAEGE (Paper 1, Revision)')
print('=' * 74)
# 9a) Azimutale Ebene: Relativphase k*phi und gemeinsame Phase delta.
#     Behauptung im Text: Axiom 1 in der x-y-Ebene erzwingt k = 1,
#     Axiom 2 (det = 1) erzwingt delta = 0 und damit die gerade Teilung.
def transverse_means(al, ph, k, dl):
    """<sigma_x>, <sigma_y> fuer (cos(a/2) e^{i(dl - k ph/2)}, sin(a/2) e^{i(dl + k ph/2)})."""
    u = math.cos(al / 2) * cmath.exp(1j * (dl - k * ph / 2))
    d = math.sin(al / 2) * cmath.exp(1j * (dl + k * ph / 2))
    ex = 2 * (u.conjugate() * d).real
    ey = 2 * (u.conjugate() * d).imag
    return ex, ey

random.seed(20260808)
ok_k1, ok_k2 = True, False
for _ in range(300):
    al, ph, dl = (random.uniform(-3, 3), random.uniform(-3, 3),
                  random.uniform(-3, 3))
    ex, ey = transverse_means(al, ph, 1.0, dl)
    if (abs(ex - math.sin(al) * math.cos(ph)) > 1e-12
            or abs(ey - math.sin(al) * math.sin(ph)) > 1e-12):
        ok_k1 = False
    ex2, _ = transverse_means(al, ph, 2.0, dl)
    if abs(ex2 - math.sin(al) * math.cos(ph)) > 1e-9:
        ok_k2 = True
check('Azimut: k = 1 reproduziert die Projektion cos(phi) in der x-y-Ebene '
      '(Axiom 1)', ok_k1)
check('Azimut: k <> 1 tut es nicht (k = 2 als Gegenprobe)', ok_k2)
ok_det = True
for _ in range(300):
    ph, dl, k = (random.uniform(-3, 3), random.uniform(-3, 3),
                 random.uniform(-3, 3))
    det = (cmath.exp(1j * (dl - k * ph / 2)) * cmath.exp(1j * (dl + k * ph / 2)))
    if abs(det - cmath.exp(2j * dl)) > 1e-12:
        ok_det = False
check('Azimut: det der Diagonalmatrix = e^{2 i delta}, unabhaengig von phi '
      '-> det = 1 erzwingt delta = 0 (gerade Teilung -+ phi/2)', ok_det)

# 9b) k > 2 Ausgaenge als symmetrisches Kompositum von 2j Zwei-Ausgangs-
#     Systemen: Binomialgesetz = |d^j_{m,j}|^2, Summe 1, <m> = j cos(alpha).
def wigner_d(j2, m2, mp2, beta):
    """d^j_{m,mp}(beta), Argumente als doppelte Quantenzahlen (j2 = 2j)."""
    j, m, mp = j2 / 2.0, m2 / 2.0, mp2 / 2.0
    pre = math.sqrt(math.factorial(int(j + m)) * math.factorial(int(j - m))
                    * math.factorial(int(j + mp)) * math.factorial(int(j - mp)))
    tot = 0.0
    for s in range(0, int(2 * j) + 1):
        a1, a2 = int(j + mp - s), int(j - m - s)
        a3 = int(m - mp + s)
        if a1 < 0 or a2 < 0 or a3 < 0:
            continue
        tot += ((-1) ** (m - mp + s)
                / (math.factorial(a1) * math.factorial(s) * math.factorial(a3)
                   * math.factorial(a2))
                * math.cos(beta / 2) ** int(2 * j + mp - m - 2 * s)
                * math.sin(beta / 2) ** int(m - mp + 2 * s))
    return pre * tot

ok_bin, ok_sum, ok_mean = True, True, True
for n in range(1, 7):                       # n = 2j, also j = 1/2 ... 3
    for beta in (0.3, 0.9, 1.7, 2.6):
        tot, mean = 0.0, 0.0
        for i in range(n + 1):              # i = Zahl der 'up' unter den 2j
            m2 = 2 * i - n                  # doppeltes m
            P = (math.comb(n, i) * math.cos(beta / 2) ** (2 * i)
                 * math.sin(beta / 2) ** (2 * (n - i)))
            tot += P
            mean += (m2 / 2.0) * P
            if abs(wigner_d(n, m2, n, beta) ** 2 - P) > 1e-11:
                ok_bin = False
        if abs(tot - 1.0) > 1e-11:
            ok_sum = False
        if abs(mean - (n / 2.0) * math.cos(beta)) > 1e-11:
            ok_mean = False
check('Spin-j: Binomialgesetz = |d^j_(m,j)(alpha)|^2 fuer 2j = 1..6 '
      '[EQ:spin_j]', ok_bin)
check('Spin-j: Summe der P_m = 1 (Axiom 2)', ok_sum)
check('Spin-j: <m> = j cos(alpha) (Axiom 1)', ok_mean)

print()
print('=' * 74)
print('ERGEBNIS: %d PASS, %d FAIL' % (NPASS[0], NFAIL[0]))
print('=' * 74)
sys.exit(0 if NFAIL[0] == 0 else 1)
