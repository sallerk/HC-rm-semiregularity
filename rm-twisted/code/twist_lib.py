"""twist_lib.py -- exact arithmetic for twisted secant classes on abelian fourfolds with real multiplication.

Model (NOTE.md Section 2).  F = Q(sqrt t), Theta = Tr E, Theta_1, Theta_2 its components at the two real places
(Theta = Theta_1 + Theta_2, L_s = sqrt t (Theta_1 - Theta_2)).  A rational class that is a polynomial in
Theta_1, Theta_2 (degree <= 2 in each) is a sigma-Hermitian 3 x 3 matrix X over F:
    x = sum_{i,j} X[i][j] Theta_1^i Theta_2^j,   X[j][i] = sigma(X[i][j]),
stored as the six entries x00, x10, x20, x11, x21, x22 (x00, x11, x22 rational).  Twisting by e^{-L_h} is
X -> E_h X sigma(E_h)^T with E_h the matrix of multiplication by e^{-h Theta} on <1, Theta, Theta^2>.
Integrality is degreewise: x00 in Z, x10 in J2, Q_{x20} + x11 T0 in H^4(Z) (the lattice L4), T_{x21} in H^6(Z)
(the lattice J6), x22 W in Z, where Q_m = m Theta_1^2 + mbar Theta_2^2, T0 = Theta_1 Theta_2,
T_k = k Theta_1^2 Theta_2 + kbar Theta_1 Theta_2^2 and W = int Theta_1^2 Theta_2^2.  The lattices are computed
from the cohomology ring of the lattice (shape_search.setup_forms: exact forms on wedge^* Hom(Lambda, Z)).

The secant space B(Theta, q) is the set of X = r C(x)C + q v T0 + g Theta_1 C_2 + gbar C_1 Theta_2 with
C_k = 1 - (q/2) Theta_k^2, r = x + z, v = z - x, g = y + w sqrt t (c = x alpha + y beta + z alpha~ +
w sqrt t beta~); P = <alpha, beta> (z = w = 0), P~ = <alpha~, beta~> (x = y = 0), the rank-one cone q r v = N(g).
"""
import importlib.util
import math
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "rm-semiregularity", "code", "shape_search.py")
_spec = importlib.util.spec_from_file_location("shape_search", os.path.normpath(SRC))
SS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SS)
flint = SS.flint


# ----------------------------------------------------------------------------- the field
class Q2:
    """u + v sqrt t, u, v rational."""
    __slots__ = ("t", "u", "v")

    def __init__(s, t, u=0, v=0):
        s.t, s.u, s.v = t, Fr(u), Fr(v)

    def _c(s, o):
        return o if isinstance(o, Q2) else Q2(s.t, o)

    def __add__(s, o):
        o = s._c(o)
        return Q2(s.t, s.u + o.u, s.v + o.v)

    __radd__ = __add__

    def __neg__(s):
        return Q2(s.t, -s.u, -s.v)

    def __sub__(s, o):
        return s + (-s._c(o))

    def __rsub__(s, o):
        return s._c(o) - s

    def __mul__(s, o):
        o = s._c(o)
        return Q2(s.t, s.u * o.u + s.t * s.v * o.v, s.u * o.v + s.v * o.u)

    __rmul__ = __mul__

    def conj(s):
        return Q2(s.t, s.u, -s.v)

    def norm(s):
        return s.u * s.u - s.t * s.v * s.v

    def tr(s):
        return 2 * s.u

    def inv(s):
        n = s.norm()
        return Q2(s.t, s.u / n, -s.v / n)

    def __truediv__(s, o):
        return s * s._c(o).inv()

    def __eq__(s, o):
        o = s._c(o)
        return s.u == o.u and s.v == o.v

    def __hash__(s):
        return hash((s.u, s.v))

    def is_rational(s):
        return s.v == 0

    def emb(s, k):
        """the real embedding k = 0 (sqrt t > 0) or 1."""
        r = math.sqrt(s.t)
        return float(s.u) + (1 if k == 0 else -1) * float(s.v) * r

    def totally_positive(s):
        return s.u > 0 and s.norm() > 0

    def __repr__(s):
        return "(%s%s%s*r%d)" % (s.u, "+" if s.v > 0 else "-", abs(s.v), s.t) if s.v else str(s.u)


def dfield(t):
    return t if t % 4 == 1 else 4 * t


# ----------------------------------------------------------------------------- Hermitian matrices
class Herm:
    """the class sum X[i][j] Theta_1^i Theta_2^j; entries x00, x11, x22 rational, x10, x20, x21 in F."""
    __slots__ = ("t", "x00", "x10", "x20", "x11", "x21", "x22")

    def __init__(s, t, x00=0, x10=0, x20=0, x11=0, x21=0, x22=0):
        s.t = t
        s.x00, s.x11, s.x22 = Fr(x00), Fr(x11), Fr(x22)
        s.x10, s.x20, s.x21 = (x if isinstance(x, Q2) else Q2(t, x) for x in (x10, x20, x21))

    def mat(s):
        t, c = s.t, lambda a: a.conj()
        z = lambda a: Q2(t, a)
        return [[z(s.x00), c(s.x10), c(s.x20)],
                [s.x10, z(s.x11), c(s.x21)],
                [s.x20, s.x21, z(s.x22)]]

    @staticmethod
    def from_mat(t, M):
        for i in range(3):
            for j in range(3):
                assert M[j][i] == M[i][j].conj(), "not sigma-Hermitian"
        for i in range(3):
            assert M[i][i].is_rational()
        return Herm(t, M[0][0].u, M[1][0], M[2][0], M[1][1].u, M[2][1], M[2][2].u)

    def __add__(s, o):
        return Herm(s.t, s.x00 + o.x00, s.x10 + o.x10, s.x20 + o.x20, s.x11 + o.x11, s.x21 + o.x21,
                    s.x22 + o.x22)

    def scale(s, a):
        a = Fr(a)
        return Herm(s.t, a * s.x00, s.x10 * a, s.x20 * a, a * s.x11, s.x21 * a, a * s.x22)

    def __neg__(s):
        return s.scale(-1)

    def key(s):
        return (s.x00, s.x10.u, s.x10.v, s.x20.u, s.x20.v, s.x11, s.x21.u, s.x21.v, s.x22)

    def is_zero(s):
        return not any(s.key())

    def __repr__(s):
        return "Herm(r=%s, x10=%s, x20=%s, x11=%s, x21=%s, x22=%s)" % (s.x00, s.x10, s.x20, s.x11, s.x21, s.x22)


def mmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(3)), Q2(A[0][0].t)) for j in range(3)] for i in range(3)]


def mconj_t(A):
    return [[A[j][i].conj() for j in range(3)] for i in range(3)]


def E_twist(t, h):
    """matrix of multiplication by e^{-h Theta} on <1, Theta, Theta^2> (columns = images)."""
    h = h if isinstance(h, Q2) else Q2(t, h)
    one, zero = Q2(t, 1), Q2(t, 0)
    return [[one, zero, zero], [-h, one, zero], [h * h * Fr(1, 2), -h, one]]


def twist(X, h):
    """e^{-L_h} x."""
    E = E_twist(X.t, h)
    return Herm.from_mat(X.t, mmul(mmul(E, X.mat()), mconj_t(E)))


def secant(t, q, x, y, z, w):
    """the class x alpha + y beta + z alpha~ + w sqrt t beta~ of B(Theta, q)."""
    q = Fr(q)
    r, v, g = Fr(x) + Fr(z), Fr(z) - Fr(x), Q2(t, y, w)
    return Herm(t, r, g, Q2(t, -q * r / 2), q * v, g.conj() * (-q / 2), r * q * q / 4)


def mvec(A, v):
    return tuple(sum((A[i][k] * v[k] for k in range(3)), Q2(A[0][0].t)) for i in range(3))


def outer(t, a, lam=1):
    """lam (a (x) sigma(a)) for a = (a0, a1, a2) in F^3 (a in <1, Theta, Theta^2>), lam rational."""
    M = [[a[i] * a[j].conj() * Fr(lam) for j in range(3)] for i in range(3)]
    return Herm.from_mat(t, M)


def chi_herm(X, W):
    """Mukai square: W [2 x00 x22 + x11^2 + 2 N(x20) - 2 Tr(x01 x21)]."""
    return W * (2 * X.x00 * X.x22 + X.x11 * X.x11 + 2 * X.x20.norm() - 2 * (X.x10.conj() * X.x21).tr())


def beta_F(a, b):
    """the F-bilinear form on <1, Theta, Theta^2>: coefficient of Theta^2 in tau(a) b."""
    return a[0] * b[2] + a[2] * b[0] - a[1] * b[1]


# ----------------------------------------------------------------------------- lattices
class Lat:
    """the cohomology of one lattice: P (products Theta^a L_s^b as forms), s0, t, W, and the degree lattices."""

    def __init__(s, t, P, s0, label=""):
        s.t, s.P, s.s0, s.label = t, P, s0, label
        T12 = SS.add((Fr(1, 4), P[(2, 0)]), (Fr(-1, 4 * t), P[(0, 2)]))
        sq = SS.wmul(T12, T12).get(255, 0)
        s.W = abs(sq)
        s.wsign = 1 if sq * s0 > 0 else -1
        # degree 2: L_g = g0 Theta + g1 L_s for g = g0 + g1 sqrt t
        s.J2 = SS.lattice2([P[(1, 0)], P[(0, 1)]])
        # degree 4: Q_{a + b sqrt t} + c T0 with Q_1 = (Theta^2 + L_s^2/t)/2, Q_sqrt t = Theta L_s, T0 = T12
        Q1 = SS.add((Fr(1, 2), P[(2, 0)]), (Fr(1, 2 * t), P[(0, 2)]))
        Qs = P[(1, 1)]
        s.L4forms = [Q1, Qs, T12]
        s.L4 = SS.lattice2(s.L4forms)
        # degree 6: T_{a + b sqrt t} with T_1 = Theta^3/3, T_sqrt t = Theta^2 L_s
        s.J6forms = [SS.add((Fr(1, 3), P[(3, 0)])), P[(2, 1)]]
        s.J6 = SS.lattice2(s.J6forms)
        s.L4inv = inv3(s.L4)
        s.J2inv, s.J6inv = inv2(s.J2), inv2(s.J6)

    def in_J2(s, g):
        return in_lat2(s.J2inv, (g.u, g.v))

    def in_J6(s, k):
        return in_lat2(s.J6inv, (k.u, k.v))

    def in_L4(s, m, n):
        u = [sum(s.L4inv[i][j] * x for j, x in enumerate((m.u, m.v, Fr(n)))) for i in range(3)]
        return all(x.denominator == 1 for x in u)

    def integral(s, X):
        return (X.x00.denominator == 1 and s.in_J2(X.x10) and s.in_L4(X.x20, X.x11) and s.in_J6(X.x21)
                and (X.x22 * s.W).denominator == 1)

    def form(s, X):
        """the class as a form (slow path, for checks)."""
        t, P = s.t, s.P
        coef = {}
        # Theta_1 = (Theta + L_s/sqrt t)/2, Theta_2 = (Theta - L_s/sqrt t)/2; lam = L_s/sqrt t
        inv_rt = Q2(t, 0, Fr(1, t))
        M = X.mat()
        for i in range(3):
            for j in range(3):
                if M[i][j] == 0:
                    continue
                # (Theta + lam)^i (Theta - lam)^j / 2^(i+j), as polynomial in Theta, lam
                poly = {(0, 0): Q2(t, 1)}
                for _ in range(i):
                    poly = pmul(t, poly, {(1, 0): Q2(t, 1), (0, 1): Q2(t, 1)})
                for _ in range(j):
                    poly = pmul(t, poly, {(1, 0): Q2(t, 1), (0, 1): Q2(t, -1)})
                for (a, b), c in poly.items():
                    cc = M[i][j] * c * Fr(1, 2 ** (i + j))
                    for _ in range(b):
                        cc = cc * inv_rt
                    coef[(a, b)] = coef.get((a, b), Q2(t, 0)) + cc
        out = {}
        for (a, b), c in coef.items():
            assert c.is_rational(), "class not rational"
            if c.u:
                for k, v in P[(a, b)].items():
                    out[k] = out.get(k, 0) + c.u * v
        return {k: v for k, v in out.items() if v != 0}

    def chi(s, X):
        return chi_herm(X, s.W)


def pmul(t, A, B):
    out = {}
    for (a1, b1), c1 in A.items():
        for (a2, b2), c2 in B.items():
            k = (a1 + a2, b1 + b2)
            out[k] = out.get(k, Q2(t, 0)) + c1 * c2
    return out


def inv2(B):
    det = B[0][0] * B[1][1] - B[0][1] * B[1][0]
    return [[B[1][1] / det, -B[0][1] / det], [-B[1][0] / det, B[0][0] / det]]


def inv3(B):
    m = flint.fmpq_mat([[flint.fmpq(x.numerator, x.denominator) for x in r] for r in B]).inv()
    return [[Fr(int(m[i, j].p), int(m[i, j].q)) for j in range(3)] for i in range(3)]


def in_lat2(Binv, z):
    return all((Binv[i][0] * z[0] + Binv[i][1] * z[1]).denominator == 1 for i in range(2))


def free_lattice(t, g1, g2, label=None):
    dl = SS.delta(t)
    P, s0 = SS.setup(t, dl * g1, dl * g1 * g2)
    return Lat(t, P, s0, label or "free(%s,%s)" % (g1, g2))


def y_lattice(t, g1, g2, A, label):
    M1, Ms, SY, P, s0, T12, W = SS.ylattice(t, g1, g2, A)
    return Lat(t, P, s0, label)


def kind(t, q, x, y, z, w):
    r, v, Ng = x + z, z - x, Fr(y) * y - t * Fr(w) * w
    if z == 0 and w == 0:
        return "P"
    if x == 0 and y == 0:
        return "P~"
    if r == 0 and y == 0 and w == 0:
        return "line"
    if Fr(q) * r * v == Ng:
        return "cone"
    return "generic"
