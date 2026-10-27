"""check_sym2.py -- exact checks of the RM-surface model of NOTE.md Section 3.

Model.  W = <w0, w1> (w0 = 1, w1 = Theta_A), V_F = Sym^2 W via 1 -> w0 w0, Theta -> w0 w1 + w1 w0,
Theta^2 -> 2 w1 w1; H = W_1 (x) W_2 with basis e_ab = w_a (x) w_b (a: place 1, b: place 2), Mukai form
b(e_ab, e_a'b') = om(a, a') om(b, b'), om(0, 1) = 1 = -om(1, 0).  A class X (Hermitian 3 x 3) goes to the symmetric
4 x 4 matrix S[(a, b), (a', b')] = sum_ij X_ij m_i[a][a'] m_j[b][b'] (entries in F, read at the place sigma_1), and
eps = S b.  Checks:
 (S1) S is symmetric and tr(eps) = 0, on random twisted secant classes c' = e^{-L_h} c (c in B(Theta, q));
 (S2) on the same classes: tr(eps^2) = 4 chi / W and det(eps) = 4 Pi / W^2, Pi = (chi/2)^2 - (W q n)^2
      (both sides scale alike in W, so W = 1 is used);
 (S3) on the same classes: det(x - eps^2) = (x - m_+)^2 (x - m_-)^2 with m_+ = 2 chi(c_P)/W, m_- = 2 chi(c_P~)/W;
 (S4) S symmetric and tr(eps) = 0 for a general Hermitian X (rank 3);
 (S5) principal free lattices (t = 5, 2, 3, 13): in the basis f = (e00, e11, e01 + e10, om e01 + om' e10) of
      H_Z = {sum c_ab e_ab : c00, c11 in Z, c01 in O_F, c10 = sigma(c01)}, det b(f_I, f_J) = d_F, and S has integral
      coordinates S(f_I^*, f_J^*) for 400 random X in L_Herm (W = 4), so det(eps) = d_F det(S_f) and
      Pi = (W^2/4) det(eps) = 4 d_F det(S_f) (exact);
 (S6) for each lattice Lambda in a list (free types (g1, g2), random lattices of index N^k over O_F^4), the quartic
      Pi(X) = (W^2/4) det(eps(X)) on L_Herm(Lambda): its content (the gcd of the finite differences Delta^alpha Pi(0),
      |alpha| <= 4, in a Z-basis of L_Herm) -- Pi is integer-valued with values in c Z exactly when c divides the
      content.  Printed: content / (4 d_F).  This is exact and covers every X in L_Herm(Lambda).
Usage: python -I -B check_sym2.py      Expected last line "ALL CHECKS PASS: True".
"""
import importlib.util
import itertools
import math
import os
import random
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("twist_lib", os.path.join(HERE, "twist_lib.py"))
TL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TL)
SS, Q2, Herm = TL.SS, TL.Q2, TL.Herm

M = [[[1, 0], [0, 0]], [[0, 1], [1, 0]], [[0, 0], [0, 2]]]
IDX = [(0, 0), (0, 1), (1, 0), (1, 1)]
OM = [[0, 1], [-1, 0]]
PERMS = [(p, (-1) ** sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]))
         for p in itertools.permutations(range(4))]


def dF(t):
    return t if t % 4 == 1 else 4 * t


def to_S(X):
    t, A = X.t, X.mat()
    S = [[Q2(t, 0) for _ in range(4)] for _ in range(4)]
    for I, (a, b) in enumerate(IDX):
        for J, (a2, b2) in enumerate(IDX):
            tot = Q2(t, 0)
            for i in range(3):
                for j in range(3):
                    c = M[i][a][a2] * M[j][b][b2]
                    if c:
                        tot = tot + A[i][j] * c
            S[I][J] = tot
    return S


def bmat(t):
    return [[Q2(t, OM[a][a2] * OM[b][b2]) for (a2, b2) in IDX] for (a, b) in IDX]


def mm(A, B):
    t = A[0][0].t
    n = len(A)
    return [[sum((A[i][k] * B[k][j] for k in range(n)), Q2(t, 0)) for j in range(n)] for i in range(n)]


def det4(A):
    t = A[0][0].t
    tot = Q2(t, 0)
    for p, sg in PERMS:
        pr = Q2(t, sg)
        for i in range(4):
            pr = pr * A[i][p[i]]
        tot = tot + pr
    return tot


def charpoly(A):
    """Faddeev-LeVerrier: coefficients c_0..c_n of det(x I - A) = sum c_k x^k."""
    t, n = A[0][0].t, len(A)
    coeffs = [Q2(t, 0)] * (n + 1)
    coeffs[n] = Q2(t, 1)
    Mk = [[Q2(t, 0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        Mk = mm(A, Mk)
        for i in range(n):
            Mk[i][i] = Mk[i][i] + coeffs[n - k + 1]
        AM = mm(A, Mk)
        coeffs[n - k] = sum((AM[i][i] for i in range(n)), Q2(t, 0)) * Fr(-1, k)
    return coeffs


def inv4(A):
    """inverse over F by Gauss-Jordan."""
    t, n = A[0][0].t, 4
    R = [list(A[i]) + [Q2(t, int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if R[r][c] != 0)
        R[c], R[p] = R[p], R[c]
        iv = R[c][c].inv()
        R[c] = [x * iv for x in R[c]]
        for r in range(n):
            if r != c and R[r][c] != 0:
                f = R[r][c]
                R[r] = [R[r][k] - f * R[c][k] for k in range(2 * n)]
    return [row[n:] for row in R]


def tr_(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A))]


def rfr(rng, den=5, num=5):
    return Fr(rng.randint(-num, num), rng.randint(1, den))


def comb(rng, B):
    k = [rng.randint(-3, 3) for _ in range(len(B[0]))]
    return [sum(B[i][j] * k[j] for j in range(len(k))) for i in range(len(B))]


def lherm_basis(L):
    """a Z-basis of L_Herm(Lambda) = Z + J2 + L4 + J6 + (1/W) Z (Hermitian matrices)."""
    t = L.t
    z = Q2(t, 0)
    out = [Herm(t, 1, z, z, 0, z, 0)]
    for j in range(2):
        out.append(Herm(t, 0, Q2(t, L.J2[0][j], L.J2[1][j]), z, 0, z, 0))
    for j in range(3):
        out.append(Herm(t, 0, z, Q2(t, L.L4[0][j], L.L4[1][j]), L.L4[2][j], z, 0))
    for j in range(2):
        out.append(Herm(t, 0, z, z, 0, Q2(t, L.J6[0][j], L.J6[1][j]), 0))
    out.append(Herm(t, 0, z, z, 0, z, Fr(1, 1) / L.W))
    return out


def frgcd(a, b):
    a, b = Fr(a), Fr(b)
    return Fr(math.gcd(a.numerator * b.denominator, b.numerator * a.denominator), a.denominator * b.denominator)


def content_Pi(L):
    """gcd of the finite differences of Pi = (W^2/4) det eps on a Z-basis of L_Herm (degree 4, 9 variables)."""
    t, W = L.t, L.W
    b = bmat(t)
    E = [mm(to_S(B), b) for B in lherm_basis(L)]
    pts = [a for a in itertools.product(range(5), repeat=9) if sum(a) <= 4]
    val = {}
    for a in pts:
        A = [[sum((E[k][i][j] * a[k] for k in range(9) if a[k]), Q2(t, 0)) for j in range(4)] for i in range(4)]
        d = det4(A)
        assert d.is_rational()
        val[a] = W * W / 4 * d.u
    g = Fr(0)
    for a in pts:
        c = Fr(0)
        for bb in itertools.product(*[range(x + 1) for x in a]):
            co = 1
            for x, y in zip(a, bb):
                co *= (-1) ** (x - y) * math.comb(x, y)
            c += co * val[bb]
        g = frgcd(g, c)
    return g


ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print("(%s) %s" % (name, bool(cond)), flush=True)


def main():
    rng = random.Random(20261025)
    g1 = g2 = g3 = g4 = True
    for t in (5, 2, 3, 13):
        b = bmat(t)
        for _ in range(25):
            q = Fr(rng.randint(1, 9), rng.randint(1, 6))
            x, y, z, w = (rfr(rng) for _ in range(4))
            h = Q2(t, rfr(rng, 4, 3), rfr(rng, 4, 3))
            X = TL.twist(TL.secant(t, q, x, y, z, w), h)
            S = to_S(X)
            eps = mm(S, b)
            g1 &= all(S[i][j] == S[j][i] for i in range(4) for j in range(4))
            g1 &= sum((eps[i][i] for i in range(4)), Q2(t, 0)) == 0
            chi1 = TL.chi_herm(X, 1)                 # chi / W
            r, v, Ng = x + z, z - x, y * y - t * w * w
            n = r * q * v - Ng
            Pi1 = (chi1 / 2) ** 2 - (q * n) ** 2     # Pi / W^2
            e2 = mm(eps, eps)
            g2 &= sum((e2[i][i] for i in range(4)), Q2(t, 0)) == 4 * chi1
            g2 &= det4(eps) == 4 * Pi1
            cP1 = 2 * q * (q * x * x + y * y)        # chi(c_P) / W
            cPt1 = 2 * q * (q * z * z + t * w * w)   # chi(c_P~) / W
            mp, mn = 2 * cP1, 2 * cPt1
            cp2 = charpoly(e2)
            want = [mp * mp * mn * mn, -2 * mp * mn * (mp + mn), mp * mp + 4 * mp * mn + mn * mn, -2 * (mp + mn), 1]
            g3 &= all(cp2[k] == want[k] for k in range(5))
            Y = Herm(t, rfr(rng), Q2(t, rfr(rng), rfr(rng)), Q2(t, rfr(rng), rfr(rng)), rfr(rng),
                     Q2(t, rfr(rng), rfr(rng)), rfr(rng))
            SY = to_S(Y)
            g4 &= all(SY[i][j] == SY[j][i] for i in range(4) for j in range(4))
            g4 &= sum((mm(SY, b)[i][i] for i in range(4)), Q2(t, 0)) == 0
    check("S1", g1)
    check("S2", g2)
    check("S3", g3)
    check("S4", g4)
    # (S5) principal free lattices
    g5 = True
    for t in (5, 2, 3, 13):
        one = SS.Fe(t, 1)
        L = TL.free_lattice(t, one, one)
        om = Q2(t, Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else Q2(t, 0, 1)
        z0, o1 = Q2(t, 0), Q2(t, 1)
        cols = [(o1, z0, z0, z0), (z0, z0, z0, o1), (z0, o1, o1, z0), (z0, om, om.conj(), z0)]
        Fm = [[cols[J][I] for J in range(4)] for I in range(4)]          # columns f_J in e-coordinates
        Fi = inv4(Fm)
        b = bmat(t)
        bf = mm(tr_(Fm), mm(b, Fm))
        g5 &= L.W == 4 and det4(bf) == dF(t)
        nint = 0
        for _ in range(100):
            m, a, k = (comb(rng, B) for B in (L.L4, L.J2, L.J6))
            X = Herm(t, rng.randint(-3, 3), Q2(t, a[0], a[1]), Q2(t, m[0], m[1]), m[2], Q2(t, k[0], k[1]),
                     Fr(rng.randint(-5, 5), L.W))
            assert L.integral(X)
            Sf = mm(Fi, mm(to_S(X), tr_(Fi)))
            g5 &= all(e.is_rational() and e.u.denominator == 1 for row in Sf for e in row)
            eps = mm(to_S(X), b)
            g5 &= det4(eps) == det4(Sf) * dF(t)
            g5 &= L.W ** 2 / 4 * det4(eps).u == 4 * dF(t) * det4(Sf).u
            nint += 1
        print("    t=%d principal: det b_f = %s, %d integral X: S_f integral" % (t, det4(bf), nint), flush=True)
    check("S5", g5)
    # (S6) content of Pi on L_Herm, lattice by lattice
    rows = []
    for t in (5, 2, 3, 13):
        gs = SS.reps(t, 9)[:3]
        lats = [TL.free_lattice(t, a, c) for a in gs for c in gs]
        for N in (2, 3):
            for k in range(2):
                a, c = rng.choice(gs), rng.choice(gs)
                lats.append(TL.y_lattice(t, a, c, SS.random_lattice(rng, N), "rand(N=%d)#%d %s,%s" % (N, k, a, c)))
        for L in lats:
            cont = content_Pi(L)
            rows.append((t, L.label, L.W, cont / (4 * dF(t))))
            print("    t=%2d %-40s W=%-8s content/(4 d_F) = %s" % (t, L.label, L.W, cont / (4 * dF(t))), flush=True)
    check("S6 (content divisible by 4 d_F on every listed lattice)", all(r[3].denominator == 1 for r in rows))
    print("ALL CHECKS PASS:", ok)


if __name__ == "__main__":
    main()
