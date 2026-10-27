"""cm_fields.py -- secant classes for every quartic CM field F(sqrt -q), q in F totally positive (NOTE.md 3.6).

For q in F totally positive, B(Theta, q) = (p^perp (x) sigma(p)^perp) n Herm, p = (1, 0, q/2), is the 4-dimensional
rational space of
    X = [[c00, c01, -c00 sigma(q)/2], [sigma(c01), c11, -sigma(c01 q)/2], [-c00 q/2, -c01 q/2, c00 N(q)/4]]
with c00, c11 in Q and c01 in F (for rational q: c00 = r, c01 = gbar, c11 = s of NOTE.md Section 2).
Checks, in exact arithmetic:
 (M1) the four basis classes lie in p^perp (x) sigma(p)^perp; for rational q they span the space of TL.secant;
 (M2) on random twisted classes X' = e^{-L_h} X (random q, h, X; N(q) a square or not): tr eps^2 = 4 chi_1 (W = 1),
      eps^4 - (tr eps^2/2) eps^2 + det(eps) I = 0, chi_1 > 0, and disc = (tr eps^2/2)^2 - 4 det eps is N(q) times a
      rational square.  So the eigenvalues m_+-, of eps^2 lie in Q(sqrt N(q)) and are conjugate when N(q) is not a
      rational square;
 (M3) on free lattices (Theta normalized as in examples.py): for q in F totally positive of small height with N(q)
      not a rational square, and h over (1/M) J2 / J2 (M <= MMAX), the least chi of a nonzero integral class of
      e^{-L_h} B(Theta, q) (LLL + Fincke-Pohst, every hit re-checked by Lat.integral).  It is compared with
      s_cm(d_F) = 2 ceil(2 sqrt d_F), the bound of Theorem 3.9; the probe itself is HEURISTIC (finite box).
Usage: python -I -B cm_fields.py [MMAX [QD [EMAX]]]      Expected last line "ALL CHECKS PASS: True".
(q = (a + b sqrt t)/d with d <= QD, |b| <= 4 and both embeddings <= EMAX; h in (1/M) J2 / J2 with M <= MMAX.)"""
import importlib.util
import math
import os
import random
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CS = load("check_sym2")
TS = load("twist_scan")
TL = CS.TL
SS, Q2, Herm = TL.SS, TL.Q2, TL.Herm


def bclass(t, q, c00, c01, c11):
    """the class of B(Theta, q) with parameters (c00, c01, c11); q, c01 in F."""
    c00, c11 = Fr(c00), Fr(c11)
    return Herm(t, c00, c01.conj(), q * (-c00 / 2), c11, c01 * q * Fr(-1, 2), c00 * q.norm() / 4)


def basis4(t, q):
    z = Q2(t, 0)
    return [bclass(t, q, 1, z, 0), bclass(t, q, 0, Q2(t, 1), 0), bclass(t, q, 0, Q2(t, 0, 1), 0),
            bclass(t, q, 0, z, 1)]


def is_sq(x):
    x = Fr(x)
    if x < 0:
        return False
    return math.isqrt(x.numerator) ** 2 == x.numerator and math.isqrt(x.denominator) ** 2 == x.denominator


def scm(dF):
    return 2 * math.ceil(2 * math.sqrt(dF) - 1e-12)


def vec9(X):
    return [X.x00, X.x10.u, X.x10.v, X.x20.u, X.x20.v, X.x11, X.x21.u, X.x21.v, X.x22]


def solve(A, b):
    """A x = b over Q (A square, invertible), by flint."""
    m = SS.flint.fmpq_mat([[SS.flint.fmpq(x.numerator, x.denominator) for x in r] for r in A])
    v = SS.flint.fmpq_mat([[SS.flint.fmpq(x.numerator, x.denominator)] for x in b])
    s = m.solve(v)
    return [Fr(int(s[i, 0].p), int(s[i, 0].q)) for i in range(len(b))]


def integral_basis(L, Ys):
    """a Z-basis of the integral classes in the Q-span of the four classes Ys."""
    Bk = CS.lherm_basis(L)
    Bm = [[vec9(B)[i] for B in Bk] for i in range(9)]
    cols = [solve(Bm, vec9(Y)) for Y in Ys]          # coordinates of Y in the basis of L_Herm
    rows = [[cols[k][i] for k in range(4)] for i in range(9)]
    D = 1
    for r in rows:
        for x in r:
            D = D * x.denominator // math.gcd(D, x.denominator)
    H = SS.flint.fmpz_mat([[int(x * D) for x in r] for r in rows]).hnf()
    R = [[Fr(int(H[i, j]), D) for j in range(4)] for i in range(4)]
    assert all(H[i, j] == 0 for i in range(4, 9) for j in range(4))
    m = SS.flint.fmpq_mat([[SS.flint.fmpq(x.numerator, x.denominator) for x in r] for r in R]).inv()
    Rinv = [[Fr(int(m[i, j].p), int(m[i, j].q)) for j in range(4)] for i in range(4)]
    out = []
    for j in range(4):
        X = None
        for k in range(4):
            if Rinv[k][j]:
                term = Ys[k].scale(Rinv[k][j])
                X = term if X is None else X + term
        out.append(X)
    return out


def lll_gram(G):
    n = len(G)
    D = 1
    for r in G:
        for x in r:
            D = D * x.denominator // math.gcd(D, x.denominator)
    _, Tm = SS.flint.fmpz_mat([[int(x * D) for x in r] for r in G]).lll(transform=True, rep="gram")
    T = [[int(Tm[i, j]) for j in range(n)] for i in range(n)]
    assert abs(int(Tm.det())) == 1
    Gr = [[sum(T[i][a] * G[a][b] * T[j][b] for a in range(n) for b in range(n)) for j in range(n)]
          for i in range(n)]
    return T, Gr


def lin(Zs, k):
    X = None
    for c, Z in zip(k, Zs):
        if c:
            term = Z.scale(Fr(c))
            X = term if X is None else X + term
    return X


def eps_of(X):
    return CS.mm(CS.to_S(X), CS.bmat(X.t))


def check_m1(rng):
    ok = True
    for t in (5, 2, 3, 13):
        for _ in range(10):
            q = Q2(t, Fr(rng.randint(1, 9), rng.randint(1, 4)), Fr(rng.randint(-3, 3), rng.randint(1, 4)))
            if not q.totally_positive():
                continue
            p, ps = (Q2(t, 1), Q2(t, 0), q * Fr(1, 2)), (Q2(t, 1), Q2(t, 0), q.conj() * Fr(1, 2))
            for X in basis4(t, q):
                M = X.mat()
                ok &= all(TL.beta_F(p, [M[i][j] for i in range(3)]) == 0 for j in range(3))
                ok &= all(TL.beta_F(ps, M[i]) == 0 for i in range(3))
        for _ in range(10):
            qr = Fr(rng.randint(1, 9), rng.randint(1, 4))
            x, y, z, w = (CS.rfr(rng) for _ in range(4))
            r, v, g = x + z, z - x, Q2(t, y, w)
            ok &= TL.secant(t, qr, x, y, z, w).key() == bclass(t, Q2(t, qr), r, g.conj(), qr * v).key()
    return ok


def check_m2(rng):
    ok, n, nsq = True, 0, 0
    for t in (5, 2, 3, 13, 7, 6):
        for trial in range(40):
            q = Q2(t, Fr(rng.randint(1, 12), rng.randint(1, 4)), Fr(rng.randint(-3, 3), rng.randint(1, 4)))
            if trial % 8 == 0:
                g = Q2(t, rng.randint(-3, 3), rng.randint(-2, 2))
                if g == 0:
                    continue
                q = g * g * Fr(rng.randint(1, 5), rng.randint(1, 3))   # N(q) a square (biquadratic control)
            if not q.totally_positive():
                continue
            h = Q2(t, CS.rfr(rng, 4, 3), CS.rfr(rng, 4, 3))
            X = TL.twist(bclass(t, q, CS.rfr(rng), Q2(t, CS.rfr(rng), CS.rfr(rng)), CS.rfr(rng)), h)
            if X.is_zero():
                continue
            e = eps_of(X)
            e2 = CS.mm(e, e)
            tr2 = sum((e2[i][i] for i in range(4)), Q2(t, 0))
            det = CS.det4(e)
            assert tr2.is_rational() and det.is_rational()
            tr2, det = tr2.u, det.u
            chi1 = TL.chi_herm(X, 1)
            e4 = CS.mm(e2, e2)
            zero = all(e4[i][j] - e2[i][j] * (tr2 / 2) + (det if i == j else 0) == 0 for i in range(4)
                       for j in range(4))
            disc = (tr2 / 2) ** 2 - 4 * det
            sq = is_sq(q.norm())
            nsq += sq
            ok &= tr2 == 4 * chi1 and zero and chi1 > 0 and disc >= 0 and is_sq(disc / q.norm())
            n += 1
    return ok, n, nsq


def qlist(t, qd, bmax=4, emax=4):
    out, seen = [], set()
    for d in range(1, qd + 1):
        for b in range(-bmax, bmax + 1):
            for a in range(1, 4 * emax * d + 1):
                for half in ((0, 1) if t % 4 == 1 else (0,)):
                    q = Q2(t, Fr(2 * a + half, 2 * d) if half else Fr(a, d), Fr(2 * b + half, 2 * d) if half
                           else Fr(b, d))
                    if not q.totally_positive() or max(q.emb(0), q.emb(1)) > emax or is_sq(q.norm()):
                        continue
                    if q in seen:
                        continue
                    seen.add(q)
                    out.append(q)
    return out


def hreps(L, mmax):
    b1, b2 = (L.J2[0][0], L.J2[1][0]), (L.J2[0][1], L.J2[1][1])
    out = set()
    for M in range(1, mmax + 1):
        for a in range(M):
            for b in range(M):
                out.add((Fr(a, M) * b1[0] + Fr(b, M) * b2[0], Fr(a, M) * b1[1] + Fr(b, M) * b2[1]))
    return sorted(out)


LATTICES = {  # t: list of (g1, g2) as (u, v) pairs; (1, 0), (1, 0) is the principal free lattice
    5: [((1, 0), (1, 0)), ((Fr(-3, 2), Fr(-3, 2)), (Fr(-3, 2), Fr(-3, 2)))],
    2: [((1, 0), (1, 0)), ((0, -2), (0, -1))],
    3: [((1, 0), (1, 0)), ((0, -1), (1, 0))],
    13: [((1, 0), (1, 0)), ((Fr(-1, 2), Fr(-1, 2)), (1, 0))],
}


def make_lattice(t, g):
    (a1, b1), (a2, b2) = g
    g1, g2 = SS.Fe(t, a1, b1), SS.Fe(t, a2, b2)
    dl = SS.delta(t)
    P, s0 = SS.setup(t, dl * g1, dl * g1 * g2)
    cth = SS.content(P[(1, 0)])
    P = {(a, b): {k: v / cth ** (a + b) for k, v in P[(a, b)].items()} for (a, b) in P}
    return TL.Lat(t, P, s0, "free(%s,%s)" % (g1, g2))


def probe(mmax, qd, emax):
    ok = True
    for t, gl in LATTICES.items():
        dF = CS.dF(t)
        bound = scm(dF) + 8
        for g in gl:
            L = make_lattice(t, g)
            best, nl, hits_at_best = None, 0, []
            qs = qlist(t, qd, emax=emax)
            for q in qs:
                for (h0, h1) in hreps(L, mmax):
                    Ys = [TL.twist(X, Q2(t, h0, h1)) for X in basis4(t, q)]
                    Zs = integral_basis(L, Ys)
                    G = [[(L.chi(Zs[i] + Zs[j]) - L.chi(Zs[i]) - L.chi(Zs[j])) / 2 for j in range(4)]
                         for i in range(4)]
                    T, Gr = lll_gram(G)
                    nl += 1
                    for v, xs in TS.short(Gr, bound):
                        k = [sum(xs[i] * T[i][j] for i in range(4)) for j in range(4)]
                        X = lin(Zs, k)
                        c = L.chi(X)
                        e = eps_of(X)
                        Pi = L.W * L.W / 4 * CS.det4(e).u
                        good = (L.integral(X) and c == v and c % 2 == 0 and Pi > 0 and Pi % (4 * dF) == 0
                                and c * c >= 16 * dF)
                        ok &= good
                        if best is None or c < best:
                            best, hits_at_best = c, []
                        if c == best and len(hits_at_best) < 3:
                            hits_at_best.append((q, (h0, h1), Pi, "cyclic" if is_sq(q.norm() / dF) else "non-Galois"))
            print("t=%2d %-34s W=%-5s #q=%d lattices=%d: least chi %s (bound s_cm(d_F) = %d; search bound %d)"
                  % (t, L.label, L.W, len(qs), nl, best, scm(dF), bound), flush=True)
            for q, h, Pi, kind in hits_at_best:
                print("      q=%s (%s field, N(q)=%s) h=(%s,%s) Pi=%s = 4 d_F * %s"
                      % (q, kind, q.norm(), h[0], h[1], Pi, Pi / (4 * dF)), flush=True)
            ok &= best is None or best >= scm(dF)
    return ok


def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    qd = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    emax = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    rng = random.Random(20261024)
    m1 = check_m1(rng)
    print("(M1) basis in p^perp (x) sigma(p)^perp; agrees with TL.secant for rational q:", m1, flush=True)
    m2, n, nsq = check_m2(rng)
    print("(M2) eps^2 has two conjugate eigenvalues in Q(sqrt N(q)), chi > 0 (%d classes, %d with N(q) a square):"
          % (n, nsq), m2, flush=True)
    m3 = probe(mmax, qd, emax)
    print("(M3) least chi >= s_cm(d_F) and every hit integral, even, 4 d_F | Pi > 0:", m3, flush=True)
    print("ALL CHECKS PASS:", m1 and m2 and m3)


if __name__ == "__main__":
    main()
