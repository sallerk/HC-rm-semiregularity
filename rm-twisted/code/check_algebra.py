"""check_algebra.py -- exact checks of the algebra of NOTE.md Section 2 (the Hermitian model of twisted secant classes).

 (A1) the Hermitian model reproduces shape_search.bclasses: the classes alpha, beta, alpha~, sqrt t beta~ as forms;
 (A2) twisting: the form of E_h X sigma(E_h)^T equals e^{-L_h} times the form of X (wedge products);
 (A3) chi_herm(X) equals the Mukai square of the form (random X, random lattices);
 (A4) degreewise integrality via the lattices J2, L4, J6, (1/W)Z equals integrality of the form (random X with
      small denominators, both outcomes occurring);
 (A5) normal form for rank r != 0: X = r w (x) sigma(w) + (n/r) t (x) sigma(t) with w = e^{-h Theta}(C + gamma Theta),
      t = e^{-h Theta} Theta, gamma = g/r, n = r q v - N(g); r^2 chi = 2 E4 + int D^2 with D = -r^2 delta_2 =
      Q_m - n T0 (m = (r^2 q + g^2)/2), E4 = W (N(m) + n N(g)); chi = (W/r^2)(4 N(m) + 2 n N(g) + n^2); and the
      twist invariance of D, r^3 delta_3 = -n T_g, E4 (computed from e^{-L_nu} X / r for the twisted X);
 (A6) chi of a pure tensor lam w (x) sigma(w) is lam^2 W N(beta(w, w)); the cone of B is q r v = N(g).

Usage: python -I -B check_algebra.py         Expected last line "ALL CHECKS PASS: True".
"""
import importlib.util
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

ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print("(%s) %s" % (name, bool(cond)), flush=True)


def expo(L):
    out, term = {0: Fr(1)}, {0: Fr(1)}
    for k in range(1, 5):
        term = {m: v / k for m, v in SS.wmul(term, L).items()}
        out = SS.add((1, out), (1, term))
    return out


def rfr(rng, den=6, num=6):
    return Fr(rng.randint(-num, num), rng.randint(1, den))


def rq2(rng, t, den=6, num=6):
    return Q2(t, rfr(rng, den, num), rfr(rng, den, num))


def comb(rng, B):
    k = [rng.randint(-2, 2) for _ in range(len(B[0]))]
    return [sum(B[i][j] * k[j] for j in range(len(k))) for i in range(len(B))]


def lattices(rng):
    out = []
    for t in (5, 2, 3, 13):
        gs = SS.reps(t, 9)[:2]
        out.append(TL.free_lattice(t, gs[0], gs[0]))
        out.append(TL.free_lattice(t, gs[0], gs[1]))
        g1, g2 = rng.choice(SS.reps(t, 9)[:5]), rng.choice(SS.reps(t, 9)[:5])
        out.append(TL.y_lattice(t, g1, g2, SS.random_lattice(rng, 3), "rand(%d)" % t))
    return out


def main():
    rng = random.Random(20261025)
    lats = lattices(rng)
    # (A1)
    good = True
    for L in lats[:4]:
        t = L.t
        for q in (Fr(1), Fr(3, 4), Fr(2, 5)):
            cl = SS.bclasses(t, q, L.P)
            for k, (x, y, z, w) in enumerate(((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))):
                X = TL.secant(t, q, x, y, z, w)
                good &= L.form(X) == {kk: v for kk, v in cl[k].items() if v != 0}
    check("A1", good)
    # (A2), (A3)
    g2 = g3 = True
    n2 = 0
    for L in lats:
        t = L.t
        for _ in range(6):
            X = Herm(t, rfr(rng), rq2(rng, t), rq2(rng, t), rfr(rng), rq2(rng, t), rfr(rng))
            h = rq2(rng, t, 4, 3)
            Lh = SS.add((h.u, L.P[(1, 0)]), (h.v, L.P[(0, 1)]))
            lhs = L.form(TL.twist(X, h))
            rhs = SS.wmul(expo({k: -v for k, v in Lh.items()}), L.form(X))
            g2 &= lhs == rhs
            fX = L.form(X)
            g3 &= L.chi(X) == SS.mukai(fX, fX, L.s0)
            n2 += 1
    check("A2", g2)
    check("A3", g3)
    # (A4) integrality, both outcomes
    g4, nyes, nno = True, 0, 0
    for L in lats:
        t = L.t
        for _ in range(60):
            # random integral combination of the degree lattices (bases = columns), then perturb half of the time
            m, a, b = (comb(rng, B) for B in (L.L4, L.J2, L.J6))
            X = Herm(t, rng.randint(-2, 2), Q2(t, a[0], a[1]), Q2(t, m[0], m[1]), m[2], Q2(t, b[0], b[1]),
                     Fr(rng.randint(-3, 3), L.W))
            if rng.random() < 0.5:
                X = X + Herm(t, 0, rq2(rng, t, 3, 1), rq2(rng, t, 3, 1), rfr(rng, 3, 1), rq2(rng, t, 3, 1),
                             rfr(rng, 3, 1))
            fast = L.integral(X)
            slow = all(v.denominator == 1 for v in L.form(X).values())
            g4 &= fast == slow
            nyes += fast
            nno += not fast
    check("A4", g4 and nyes > 20 and nno > 20)
    print("    integral/non-integral samples:", nyes, nno)
    # (A5) normal form and invariants for r != 0
    g5 = True
    for L in lats:
        t, W = L.t, L.W
        for _ in range(8):
            q = Fr(rng.randint(1, 9), rng.randint(1, 6))
            x, y, z, w = rfr(rng), rfr(rng), rfr(rng), rfr(rng)
            if x + z == 0:
                continue
            h = rq2(rng, t, 4, 3)
            r, v, g = x + z, z - x, Q2(t, y, w)
            n = r * q * v - g.norm()
            gam = g / r
            X = TL.twist(TL.secant(t, q, x, y, z, w), h)
            E = TL.E_twist(t, h)
            wv = TL.mvec(E, (Q2(t, 1), gam, Q2(t, -q / 2)))
            tv = TL.mvec(E, (Q2(t, 0), Q2(t, 1), Q2(t, 0)))
            Y = TL.outer(t, wv, r) + TL.outer(t, tv, n / r)
            g5 &= Y.key() == X.key()
            m = (Q2(t, r * r * q) + g * g) * Fr(1, 2)
            chi = L.chi(X)
            E4 = W * (m.norm() + n * g.norm())
            ID2 = W * (2 * m.norm() + n * n)
            g5 &= r * r * chi == 2 * E4 + ID2
            g5 &= chi == W / (r * r) * (4 * m.norm() + 2 * n * g.norm() + n * n)
            # twist invariance: N = e^{-L_nu} X / r with nu = c_1/r
            nu = X.x10 * Fr(1) / r
            N = TL.twist(X, nu).scale(Fr(1) / r)
            g5 &= N.x00 == 1 and N.x10 == 0
            g5 &= (N.x20 * (-r * r)) == m and N.x11 * (-r * r) == -n          # D = Q_m - n T0
            g5 &= N.x21 * (r ** 3) == -n * g                                  # r^3 delta_3 = -n T_g
            g5 &= N.x22 * (r ** 4) * W == E4
    check("A5", g5)
    # (A6) pure tensors
    g6 = True
    for L in lats[:6]:
        t = L.t
        for _ in range(10):
            wv = (rq2(rng, t), rq2(rng, t), rq2(rng, t))
            lam = rfr(rng)
            P = TL.outer(t, wv, lam)
            g6 &= L.chi(P) == lam * lam * L.W * TL.beta_F(wv, wv).norm()
        for _ in range(10):
            q = Fr(rng.randint(1, 9), rng.randint(1, 6))
            r = rfr(rng)
            if r == 0:
                continue
            g = rq2(rng, t)
            v = g.norm() / (q * r)
            x, z = (r - v) / 2, (r + v) / 2
            X = TL.secant(t, q, x, g.u, z, g.v)
            M = X.mat()
            minors = [M[i][k] * M[j][l] - M[i][l] * M[j][k] for i in range(3) for j in range(3) for k in range(3)
                      for l in range(3)]
            g6 &= all(mm == 0 for mm in minors)
    check("A6", g6)
    print("ALL CHECKS PASS:", ok)


if __name__ == "__main__":
    main()
