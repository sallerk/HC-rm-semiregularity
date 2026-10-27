"""Scratch probe (heuristic; not part of rm-semiregularity): integral points of the line-bundle twists e^{-L_h} B of
Markman's secant space on abelian fourfolds with real multiplication, in the model of
rm-semiregularity/code/shape_search.py (Lambda = O_F^4, e = d1 J + d2 J, Theta = Tr e,
L_s = Tr(sqrt t e)).

Question.  NOTE Theorem 7.4(b) and Corollary 9.10 concern a translation-equivariant sheaf F with ch(F) in B.  In
[M2] the translation-equivariant sheaf is I = I(Theta) (x) O(-Theta), whose Chern character is a line-bundle twist
of a secant class.  The count of 7.4 is twist-invariant (chi and the contraction ranks do not change), but the
integrality of Section 9 is not.  So: does e^{-L} B, with L = L_h rational (h in F), contain an integral class c'
with chi(c') <= 6, off P, P~ (where 7.4(a) excludes everything anyway)?  For K(c) = K(B) the count needs
chi_Y in {2, 4, 6}; on the rank-one cone (q r v = N(g)) it needs chi_Y in {2, 4}.

Usage: python -I -B twist_scan.py OUT.json [NMAX] [QDEN]
"""
import importlib.util
import json
import math
import sys
import time
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8")
SRC = __file__.replace("twist_scan.py", "../../rm-semiregularity/code/shape_search.py")
spec = importlib.util.spec_from_file_location("shape_search", SRC)
SS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SS)


def expo(L):
    out, term = {0: Fr(1)}, {0: Fr(1)}
    for k in range(1, 5):
        term = {m: v / k for m, v in SS.wmul(term, L).items()}
        out = SS.add((1, out), (1, term))
    return out


def integral(A):
    return all(Fr(v).denominator == 1 for v in A.values())


def gs(G):
    n = len(G)
    mu = [[Fr(0)] * n for _ in range(n)]
    B = [Fr(0)] * n
    for i in range(n):
        for j in range(i):
            mu[i][j] = (G[i][j] - sum(mu[j][l] * mu[i][l] * B[l] for l in range(j))) / B[j]
        B[i] = G[i][i] - sum(mu[i][l] ** 2 * B[l] for l in range(i))
    return mu, B


def lll_gram(G):
    """exact LLL (delta = 3/4) on a positive definite rational Gram matrix; returns (T, G') with G' = T G T^t."""
    n = len(G)
    T = [[int(i == j) for j in range(n)] for i in range(n)]

    def gram(T):
        return [[sum(T[i][a] * G[a][b] * T[j][b] for a in range(n) for b in range(n)) for j in range(n)]
                for i in range(n)]

    k = 1
    while k < n:
        Gp = gram(T)
        mu, B = gs(Gp)
        for j in range(k - 1, -1, -1):
            qq = math.floor(mu[k][j] + Fr(1, 2))
            if qq:
                T[k] = [T[k][a] - qq * T[j][a] for a in range(n)]
                Gp = gram(T)
                mu, B = gs(Gp)
        if B[k] < (Fr(3, 4) - mu[k][k - 1] ** 2) * B[k - 1]:
            T[k], T[k - 1] = T[k - 1], T[k]
            k = max(k - 1, 1)
        else:
            k += 1
    return T, gram(T)


def short(G, bound):
    """all x != 0 in Z^n with x^t G x <= bound (Fincke-Pohst on the GS data, exact final test)."""
    n = len(G)
    mu, B = gs(G)
    out = []
    x = [0] * n

    def rec(i, partial):
        if i < 0:
            if any(x):
                v = sum(G[a][b] * x[a] * x[b] for a in range(n) for b in range(n))
                if 0 < v <= bound:
                    out.append((v, tuple(x)))
            return
        c = -sum(mu[j][i] * x[j] for j in range(i + 1, n))
        rem = (bound - partial) / B[i]
        if rem < 0:
            return
        r = math.sqrt(float(rem)) + 1e-9
        for xi in range(math.ceil(float(c) - r), math.floor(float(c) + r) + 1):
            x[i] = xi
            t = Fr(xi) - c
            np_ = partial + B[i] * t * t
            if np_ <= bound:
                rec(i - 1, np_)
        x[i] = 0

    rec(n - 1, Fr(0))
    return out


def kind(t, q, u):
    x, y, z, w = u
    r, v = x + z, z - x
    Ng = y * y - t * w * w
    if z == 0 and w == 0:
        return "P"
    if x == 0 and y == 0:
        return "P~"
    if r == 0 and y == 0 and w == 0:
        return "line"
    if q * r * v == Ng:
        return "cone"
    return "generic"


def scan(t, d1, d2, qs, hs, bound, label):
    P, s0 = SS.setup(t, d1, d2)
    found, nlat = [], 0
    for q in qs:
        cl = SS.bclasses(t, q, P)
        diag = [SS.mukai(c, c, s0) for c in cl]
        assert all(d > 0 for d in diag)
        for (h0, h1) in hs:
            E = expo(SS.add((-h0, P[(1, 0)]), (-h1, P[(0, 1)])))
            tw = [SS.wmul(E, c) for c in cl]
            Bas = SS.lattice2(tw)
            G = [[sum(Bas[a][i] * diag[a] * Bas[a][j] for a in range(4)) for j in range(4)] for i in range(4)]
            T, Gr = lll_gram(G)
            nlat += 1
            for v, xs in short(Gr, bound):
                k = [sum(xs[i] * T[i][j] for i in range(4)) for j in range(4)]
                u = tuple(sum(Bas[a][j] * k[j] for j in range(4)) for a in range(4))
                cprime = SS.lin(tw, u)
                assert integral(cprime)
                assert SS.mukai(SS.lin(cl, u), SS.lin(cl, u), s0) == v == SS.mukai(cprime, cprime, s0)
                found.append({"type": label, "q": str(q), "h": "%s+%s*sqrt%d" % (h0, h1, t), "chi": int(v),
                              "u": [str(a) for a in u], "kind": kind(t, Fr(q), u), "rank": str(u[0] + u[2])})
    return found, nlat


def main():
    out = sys.argv[1]
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    qden = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    qs = sorted({Fr(a, b) for b in range(1, qden + 1) for a in range(1, 4 * b + 1)})
    res, t0 = {}, time.time()
    for t in (5, 2, 3, 13):
        dl, w = SS.delta(t), SS.omega(t)
        hs = sorted({(Fr(a, N) + Fr(b, N) * w.u, Fr(b, N) * w.v) for N in range(1, nmax + 1)
                     for a in range(N) for b in range(N)})
        types = [("principal", dl, dl)]
        for g in SS.reps(t, 9)[1:3]:
            types.append(("(1,%s)" % g, dl, dl * g))
        rows = []
        for label, d1, d2 in types:
            f, nl = scan(t, d1, d2, qs, hs, 6, label)
            rows += f
            by = {}
            for r in f:
                by.setdefault(r["kind"], []).append(r)
            summ = {k: (min(r["chi"] for r in v), len(v)) for k, v in by.items()}
            unt = [r for r in f if r["h"].startswith("0+0*")]
            print("t=%d %-12s lattices %d: (least chi, count) by kind: %s; untwisted hits: %d  [%.0fs]"
                  % (t, label, nl, summ, len(unt), time.time() - t0), flush=True)
            for k in ("generic", "cone"):
                ex = sorted(by.get(k, []), key=lambda r: (r["chi"], r["q"]))[:3]
                for r in ex:
                    print("     %s: chi=%d q=%s h=%s u=%s rank=%s" % (k, r["chi"], r["q"], r["h"], r["u"], r["rank"]))
        res[str(t)] = rows
    with open(out, "w", newline="\n") as fh:
        json.dump({"nmax": nmax, "qden": qden, "results": res}, fh, indent=1)
    print("done in %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
