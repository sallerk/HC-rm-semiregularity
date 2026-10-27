"""p88_shape.py -- Phase 88: the one shape left open in "Cor 5.6 for every quartic CM field".

For q in F totally positive, an integral class c of B(Theta, q) on an abelian fourfold with real multiplication by
F = Q(sqrt t) and any lattice has chi(c) >= 8, except possibly for
    c = +-(E1^q + e Theta_1 Theta_2),   E1^q = 1 - Q_q/2 + (N(q)/4) Theta_1^2 Theta_2^2,
with W N(q) = 4 and W e^2 = 2 (W = int Theta_1^2 Theta_2^2), where chi(c) = 6 and Pi(c) = 1.  In the Hermitian model
c = Herm(r = +-1, x20 = m, x11 = n, x22 = r/W) with m = -r q/2, so c is integral iff (m, n) lies in the lattice L4 of
integral degree-4 classes Q_m + n Theta_1 Theta_2, with W N(m) = 1 and W n^2 = 2 (degrees 0 and 8 are then integral).
So the shape needs W/2 to be a rational square and 2 to be a norm from F.

For each lattice (Lambda, E): W, whether W/2 is a square, and all m with N(m) = 1/W and (m, +-sqrt(2/W)) in L4:
D m lies in Z[sqrt t] (D = common denominator of L4 in F), N(D m) = D^2/W, and the solution set is stable under
m -> u^4 m for the units u of the multiplier ring O(Lambda) (pull-back by x -> u x), so it suffices to test
+- eps1^i alpha / D, alpha of norm D^2/W in a fundamental domain of eps1 (the unit of norm 1 generating the
norm-one units of Z[sqrt t]), 0 <= i < s, u^4 = eps1^s.  Controls: the solver finds the degree-4 part
(-1/2, -1) of alpha (q = 1) on the principal free lattice; Pi(c) = 1 and chi(c) = 6 for the shape (exact,
via check_sym2.to_S); a brute-force box search agrees with the solver on a few lattices.
Lattices: free types (split, control), the non-split lattices of Lemma 6.7 (O^4, O = Z + f O_F) and superlattices,
random lattices O_F^4 + sum Z gamma_i / N; fields t = 2, 7, 14, 17 (2 is a norm from these).
Usage: python -I -B p88_shape.py OUT.txt
"""
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
NP = load("nonproduct")
TL = CS.TL
SS = TL.SS
Q2, Herm = TL.Q2, TL.Herm


def is_sq(x):
    x = Fr(x)
    if x < 0:
        return None
    a, b = math.isqrt(x.numerator), math.isqrt(x.denominator)
    return Fr(a, b) if a * a == x.numerator and b * b == x.denominator else None


def fund_unit(t):
    """fundamental unit of Z[sqrt t] (a + b sqrt t, b > 0 least)."""
    b = 1
    while True:
        for s in (-1, 1):
            a2 = t * b * b + s
            a = math.isqrt(a2)
            if a * a == a2:
                return Q2(t, a, b)
        b += 1


def qpow(x, k):
    out = Q2(x.t, 1)
    for _ in range(k):
        out = out * x
    return out


class Lattice:
    """Lat plus the action of sqrt t in lattice coordinates (for the multiplier ring)."""

    def __init__(s, t, g1, g2, A, label):
        if A is None:
            A = [[Fr(int(i == j)) for j in range(8)] for i in range(8)]
        M1, Ms, SY, P, s0, T12, W = SS.ylattice(t, g1, g2, A)
        s.L = TL.Lat(t, P, s0, label)
        s.SY, s.t, s.label = SY, t, label

    def preserves(s, a):
        """does multiplication by a = u + v sqrt t map Lambda into Lambda?"""
        return all((a.u * (i == j) + a.v * s.SY[i][j]).denominator == 1 for i in range(8) for j in range(8))


def unit_period(lat, eps0):
    j = 1
    while not lat.preserves(qpow(eps0, j)):
        j += 1
    return j


def norm_reps(t, N0, eps1):
    """alpha = a + b sqrt t in Z[sqrt t] with N(alpha) = N0 > 0 and 1 <= |s1/s2| < eps1^2 (with a margin)."""
    e1 = eps1.emb(0)
    bmax = int((e1 + 1) * math.sqrt(N0) / (2 * math.sqrt(t))) + 2
    out = []
    for b in range(-bmax, bmax + 1):
        a2 = N0 + t * b * b
        a = math.isqrt(a2)
        if a * a != a2:
            continue
        for aa in {a, -a}:
            al = Q2(t, aa, b)
            r = abs(al.emb(0) / al.emb(1))
            if 1 - 1e-9 <= r < e1 * e1 * (1 + 1e-9):
                out.append(al)
    return out


def solve(lat, Nm, ns, eps0):
    """all m with N(m) = Nm and (m, n) in L4 for some n in ns, up to m -> u^4 m (u in O(Lambda)^x)."""
    L, t = lat.L, lat.t
    D = 1
    for i in range(2):
        for j in range(3):
            D = D * L.L4[i][j].denominator // math.gcd(D, L.L4[i][j].denominator)
    N0 = Fr(D * D) * Nm
    if N0.denominator != 1:
        return [], 0, None
    N0 = int(N0)
    eps1 = eps0 if eps0.norm() == 1 else eps0 * eps0
    j = unit_period(lat, eps0)
    s = 4 * j if eps0.norm() == 1 else 2 * j
    reps = norm_reps(t, N0, eps1)
    sols, ntest = [], 0
    for al in reps:
        x = al * Fr(1, D)
        for i in range(s):
            for sg in (1, -1):
                m = x * sg
                for n in ns:
                    ntest += 1
                    if L.in_L4(m, n):
                        sols.append((m, n))
            x = x * eps1
    return sols, ntest, (D, N0, j, s, len(reps))


def brute(lat, Nm, ns, box):
    """box search over the L4 basis coefficients (cross-check of solve)."""
    L = lat.L
    out = []
    rng = range(-box, box + 1)
    for k0 in rng:
        for k1 in rng:
            for k2 in rng:
                a = sum(L.L4[0][j] * k for j, k in enumerate((k0, k1, k2)))
                b = sum(L.L4[1][j] * k for j, k in enumerate((k0, k1, k2)))
                c = sum(L.L4[2][j] * k for j, k in enumerate((k0, k1, k2)))
                m = Q2(lat.t, a, b)
                if c in ns and m.norm() == Nm:
                    out.append((m, c))
    return out


def shape_class(t, r, m, n, W):
    """r E1^q + n Theta_1 Theta_2 with x20 = m = -r q/2 (q = -2 m / r)."""
    return Herm(t, r, Q2(t, 0), m, n, Q2(t, 0), Fr(r) / W)


def Pi_of(X, W):
    b = CS.bmat(X.t)
    d = CS.det4(CS.mm(CS.to_S(X), b))
    assert d.is_rational()
    return W * W / 4 * d.u


def hnf_cols(gens):
    D = 1
    for v in gens:
        for x in v:
            D = D * x.denominator // math.gcd(D, x.denominator)
    H = SS.flint.fmpz_mat([[int(x * D) for x in v] for v in gens]).hnf()
    rows = [[Fr(int(H[r, c]), D) for c in range(8)] for r in range(H.nrows()) if any(H[r, c] != 0 for c in range(8))]
    assert len(rows) == 8
    return [[rows[c][r] for c in range(8)] for r in range(8)]


def superlattice(A, gam):
    gens = [[A[i][j] for i in range(8)] for j in range(8)]
    gens.append([sum(A[i][k] * gam[k] for k in range(8)) for i in range(8)])
    return hnf_cols(gens)


NORM2 = {2: Q2(2, 0, 1), 7: Q2(7, 3, 1), 14: Q2(14, 4, 1), 17: Q2(17, Fr(5, 2), Fr(1, 2))}


def fe(t, x):
    """Q2 -> shape_search.Fe."""
    return SS.Fe(t, x.u, x.v)


def lattices(t, rng, nrand):
    one = Q2(t, 1)
    g2s = [one, NORM2[t]]
    out = []
    for g2 in g2s:
        out.append((one, g2, None, "free(1,%s)" % g2))
    for f in (2, 4):
        for e24 in (False, True):
            A = NP.lattice_matrix(t, f, e24)
            for g2 in g2s:
                out.append((one, g2, A, "O^4 f=%d%s (1,%s)" % (f, " +e24" if e24 else "", g2)))
            for _ in range(2):
                gam = [Fr(rng.randrange(f), f) for _ in range(8)]
                A2 = superlattice(A, gam)
                for g2 in g2s:
                    out.append((one, g2, A2, "O^4 f=%d%s + Zgam/%d (1,%s)" % (f, " +e24" if e24 else "", f, g2)))
    Ns = {2: (2, 4, 8), 7: (2, 4, 7, 14), 14: (2, 4, 7), 17: (2, 4, 17)}[t]
    for N in Ns:
        for _ in range(nrand):
            A = SS.random_lattice(rng, N)
            for g2 in g2s:
                out.append((one, g2, A, "rand N=%d (1,%s)" % (N, g2)))
    return out


def main():
    out = open(sys.argv[1], "w", encoding="utf-8") if len(sys.argv) > 1 else sys.stdout

    def say(*a):
        msg = " ".join(str(x) for x in a)
        print(msg, flush=True)
        if out is not sys.stdout:
            out.write(msg + "\n")
            out.flush()

    ok = True
    # controls
    t = 2
    eps0 = fund_unit(t)
    lat = Lattice(t, SS.Fe(t, 1), SS.Fe(t, 1), None, "principal")
    sols, ntest, info = solve(lat, Fr(1, 4), [Fr(1), Fr(-1)], eps0)
    c1 = any(m == Q2(t, Fr(-1, 2)) and n == -1 for m, n in sols)
    say("(control 1) principal t=2, W=%s: solver finds (-1/2, -1) (degree-4 part of alpha, q=1): %s; %d solutions, "
        "info %s" % (lat.L.W, c1, len(sols), info))
    ok &= c1
    for (tt, m, W) in ((2, Q2(2, Fr(1, 2), Fr(1, 4)), Fr(8)), (7, Q2(7, Fr(3, 4), Fr(1, 4)), Fr(8)),
                       (17, Q2(17, Fr(5, 8), Fr(1, 8)), Fr(8)), (2, Q2(2, Fr(1, 4), Fr(1, 8)), Fr(32))):
        assert m.norm() * W == 1 and m.totally_positive()
        n = is_sq(2 / W)
        for r, nn in ((-1, n), (-1, -n), (1, n)):
            X = shape_class(tt, r, m, nn, W)
            chi, Pi = TL.chi_herm(X, W), Pi_of(X, W)
            say("(control 2) t=%d r=%d m=%s (q=%s) W=%s n=%s: chi=%s Pi=%s" % (tt, r, m, m * Fr(-2, r), W, nn, chi,
                                                                             Pi))
            if r == -1:
                ok &= chi == 6 and Pi == 1
    rng = random.Random(20261026)
    total = {"lattices": 0, "W/2 square": 0, "tests": 0, "solutions": 0}
    brute_done = 0
    for t in (2, 7, 14, 17):
        eps0 = fund_unit(t)
        for g1, g2, A, label in lattices(t, rng, {2: 25, 7: 10, 14: 8, 17: 8}[t]):
            try:
                lat = Lattice(t, fe(t, g1), fe(t, g2), A, label)
            except AssertionError:
                say("t=%d %s: skipped (setup assertion)" % (t, label))
                continue
            W = lat.L.W
            total["lattices"] += 1
            n0 = is_sq(2 / W)
            if n0 is None:
                say("t=%2d %-36s W=%-6s W/2 not a square" % (t, label, W))
                continue
            total["W/2 square"] += 1
            sols, ntest, info = solve(lat, Fr(1) / W, [n0, -n0], eps0)
            total["tests"] += ntest
            total["solutions"] += len(sols)
            line = "t=%2d %-36s W=%-6s n0=%-5s (D,N0,j,s,#reps)=%s tests=%d SOLUTIONS=%s" % (
                t, label, W, n0, info, ntest, [(str(m), str(n)) for m, n in sols])
            if brute_done < 6 and info is not None:
                b = brute(lat, Fr(1) / W, [n0, -n0], 6)
                agree = (len(b) == 0) == (len(sols) == 0)
                line += " brute(box 6): %d found, agree %s" % (len(b), agree)
                ok &= agree
                brute_done += 1
            say(line)
    say("TOTAL %s" % total)
    say("ALL CONTROLS PASS: %s" % ok)


if __name__ == "__main__":
    main()
