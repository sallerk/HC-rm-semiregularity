"""
shape_search.py -- checks for ../NOTE.txt, Section 9: a search for B-secant Chern characters with the numerical
shape of Theorem 7.4(b), i.e. integral classes c in Markman's secant space B, generic (nonzero components in both
the K_0-plane P and the K_1-plane), with chi = (c, c) in {2, 4, 6} on some abelian fourfold Y with real
multiplication.  A Gamma-equivariant sheaf on X (Gamma a group of translations) is a sheaf on Y = X/Gamma, and
chi_Y(c) = chi_X(c)/|Gamma|; every such Y is X/Gamma for X in the model below (H_1(X) = O_F^4 with a symplectic
O_F-basis, e = d1 J + d2 J, d_i in the inverse different, class number one).

Model.  F = Q(sqrt t), Lambda = O_F^4 with Z-basis {1, omega} in each coordinate, Theta = Tr(e), L_s = Tr(sqrt t e);
the classes alpha, beta, alpha~, sqrt t beta~ of B are polynomials in Theta and L_s (Theta_1 Theta_2 =
(Theta^2 - L_s^2/t)/4); H^*(X, Z) = wedge^* Hom(Lambda, Z) (exact rational arithmetic on bitmask forms).  The
integral classes of B form L_ev + L_odd (degrees 0, 4, 8 and 2, 6, integrality being degreewise), and the Mukai
form is diagonal in x, y, z, w (c = x alpha + y beta + z alpha~ + w sqrt t beta~).

 (H) Free O_F-lattices: for t = 5, 2, 3, 13, 25 polarization types (d_i = delta g_i, g_i of norm <= 9) and 552
     values q = a/b (b <= 12, q <= 12): the least chi of a generic integral class (Gauss-reduced enumeration of the
     rank-2 lattices L_ev, L_odd).  Expected: 12 for t = 5 (the class (beta + sqrt 5 beta~)/2 type) and nothing
     <= 14 for t = 2, 3, 13; in any case >= 8 (NOTE Theorem 9.2).
 (I) Cyclic quotients: for every generic integral class c with chi_X(c) <= 150 (16 types, 96 values of q, t = 5,
     2, 3, 13), the largest cyclic Gamma with c integral on X/Gamma has order s_max(c), the largest elementary
     divisor of lambda -> i_lambda c (Smith normal form); least chi_X(c)/s_max(c) > 6.
 (J) Elementary abelian quotients Gamma = (Z/p)^r in X[p], p = 2, 3, 5: c is integral on X/Gamma iff every
     |S|-fold contraction of c with a set S of generators is divisible by p^|S|; exact depth-first search for the
     largest Gamma (generators restricted to the kernel mod p of lambda -> i_lambda c), for every generic integral
     c with chi_X(c) <= 100 (9 types, 16 values of q, t = 5); least chi_X(c)/p^r > 6.
 (L) Odd classes on random non-free lattices Lambda_Y = Lambda_X + sum Z gamma_i/N (t = 5 with N = 2, 3, 4, 5;
     t = 2, 3, 13), q = 1, 1/2, 2/3, 3, 5/4: y beta + w sqrt t beta~ is integral iff g Lambda lies in Lambda^Theta
     and (qW/4) g Lambda^Theta in Lambda (g = y + w sqrt t, W = int Theta_1^2 Theta_2^2; both directions, with
     non-integral multiples as negative tests); then a = (qW/4) g^2 is integral on Lambda (a in O(Lambda)),
     chi = 4 Tr(a), Tr(a) >= 2, and Tr(a) >= 3 for generic g (NOTE Theorem 9.4).
 (M) The determinant lattice D = span det_F(Lambda^4) of random lattices: vol(Lambda) >= n(D) f_D (n = covolume
     relative to O_F, f_D the conductor of the multiplier ring of D), while an integral E1 would need
     vol(Lambda) <= n(D) f_D / d_F (NOTE Proposition 9.5).
 (N) On 430 random lattices (t = 5, 2, 3, 13, N up to 13): 4 d_F divides chi of the primitive class on the line
     Q Theta_1 Theta_2, and B has no nonzero integral class with chi <= 6 for q = 1, 1/2, 1/3, 2, 1/6 (NOTE
     Proposition 9.6, Theorem 9.8).
 (O) The proof of Theorem 9.8 step by step on random lattices with conductor a power of l (l = 2, 3, 5, 13): a pair
     x, y with E(x, y) generating J at l, P = <x, y>_F, Q = P^perp, the glue index iota = [Lambda : Lambda_P +
     Lambda_Q]; Lemma 9.9 v(iota) <= v(c_P) - v(d_F) - v n(J_P,m) (J_P,m = span E(pi_P Lambda, Lambda_P)), the trace
     bound 2 v(g) <= v n(J_P,m) + v n(J_Q) + v(d_F) (g the content of Theta_1 Theta_2), and v_l(chi(E2)) >=
     v_l(4 d_F).
 (P) Dimensions n = 6, 8 (NOTE Proposition 9.11), on principal and random lattices (t = 5, 57, 13 for n = 6;
     t = 11, 5, 2, 3 for n = 8), q = 1, 1/3, 3/2: a = q^(m-1) P g^2 is integral on Lambda and chi_odd =
     2^(n-2) Tr(a) (m = n/2 even) or -2^(n-1) q^(m-1) P N(g) (m odd); for n = 6 every chi on the integral lattices
     of B is divisible by 32; for n = 8 even classes of nonzero rank have chi >= 64 and every chi < 64 is 0 mod 4.
     n = 8, the line r = g = 0 (NOTE Corollary 9.12, Theorem 9.14): with T = Theta_1 Theta_2 and g_4, g_12 the
     contents of T and T^3, every integral class there has chi = 4 int c_4 c_12 >= 4 N0, N0 = int T^4/(g_4 g_12)
     (an integer, any q); N0 = 16/(g_4 g_4^v) with g_4^v the content of pi_1 pi_2 on the dual lattice (pi =
     Theta^-1), and 16 | N0, on every tested lattice (N0 >= 256 there).
 (R) NOTE Proposition 8.1(e) on the principal free lattice (t = 5, 2, 3, 13): for q = 1 the nonzero integral
     classes of B with c_1 in Q.Theta and chi <= 14 are exactly +-alpha, +-beta; alpha~ is integral first at q = t.
Usage:  python -I shape_search.py OUT.json      Expected last line "ALL CHECKS PASS: True".  About 7 minutes.
"""

import itertools
import json
import math
import random
import sys
import time
from fractions import Fraction as Fr

import flint


def popc(x):
    return bin(x).count("1")


def wmul(A, B):
    out = {}
    for m1, c1 in A.items():
        for m2, c2 in B.items():
            if m1 & m2:
                continue
            s, x = 0, m2
            while x:
                j = (x & -x).bit_length() - 1
                s += popc(m1 >> (j + 1))
                x &= x - 1
            k = m1 | m2
            out[k] = out.get(k, 0) + (-c1 * c2 if s & 1 else c1 * c2)
    return {k: v for k, v in out.items() if v != 0}


def add(*terms):
    out = {}
    for c, A in terms:
        for k, v in A.items():
            out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


class Fe:
    """u + v sqrt t."""

    def __init__(s, t, u, v=0):
        s.t, s.u, s.v = t, Fr(u), Fr(v)

    def __mul__(s, o):
        return Fe(s.t, s.u * o.u + s.t * s.v * o.v, s.u * o.v + s.v * o.u)

    def tr(s):
        return 2 * s.u

    def __repr__(s):
        return "%s+%sr%d" % (s.u, s.v, s.t)


def omega(t):
    return Fe(t, Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else Fe(t, 0, 1)


def delta(t):
    """generator of the inverse different."""
    return Fe(t, 0, Fr(1, t)) if t % 4 == 1 else Fe(t, 0, Fr(1, 2 * t))


def reps(t, maxnorm):
    w, seen, out = omega(t), set(), []
    for a in range(0, 12):
        for b in range(-12, 13):
            g = Fe(t, a + b * w.u, b * w.v)
            N = g.u * g.u - t * g.v * g.v
            if N == 0 or abs(N) > maxnorm or (abs(N), abs(g.v)) in seen:
                continue
            seen.add((abs(N), abs(g.v)))
            out.append(g)
    return out


def twoform(t, d1, d2, f):
    bs = [Fe(t, 1), omega(t)]
    out = {}
    for (a, b, d) in ((0, 1, d1), (2, 3, d2)):
        for k in range(2):
            for l in range(2):
                val = (f * d * bs[k] * bs[l]).tr()
                i, j = 2 * a + k, 2 * b + l
                if val != 0:
                    key = (1 << i) | (1 << j)
                    out[key] = out.get(key, 0) + (1 if i < j else -1) * val
    return {k: v for k, v in out.items() if v != 0}


def setup(t, d1, d2):
    L1, Ls = twoform(t, d1, d2, Fe(t, 1)), twoform(t, d1, d2, Fe(t, 0, 1))
    assert all(v.denominator == 1 for v in list(L1.values()) + list(Ls.values()))
    return setup_forms(L1, Ls)


def setup_forms(L1, Ls):
    P = {}
    for a in range(5):
        for b in range(5 - a):
            X = {0: Fr(1)}
            for _ in range(a):
                X = wmul(X, L1)
            for _ in range(b):
                X = wmul(X, Ls)
            P[(a, b)] = X
    s0 = 1 if P[(4, 0)].get(255, 0) > 0 else -1
    return P, s0


def bclasses(t, q, P):
    q = Fr(q)
    al = add((1, P[(0, 0)]), (-q / 2, P[(2, 0)]), (q * q / 24, P[(4, 0)]))
    be = add((1, P[(1, 0)]), (-q / 6, P[(3, 0)]))
    alt = add((1, P[(0, 0)]), (-q / (2 * t), P[(0, 2)]), (q * q / (24 * t * t), P[(0, 4)]))
    sbt = add((1, P[(0, 1)]), (-q / (6 * t), P[(0, 3)]))
    return [al, be, alt, sbt]


def tau(A):
    return {k: (-v if (popc(k) // 2) % 2 else v) for k, v in A.items()}


def mukai(A, B, s0):
    return s0 * wmul(tau(A), B).get(255, 0)


def lattice2(classes):
    """basis (columns) of {u : sum u_i classes_i integral}."""
    k = len(classes)
    keys = sorted(set(kk for c in classes for kk in c))
    rows = [[Fr(c.get(kk, 0)) for c in classes] for kk in keys]
    D = 1
    for r in rows:
        for x in r:
            D = D * x.denominator // math.gcd(D, x.denominator)
    R = flint.fmpz_mat([[int(x * D) for x in r] for r in rows]).hnf()
    M = flint.fmpq_mat([[flint.fmpq(int(R[i, j]), D) for j in range(k)] for i in range(k)])
    Mi = M.inv()
    return [[Fr(int(Mi[i, j].p), int(Mi[i, j].q)) for j in range(k)] for i in range(k)]


def enum2(Bas, diag, bound):
    """vectors u = Bas k with diag[0] u0^2 + diag[1] u1^2 <= bound (Gauss-reduced enumeration)."""
    def ip(u, v):
        return diag[0] * u[0] * v[0] + diag[1] * u[1] * v[1]
    b1, b2 = (Bas[0][0], Bas[1][0]), (Bas[0][1], Bas[1][1])
    while True:
        if ip(b1, b1) > ip(b2, b2):
            b1, b2 = b2, b1
        m = math.floor(ip(b1, b2) / ip(b1, b1) + Fr(1, 2))
        if m == 0:
            break
        b2 = (b2[0] - m * b1[0], b2[1] - m * b1[1])
    n1, n12 = ip(b1, b1), ip(b1, b2)
    d2 = ip(b2, b2) - n12 * n12 / n1
    out = []
    K2 = int(math.isqrt(int(bound / d2) + 1)) + 1
    for k2 in range(-K2, K2 + 1):
        rem = bound - k2 * k2 * d2
        if rem < 0:
            continue
        c, r = -k2 * n12 / n1, math.sqrt(float(rem / n1)) + 1
        for k1 in range(math.floor(float(c) - r), math.ceil(float(c) + r) + 1):
            u = (k1 * b1[0] + k2 * b2[0], k1 * b1[1] + k2 * b2[1])
            v = diag[0] * u[0] ** 2 + diag[1] * u[1] ** 2
            if 0 < v <= bound:
                out.append((v, u))
    return out


def generic_classes(t, q, P, s0, bound, generic=True):
    """generic (or, with generic=False, all nonzero) integral classes of B with chi <= bound: (chi, (x, y, z, w),
    class)."""
    cl = bclasses(t, q, P)
    al, be, alt, sbt = cl
    dev = [mukai(al, al, s0), mukai(alt, alt, s0)]
    dod = [mukai(be, be, s0), mukai(sbt, sbt, s0)]
    for i in range(4):
        for j in range(4):
            if i != j:
                assert mukai(cl[i], cl[j], s0) == 0
    ev = enum2(lattice2([al, alt]), dev, bound) + [(0, (Fr(0), Fr(0)))]
    od = enum2(lattice2([be, sbt]), dod, bound) + [(0, (Fr(0), Fr(0)))]
    out = []
    for v1, u1 in ev:
        for v2, u2 in od:
            v = v1 + v2
            x, z, y, w = u1[0], u1[1], u2[0], u2[1]
            if 0 < v <= bound and (not generic or ((x or y) and (z or w))):
                out.append((v, (x, y, z, w), cl))
    return out


def lin(cl, coeffs):
    out = {}
    for c, A in zip(coeffs, cl):
        if c:
            for k, v in A.items():
                out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def contr(lam, A):
    out = {}
    for m, c in A.items():
        x = m
        while x:
            i = (x & -x).bit_length() - 1
            if lam[i]:
                sgn = -1 if popc(m & ((1 << i) - 1)) & 1 else 1
                k = m ^ (1 << i)
                out[k] = out.get(k, 0) + sgn * lam[i] * c
            x &= x - 1
    return {k: v for k, v in out.items() if v != 0}


def contraction_rows(c):
    rows = {}
    for i in range(8):
        lam = [0] * 8
        lam[i] = 1
        for k, v in contr(lam, c).items():
            rows.setdefault(k, [0] * 8)[i] = v
    return list(rows.values())


def s_max(c):
    R = contraction_rows(c)
    assert all(Fr(x).denominator == 1 for r in R for x in r)
    S = flint.fmpz_mat([[int(x) for x in r] for r in R]).snf()
    d = [int(S[i, i]) for i in range(min(S.nrows(), S.ncols()))]
    return None if 0 in d else max(d)


def kernel_mod(rows, p):
    """F_p-basis of {lambda in F_p^8 : rows . lambda = 0 mod p}."""
    if not rows:
        return [[int(i == j) for i in range(8)] for j in range(8)]
    K = flint.nmod_mat([[int(x) % p for x in r] for r in rows], p).nullspace()[0]
    return [[int(K[i, j]) for i in range(8)] for j in range(K.ncols()) if any(int(K[i, j]) for i in range(8))]


def cond_rows(A, p, k):
    """rows of the F_p-linear map lambda -> (i_lambda A)/p^k mod p (A an integral form divisible as required)."""
    rows = {}
    for i in range(8):
        lam = [0] * 8
        lam[i] = 1
        for key, v in contr(lam, A).items():
            assert int(v) % p ** k == 0
            rows.setdefault(key, [0] * 8)[i] = (int(v) // p ** k) % p
    return list(rows.values())


def max_gamma(c, p, maxr=6):
    """largest r such that c is integral on X/Gamma for some Gamma = (Z/p)^r in X[p] (exact search).  A new
    generator gamma must satisfy, for every set S of earlier generators, i_gamma i_S c = 0 mod p^(|S|+1); given
    the earlier conditions this is F_p-linear in gamma, so candidates form a subspace."""
    best = [0]

    def span(gens):
        out = set()
        for co in itertools.product(range(p), repeat=len(gens)):
            out.add(tuple(sum(a * g[i] for a, g in zip(co, gens)) % p for i in range(8)))
        return out

    def dfs(gens, contracted, last):
        # contracted: dict S (tuple of indices) -> i_S c (exactly divisible by p^|S|)
        best[0] = max(best[0], len(gens))
        if len(gens) >= maxr:
            return
        rows = []
        for S, A in contracted.items():
            if A:
                rows += cond_rows(A, p, len(S))
        kb = kernel_mod(rows, p)
        sp = span(gens)
        cands = set()
        for co in itertools.product(range(p), repeat=len(kb)):
            v = tuple(sum(a * b[i] for a, b in zip(co, kb)) % p for i in range(8))
            if any(v) and v not in sp:
                f = pow(v[next(i for i in range(8) if v[i])], -1, p)
                cands.add(tuple((x * f) % p for x in v))
        for v in sorted(cands):
            if last is not None and v <= last:
                continue
            new = dict(contracted)
            okv = True
            for S, A in contracted.items():
                B = contr(list(v), A) if A else {}
                if any(int(x) % p ** (len(S) + 1) for x in B.values()):
                    okv = False
                    break
                new[S + (len(gens),)] = B
            if okv:
                dfs(gens + [list(v)], new, v)

    dfs([], {(): c}, None)
    return best[0]


def check_H():
    rows, ok = [], True
    for t in (5, 2, 3, 13):
        dl = delta(t)
        gs = reps(t, 9)[:5]
        qs = sorted({Fr(a, b) for b in range(1, 13) for a in range(1, 12 * b + 1) if Fr(a, b) <= 12})
        best = None
        for g1 in gs:
            for g2 in gs:
                P, s0 = setup(t, dl * g1, dl * g1 * g2)
                for q in qs:
                    for v, u, _ in generic_classes(t, q, P, s0, 14):
                        if best is None or v < best[0]:
                            best = (v, str(g1), str(g2), str(q), [str(x) for x in u])
        good = best is None or best[0] >= 8
        ok &= good
        rows.append({"t": t, "types": len(gs) ** 2, "q": len(qs), "least_generic_chi_le_14": best})
        print("    (H) t=%d: %d types x %d q: least generic chi <= 14: %s: %s"
              % (t, len(gs) ** 2, len(qs), best, good), flush=True)
    return rows, ok


def check_I():
    rows, ok = [], True
    for t in (5, 2, 3, 13):
        dl = delta(t)
        gs = reps(t, 9)[:4]
        qs = sorted({Fr(a, b) for b in (1, 2, 3, 4, 6, 12) for a in range(1, 8 * b + 1)})
        best, n = None, 0
        for g1 in gs:
            for g2 in gs:
                P, s0 = setup(t, dl * g1, dl * g1 * g2)
                for q in qs:
                    for v, u, cl in generic_classes(t, q, P, s0, 150):
                        sm = s_max(lin(cl, u))
                        n += 1
                        if sm and (best is None or Fr(v, sm) < best[0]):
                            best = (Fr(v, sm), v, sm, str(g1), str(g2), str(q), [str(x) for x in u])
        good = best[0] > 6
        ok &= good
        rows.append({"t": t, "classes": n, "least_chi_over_smax": [str(x) for x in best]})
        print("    (I) t=%d: %d generic classes with chi_X <= 150: least chi_X/s_max = %s (chi_X %s, s_max %s): %s"
              % (t, n, best[0], best[1], best[2], good), flush=True)
    return rows, ok


def check_J():
    rows, ok = [], True
    t = 5
    dl = delta(t)
    gs = reps(t, 9)[:3]
    qs = sorted({Fr(a, b) for b in (1, 2, 3, 4) for a in range(1, 4 * b + 1)})
    for p in (2, 3, 5):
        best, n, seen = None, 0, set()
        for g1 in gs:
            for g2 in gs:
                P, s0 = setup(t, dl * g1, dl * g1 * g2)
                for q in qs:
                    for v, u, cl in generic_classes(t, q, P, s0, 100):
                        key = (str(g1), str(g2), q, tuple(abs(x) for x in u))
                        if key in seen:
                            continue
                        seen.add(key)
                        r = max_gamma(lin(cl, u), p)
                        n += 1
                        if best is None or Fr(v, p ** r) < best[0]:
                            best = (Fr(v, p ** r), v, r, str(g1), str(g2), str(q), [str(x) for x in u])
        good = best[0] > 6
        ok &= good
        rows.append({"t": t, "p": p, "classes": n, "least": [str(x) for x in best]})
        print("    (J) t=5, p=%d: %d generic classes with chi_X <= 100: least chi_X/p^r = %s (chi_X %s, r %s): %s"
              % (p, n, best[0], best[1], best[2], good), flush=True)
    # the two quotients quoted in NOTE (J): gs[1] = -1 - sqrt 5 (norm 4), gs[2] = -phi (a unit)
    found = {}
    for key, g, q, u0, bound, mult in (("norm4", gs[1], Fr(1, 4), (0, Fr(1, 4), 0, Fr(1, 4)), 48, 1),
                                       ("2b_phi", gs[2], Fr(1), (0, Fr(1, 2), 0, Fr(1, 2)), 12, 2)):
        P, s0 = setup(t, dl * g, dl * g * g)
        for v, u, cl in generic_classes(t, q, P, s0, bound):
            if tuple(u) == u0:
                found[key] = (int(v) * mult ** 2, max_gamma(lin(cl, [mult * x for x in u]), 2))
    good = found == {"norm4": (48, 2), "2b_phi": (48, 1)}
    ok &= good
    rows.append({"examples": {k: [str(x) for x in v] for k, v in found.items()}})
    print("    (J) t=5, p=2: type (1 + sqrt 5, 1 + sqrt 5), q = 1/4, (beta + sqrt 5 beta~)/4: (chi_X, r) = %s; "
          "principal, q = 1, 2 b_phi: (chi_X, r) = %s: %s" % (found.get("norm4"), found.get("2b_phi"), good),
          flush=True)
    return rows, ok


# ----------------------------------------------------------------------------- general lattices, checks (L)-(N)
def mat_of(form):
    M = [[Fr(0)] * 8 for _ in range(8)]
    for k, v in form.items():
        i, j = [b for b in range(8) if k >> b & 1]
        M[i][j] += v
        M[j][i] -= v
    return M


def form_of(M):
    return {(1 << i) | (1 << j): M[i][j] for i in range(8) for j in range(i + 1, 8) if M[i][j]}


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(8)) for j in range(8)] for i in range(8)]


def tr8(A):
    return [list(r) for r in zip(*A)]


def qmat(A):
    return flint.fmpq_mat([[flint.fmpq(x.numerator, x.denominator) for x in r] for r in A])


def frmat(m):
    return [[Fr(int(m[i, j].p), int(m[i, j].q)) for j in range(m.ncols())] for i in range(m.nrows())]


def integral(A):
    return all(Fr(x).denominator == 1 for r in A for x in r)


def content(form):
    g = Fr(0)
    for v in form.values():
        v = abs(Fr(v))
        g = Fr(math.gcd(g.numerator * v.denominator, v.numerator * g.denominator), g.denominator * v.denominator) \
            if g else v
    return g


def dfield(t):
    return t if t % 4 == 1 else 4 * t


def sqrt_t_matrix(t):
    """multiplication by sqrt t on O_F^4 in the Z-basis (1, omega) of each coordinate."""
    w, cols = omega(t), []
    for fe in (Fe(t, 0, 1), Fe(t, 0, 1) * w):
        b = fe.v / w.v
        cols.append((fe.u - b * w.u, b))
    M = [[Fr(0)] * 8 for _ in range(8)]
    for k in range(4):
        for c in range(2):
            M[2 * k][2 * k + c], M[2 * k + 1][2 * k + c] = cols[c]
    return M


def random_lattice(rng, N):
    """Z-basis (columns, coordinates of Lambda_X = O_F^4) of Lambda_Y = Lambda_X + sum Z gamma_i / N (1 to 5
    random gamma_i); the multiplier ring of Lambda_Y contains Z + N O_F."""
    gens = [[Fr(int(i == j)) for i in range(8)] for j in range(8)]
    for _ in range(rng.randint(1, 5)):
        v = [rng.randrange(N) for _ in range(8)]
        if rng.random() < 0.5:
            v = [x if rng.random() < 0.4 else 0 for x in v]
        gens.append([Fr(x, N) for x in v])
    D = 1
    for v in gens:
        for x in v:
            D = D * x.denominator // math.gcd(D, x.denominator)
    H = flint.fmpz_mat([[int(x * D) for x in v] for v in gens]).hnf()
    rows = [[Fr(int(H[r, c]), D) for c in range(8)] for r in range(H.nrows())
            if any(H[r, c] != 0 for c in range(8))]
    assert len(rows) == 8
    return [[rows[c][r] for c in range(8)] for r in range(8)]


def ylattice(t, g1, g2, A):
    """Lambda_Y with basis A inside the free type (g1, g2): Theta and L_sqrt t (matrices), the action of sqrt t,
    the products P, the orientation, Theta_1 Theta_2 and W = int Theta_1^2 Theta_2^2 (all in Y coordinates)."""
    dl = delta(t)
    M1x = mat_of(twoform(t, dl * g1, dl * g1 * g2, Fe(t, 1)))
    Msx = mat_of(twoform(t, dl * g1, dl * g1 * g2, Fe(t, 0, 1)))
    Sx = sqrt_t_matrix(t)
    assert mm(tr8(Sx), M1x) == Msx
    M1, Ms = mm(tr8(A), mm(M1x, A)), mm(tr8(A), mm(Msx, A))
    SY = mm(frmat(qmat(A).inv()), mm(Sx, A))
    P, s0 = setup_forms(form_of(M1), form_of(Ms))
    T12 = add((Fr(1, 4), P[(2, 0)]), (Fr(-1, 4 * t), P[(0, 2)]))
    W = abs(wmul(T12, T12).get(255, 0))
    return M1, Ms, SY, P, s0, T12, W


def check_L():
    """the odd criterion of Theorem 9.4 on random non-free lattices."""
    rng = random.Random(20261019)
    rows, ok = [], True
    for t, N, ntr in ((5, 2, 20), (5, 3, 12), (5, 4, 12), (5, 5, 12), (2, 2, 12), (3, 3, 12), (13, 2, 8)):
        gs, n, least, good_all = reps(t, 9)[:3], 0, None, True
        for _ in range(ntr):
            g1, g2 = rng.choice(gs), rng.choice(gs)
            A = random_lattice(rng, N)
            M1, Ms, SY, P, s0, T12, W = ylattice(t, g1, g2, A)
            M1i = frmat(qmat(M1).inv())
            for q in (Fr(1), Fr(1, 2), Fr(2, 3), Fr(3), Fr(5, 4)):
                cl = bclasses(t, q, P)
                bas = lattice2([cl[1], cl[3]])
                dod = [mukai(cl[1], cl[1], s0), mukai(cl[3], cl[3], s0)]
                vs = [(bas[0][0], bas[1][0]), (bas[0][1], bas[1][1]),
                      (bas[0][0] + bas[0][1], bas[1][0] + bas[1][1])]
                tests = [(y, w, True) for y, w in vs] + [(vs[0][0] / 2, vs[0][1] / 2, None),
                                                          (vs[1][0] / 3, vs[1][1] / 3, None)]
                for y, w, known in tests:
                    G = [[y * (i == j) + w * SY[i][j] for j in range(8)] for i in range(8)]
                    crit = (integral([[y * M1[i][j] + w * Ms[i][j] for j in range(8)] for i in range(8)])
                            and integral([[q * W / 4 * x for x in r] for r in mm(G, M1i)]))
                    integ = all(v.denominator == 1 for v in lin(cl, [0, y, 0, w]).values())
                    good = crit == integ and (known is None or integ)
                    if integ:
                        a = [[q * W / 4 * x for x in r] for r in mm(G, G)]
                        tra = q * W / 4 * 2 * (y * y + t * w * w)
                        chi = dod[0] * y * y + dod[1] * w * w
                        good = good and integral(a) and chi == 4 * tra and tra >= 2 and (not (y and w) or tra >= 3)
                        if y and w and (least is None or chi < least):
                            least = chi
                    good_all &= good
                    n += 1
        ok &= good_all
        rows.append({"t": t, "N": N, "lattices": ntr, "tests": n, "least_generic_odd_chi": str(least)})
        print("    (L) t=%d N=%d: %d lattices, %d tests of the odd criterion and of chi = 4 Tr(a), a in O(Lambda);"
              " least generic odd chi %s: %s" % (t, N, ntr, n, least, good_all), flush=True)
    return rows, ok


def fdet4(t, cols):
    """det over F of a 4 x 4 matrix of Fe entries (cols[j][i] = entry (i, j))."""
    tot = Fe(t, 0)
    for perm in itertools.permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if perm[i] > perm[j]:
                    sg = -sg
        pr = Fe(t, sg)
        for j in range(4):
            pr = pr * cols[j][perm[j]]
        tot = Fe(t, tot.u + pr.u, tot.v + pr.v)
    return tot


def in_lattice2(bas, z):
    """is z (coordinates in the basis 1, omega) in the Z-span of the two vectors bas?"""
    det = bas[0][0] * bas[1][1] - bas[0][1] * bas[1][0]
    c0 = (z[0] * bas[1][1] - z[1] * bas[1][0]) / det
    c1 = (bas[0][0] * z[1] - bas[0][1] * z[0]) / det
    return c0.denominator == 1 and c1.denominator == 1


def omega_times(t, z, f):
    """f omega (z0 + z1 omega) in the basis 1, omega (omega^2 = Tr(omega) omega - N(omega))."""
    w = omega(t)
    tw, nw = 2 * w.u, w.u * w.u - t * w.v * w.v
    return (-f * z[1] * nw, f * (z[0] + z[1] * tw))


def check_M():
    """E1: the determinant lattice D = span det_F(Lambda^4) satisfies n(D) f_D <= vol(Lambda) (Proposition 9.5)."""
    rng = random.Random(1019)
    rows, ok = [], True
    for t, N, ntr in ((5, 2, 25), (5, 3, 15), (5, 5, 15), (5, 6, 10), (2, 2, 20), (3, 2, 10), (13, 13, 10)):
        w, worst, good_all = omega(t), None, True
        for _ in range(ntr):
            A = random_lattice(rng, N)
            fc = [[Fe(t, A[2 * k][c] + A[2 * k + 1][c] * w.u, A[2 * k + 1][c] * w.v) for k in range(4)]
                  for c in range(8)]
            coords = []
            for S in itertools.combinations(range(8), 4):
                dd = fdet4(t, [fc[c] for c in S])
                b = dd.v / w.v
                coords.append((dd.u - b * w.u, b))
            den = 1
            for x in coords:
                for y in x:
                    den = den * y.denominator // math.gcd(den, y.denominator)
            H = flint.fmpz_mat([[int(x * den) for x in c] for c in coords]).hnf()
            bas = [(Fr(int(H[r, 0]), den), Fr(int(H[r, 1]), den)) for r in range(2)]
            nD = abs(bas[0][0] * bas[1][1] - bas[0][1] * bas[1][0])
            fD = next(f for f in range(1, N + 1) if N % f == 0 and all(in_lattice2(bas, omega_times(t, z, f))
                                                                          for z in bas))
            dA = qmat(A).det()
            vol = abs(Fr(int(dA.p), int(dA.q)))
            r = vol / (nD * fD)
            good = r >= 1
            good_all &= good
            if worst is None or r < worst:
                worst = r
        ok &= good_all
        rows.append({"t": t, "N": N, "lattices": ntr, "least_vol_over_nD_fD": str(worst)})
        print("    (M) t=%d N=%d: %d lattices: least vol(Lambda)/(n(D) f_D) = %s >= 1 (an integral E1 needs"
              " <= 1/d_F = 1/%d): %s" % (t, N, ntr, worst, dfield(t), good_all), flush=True)
    return rows, ok


def check_N():
    """chi of the primitive class on the line Q Theta_1 Theta_2, and no nonzero integral class of B with chi <= 6,
    on random lattices."""
    rng = random.Random(191019)
    rows, ok = [], True
    for t, N, ntr in ((5, 2, 60), (5, 3, 40), (5, 4, 40), (5, 5, 60), (5, 10, 30), (2, 2, 60), (3, 2, 40),
                      (3, 3, 40), (13, 2, 30), (13, 13, 30)):
        gs, dF = reps(t, 9)[:5], dfield(t)
        least, div, small, good_all = None, True, 0, True
        for _ in range(ntr):
            g1, g2 = rng.choice(gs), rng.choice(gs)
            A = random_lattice(rng, N)
            M1, Ms, SY, P, s0, T12, W = ylattice(t, g1, g2, A)
            chi = W / content(T12) ** 2
            div &= chi.denominator == 1 and int(chi) % (4 * dF) == 0
            if least is None or chi < least:
                least = chi
            for q in (Fr(1), Fr(1, 2), Fr(1, 3), Fr(2), Fr(1, 6)):
                small += len(generic_classes(t, q, P, s0, 6, generic=False))
        good_all = div and small == 0
        ok &= good_all
        rows.append({"t": t, "N": N, "lattices": ntr, "least_chi_E2": str(least), "4dF_divides": div,
                     "classes_chi_le_6": small})
        print("    (N) t=%d N=%d: %d lattices: least chi(E2) = %s, 4 d_F = %d divides every chi(E2): %s; nonzero"
              " integral classes of B with chi <= 6 (5 values of q): %d: %s"
              % (t, N, ntr, least, 4 * dF, div, small, good_all), flush=True)
    return rows, ok


def vval(x, l):
    x = Fr(x)
    if x == 0:
        return 10 ** 6
    n = 0
    while x.numerator % l == 0:
        x, n = x / l, n + 1
    while x.denominator % l == 0:
        x, n = x * l, n - 1
    return n


def int_kernel(rows):
    """Z-basis of {z in Z^8 : rows . z = 0} for rational rows."""
    D = 1
    for r in rows:
        for x in r:
            D = D * Fr(x).denominator // math.gcd(D, Fr(x).denominator)
    m = len(rows)
    A = flint.fmpz_mat([[int(Fr(rows[i][j]) * D) for i in range(m)] + [int(j == k) for k in range(8)]
                        for j in range(8)])
    H = A.hnf()
    return [[int(H[r, m + k]) for k in range(8)] for r in range(8) if all(H[r, c] == 0 for c in range(m))]


def fe_div(a, b):
    n = b.u * b.u - b.t * b.v * b.v
    return a * Fe(a.t, b.u / n, -b.v / n)


def span_n(t, vals):
    """n(L) = covol(L)/covol(O_F) for L = Z-span of the F-elements vals, and the conductor of O(L)."""
    w = omega(t)
    co = [((z.u - (z.v / w.v) * w.u), z.v / w.v) for z in vals]
    den = 1
    for x in co:
        for y in x:
            den = den * y.denominator // math.gcd(den, y.denominator)
    H = flint.fmpz_mat([[int(x * den) for x in c] for c in co]).hnf()
    bas = [(Fr(int(H[r, 0]), den), Fr(int(H[r, 1]), den)) for r in range(2)]
    n = abs(bas[0][0] * bas[1][1] - bas[0][1] * bas[1][0])
    f = 1
    while not all(in_lattice2(bas, omega_times(t, z, f)) for z in bas):
        f += 1
    return n, f


def check_O():
    """the proof of Theorem 9.8 step by step: a plane P = <x, y> with E(x, y) generating J, the glue index iota,
    Lemma 9.9 (v(iota) <= v(c_P) - v(d) - v n(J_P,m)), the trace bound 2 v(g) <= v n(J_P,m) + v n(J_Q) + v(d),
    and v(chi(E2)) >= v(4 d_F), at the prime l of the conductor."""
    rng = random.Random(9819)
    rows, ok = [], True
    for t, N, ntr in ((5, 5, 25), (5, 25, 12), (5, 3, 12), (13, 13, 10), (3, 3, 12), (3, 2, 12), (2, 2, 12),
                      (2, 4, 8)):
        l = next(p for p in range(2, N + 1) if N % p == 0)
        dF, w, gs = dfield(t), omega(t), reps(t, 9)[:4]
        good_all, tight = True, 0
        for _ in range(ntr):
            g1, g2 = rng.choice(gs), rng.choice(gs)
            A = random_lattice(rng, N)
            dl = delta(t)
            Eb = [(0, 1, dl * g1), (2, 3, dl * g1 * g2)]

            def toF(v):
                return [Fe(t, v[2 * k] + v[2 * k + 1] * w.u, v[2 * k + 1] * w.v) for k in range(4)]

            def EF(a, b):
                tot = Fe(t, 0)
                for i, j, dd in Eb:
                    p1, p2 = dd * (a[i] * b[j]), dd * (a[j] * b[i])
                    tot = Fe(t, tot.u + p1.u - p2.u, tot.v + p1.v - p2.v)
                return tot

            basis = [[A[i][c] for i in range(8)] for c in range(8)]
            bF = [toF(v) for v in basis]
            nJ, hJ = span_n(t, [EF(bF[i], bF[j]) for i in range(8) for j in range(i + 1, 8)])
            cands = [[int(i == j) for j in range(8)] for i in range(8)]
            cands += [[a + b for a, b in zip(cands[i], cands[j])] for i in range(8) for j in range(i + 1, 8)]
            pair = None
            for cx, cy in itertools.combinations(cands, 2):
                xv = [sum(cx[c] * basis[c][i] for c in range(8)) for i in range(8)]
                yv = [sum(cy[c] * basis[c][i] for c in range(8)) for i in range(8)]
                e = EF(toF(xv), toF(yv))
                # E(x, y) generates J at l iff n(E(x, y) O(J)) = |N E(x, y)| f_J has the l-valuation of n(J)
                if (e.u or e.v) and vval(abs(e.u * e.u - t * e.v * e.v) * hJ, l) == vval(nJ, l):
                    pair = (xv, yv, e)
                    break
            if pair is None:
                good_all = False
                continue
            xv, yv, e = pair
            xF, yF = toF(xv), toF(yv)
            Sx = sqrt_t_matrix(t)
            Ai = frmat(qmat(A).inv())

            def toL(v):
                return [sum(Ai[i][j] * v[j] for j in range(8)) for i in range(8)]

            def toX(c):
                return [sum(A[i][j] * c[j] for j in range(8)) for i in range(8)]
            Pv = [toL(z) for z in (xv, [sum(Sx[i][j] * xv[j] for j in range(8)) for i in range(8)], yv,
                                   [sum(Sx[i][j] * yv[j] for j in range(8)) for i in range(8)])]
            LP = int_kernel(int_kernel(Pv))
            M1, Ms, SY, P, s0, T12, W = ylattice(t, g1, g2, A)
            LQ = int_kernel([[sum(pv[i] * M1[i][j] for i in range(8)) for j in range(8)] for pv in Pv])
            iota = abs(int(flint.fmpz_mat(LP + LQ).det()))
            LPx, LQx = [toF(toX(c)) for c in LP], [toF(toX(c)) for c in LQ]
            proj = []
            for zF in bF:
                a_, b_ = fe_div(EF(zF, yF), e), fe_div(EF(zF, xF), e)
                proj.append([Fe(t, (a_ * xF[k]).u - (b_ * yF[k]).u, (a_ * xF[k]).v - (b_ * yF[k]).v)
                             for k in range(4)])
            nJPm, _ = span_n(t, [EF(a, b) for a in proj for b in LPx])
            nJQ, _ = span_n(t, [EF(a, b) for a in LQx for b in LQx])
            GP = flint.fmpq_mat([[flint.fmpq(int((2 * EF(a, b).u).numerator), int((2 * EF(a, b).u).denominator))
                                  for b in LPx] for a in LPx]).det()
            vcP = vval(Fr(int(GP.p), int(GP.q)), l) // 2
            g = content(T12)
            chi = W / g ** 2
            lemma = vval(iota, l) <= vcP - vval(dF, l) - vval(nJPm, l)
            trace = 2 * vval(g, l) <= vval(nJPm, l) + vval(nJQ, l) + vval(dF, l)
            final = vval(chi, l) >= vval(4 * dF, l)
            good_all &= lemma and trace and final
            tight += vval(iota, l) == vcP - vval(dF, l) - vval(nJPm, l) and vval(iota, l) > 0
        ok &= good_all
        rows.append({"t": t, "N": N, "lattices": ntr, "tight_with_glue": tight})
        print("    (O) t=%d N=%d (l=%d): %d lattices: generator pair found, Lemma 9.9, trace bound and"
              " v_l(chi(E2)) >= v_l(4 d_F): %s (Lemma 9.9 tight with nonzero glue: %d)"
              % (t, N, l, ntr, good_all, tight), flush=True)
    return rows, ok


# ----------------------------------------------------------------------------- dimensions n = 6, 8 (check (P))
def twoform_gen(t, ds, f):
    """Tr(f E) on O_F^n with E = sum_i d_i e_2i^* ^ e_2i+1^* (n = 2 len(ds)); Z-basis (1, omega) per coordinate."""
    bs, out = [Fe(t, 1), omega(t)], {}
    for blk, d in enumerate(ds):
        a, b = 2 * blk, 2 * blk + 1
        for k in range(2):
            for l in range(2):
                val = (f * d * bs[k] * bs[l]).tr()
                i, j = 2 * a + k, 2 * b + l
                if val != 0:
                    out[(1 << i) | (1 << j)] = out.get((1 << i) | (1 << j), 0) + (1 if i < j else -1) * val
    return {k: v for k, v in out.items() if v != 0}


def gmat(form, D):
    M = [[Fr(0)] * D for _ in range(D)]
    for k, v in form.items():
        i, j = [b for b in range(D) if k >> b & 1]
        M[i][j] += v
        M[j][i] -= v
    return M


def gmm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def gform(M):
    D = len(M)
    return {(1 << i) | (1 << j): M[i][j] for i in range(D) for j in range(i + 1, D) if M[i][j]}


def glattice(rng, N, D):
    gens = [[Fr(int(i == j)) for i in range(D)] for j in range(D)]
    if N > 1:
        for _ in range(rng.randint(1, 4)):
            v = [rng.randrange(N) for _ in range(D)]
            if rng.random() < 0.5:
                v = [x if rng.random() < 0.4 else 0 for x in v]
            gens.append([Fr(x, N) for x in v])
    den = 1
    for v in gens:
        for x in v:
            den = den * x.denominator // math.gcd(den, x.denominator)
    H = flint.fmpz_mat([[int(x * den) for x in v] for v in gens]).hnf()
    rows = [[Fr(int(H[r, c]), den) for c in range(D)] for r in range(H.nrows())
            if any(H[r, c] != 0 for c in range(D))]
    return [[rows[c][r] for c in range(D)] for r in range(D)]


def check_P():
    """NOTE Proposition 9.11 on principal and random lattices in dimensions n = 6, 8: for every integral class of B
    with odd part g = y + w sqrt t, a = q^(m-1) P g^2 (P = Pf(Theta|Lambda), m = n/2) is integral on Lambda and
    chi_odd = 2^(n-2) Tr(a) (m even) or -2^(n-1) q^(m-1) P N(g) (m odd); for m odd every value of chi on the
    integral lattices of B is divisible by 2^(n-1); for n = 8 the even classes have chi >= 64 when the rank is
    nonzero, the rank-0 line has chi = 0 mod 4, and no class has chi = 6."""
    rng = random.Random(611)
    rows, ok = [], True
    for n, t, N, ntr in ((6, 5, 1, 1), (6, 5, 2, 3), (6, 57, 3, 2), (6, 13, 13, 2), (8, 11, 1, 1), (8, 5, 1, 1),
                         (8, 11, 2, 1), (8, 2, 1, 1), (8, 3, 1, 1), (8, 5, 5, 2), (8, 3, 3, 2), (8, 2, 2, 2)):
        m, D, dl, gs = n // 2, 2 * n, delta(t), reps(t, 9)[:3]
        good_all, least, n0s = True, None, []
        for tr in range(ntr):
            ds = [dl * rng.choice(gs) for _ in range(m)] if N > 1 else [dl] * m
            A = glattice(rng, N, D)
            AT = [list(r) for r in zip(*A)]
            M1 = gmm(AT, gmm(gmat(twoform_gen(t, ds, Fe(t, 1)), D), A))
            Ms = gmm(AT, gmm(gmat(twoform_gen(t, ds, Fe(t, 0, 1)), D), A))
            L1, Ls = gform(M1), gform(Ms)
            Pw = {}
            for aa in range(n + 1):
                X = {0: Fr(1)}
                for _ in range(aa):
                    X = wmul(X, L1)
                for bb in range(n + 1 - aa):
                    Pw[(aa, bb)] = X
                    X = wmul(X, Ls)
            top = (1 << D) - 1
            s0 = 1 if Pw[(n, 0)].get(top, 0) > 0 else -1
            Pf = abs(Pw[(n, 0)].get(top, 0)) / math.factorial(n)
            if m == 4:
                # the line r = g = 0: chi = 4 int c_4 c_12 >= 4 N0, N0 = int (T/g_4)(T^3/g_12), T = Theta_1 Theta_2
                tt = Fr(t)
                T1 = add((Fr(1, 4), Pw[(2, 0)]), (-1 / (4 * tt), Pw[(0, 2)]))
                T3 = add((Fr(1, 64), Pw[(6, 0)]), (-3 / (64 * tt), Pw[(4, 2)]), (3 / (64 * tt ** 2), Pw[(2, 4)]),
                         (-1 / (64 * tt ** 3), Pw[(0, 6)]))
                T4 = sum(Fr(math.comb(4, j)) * (-1 / tt) ** j * Pw[(8 - 2 * j, 2 * j)].get(top, 0)
                         for j in range(5)) / 256
                n0 = abs(T4) / (content(T1) * content(T3))
                # Theorem 9.14: N0 = 16/(g_4 g_4^v), g_4^v the content of pi_1 pi_2 on the dual lattice
                # (pi = Theta^-1, pi_sqrt t = t L_sqrt t^-1), and 16 | N0
                p1 = gform([list(r) for r in zip(*frmat(qmat(M1).inv()))])
                ps = {k_: tt * v_ for k_, v_ in gform([list(r) for r in zip(*frmat(qmat(Ms).inv()))]).items()}
                Tv = add((Fr(1, 4), wmul(p1, p1)), (-1 / (4 * tt), wmul(ps, ps)))
                good_all &= n0.denominator == 1 and n0 == 16 / (content(T1) * content(Tv)) and n0 % 16 == 0
                n0s.append(int(n0))
            # sqrt t acting on Lambda: Ms = St^T M1 (L_sqrt t = Theta(sqrt t ., .)), so St = M1^-1 Ms
            St = gmm(frmat(qmat(M1).inv()), Ms)
            assert gmm([list(r) for r in zip(*St)], M1) == Ms
            for q in (Fr(1), Fr(1, 3), Fr(3, 2)):
                cl = [{}, {}, {}, {}]
                for k in range(n + 1):
                    c = Fr((-q) ** (k // 2), math.factorial(k))
                    if k % 2 == 0:
                        cl[0] = add((1, cl[0]), (c, Pw[(k, 0)]))
                        cl[2] = add((1, cl[2]), (c / Fr(t) ** (k // 2), Pw[(0, k)]))
                    else:
                        cl[1] = add((1, cl[1]), (c, Pw[(k, 0)]))
                        cl[3] = add((1, cl[3]), (c / Fr(t) ** (k // 2), Pw[(0, k)]))

                def muk(A_, B_):
                    tA = {kk: (-v if (popc(kk) // 2) % 2 else v) for kk, v in A_.items()}
                    return s0 * wmul(tA, B_).get(top, 0)
                gram = [muk(cl[i], cl[i]) for i in range(4)]
                ev, od = lattice2([cl[0], cl[2]]), lattice2([cl[1], cl[3]])
                evv = [(ev[0][0], ev[1][0]), (ev[0][1], ev[1][1])]
                odv = [(od[0][0], od[1][0]), (od[0][1], od[1][1])]
                evv.append((evv[0][0] + evv[1][0], evv[0][1] + evv[1][1]))
                odv.append((odv[0][0] + odv[1][0], odv[0][1] + odv[1][1]))
                for y, w in odv:
                    G = [[y * (i == j) + w * St[i][j] for j in range(D)] for i in range(D)]
                    a_ = [[q ** (m - 1) * Pf * x for x in r] for r in gmm(G, G)]
                    chi = gram[1] * y * y + gram[3] * w * w
                    want = (2 ** (n - 2) * q ** (m - 1) * Pf * 2 * (y * y + t * w * w) if m % 2 == 0 else
                            -2 ** (n - 1) * q ** (m - 1) * Pf * (y * y - t * w * w))
                    good_all &= integral(a_) and chi == want
                    if m % 2:
                        good_all &= chi % 2 ** (n - 1) == 0
                for x, z in evv:
                    chi = gram[0] * x * x + gram[2] * z * z
                    if m % 2:
                        good_all &= chi % 2 ** (n - 1) == 0
                    else:
                        good_all &= (x + z == 0 and chi % 4 == 0) or (x + z != 0 and chi >= 64)
                if m % 2 == 0:
                    small = [v for v, u in enum2(ev, [gram[0], gram[2]], 63)] + \
                        [v for v, u in enum2(od, [gram[1], gram[3]], 63)]
                    good_all &= all(v % 4 == 0 for v in small)
                    vals = [v for v, u in enum2(ev, [gram[0], gram[2]], 300)] + \
                        [v for v, u in enum2(od, [gram[1], gram[3]], 300)]
                    if vals:
                        least = min(vals) if least is None else min(least, min(vals))
        ok &= good_all
        rows.append({"n": n, "t": t, "N": N, "lattices": ntr, "least_chi_le_300": str(least), "N0": n0s})
        print("    (P) n=%d t=%d N=%d: %d lattices, q = 1, 1/3, 3/2: Proposition 9.11 identities and %s: %s"
              % (n, t, N, ntr, "32 | chi" if m % 2 else "chi in 4Z or >= 64, least chi <= 300: %s; line r = g = 0: "
                 "N0 = %s = 16/(g_4 g_4^v), 16 | N0, so chi >= 4 N0 there" % (least, n0s), good_all), flush=True)
    return rows, ok


def check_R():
    """NOTE Proposition 8.1(e): on the principal free lattice (g_1 = g_2 = 1), the nonzero integral classes of B
    with c_1 in Q.Theta (w = 0) and chi <= 14 for q = 1 are exactly +-alpha, +-beta; and the least integer q for
    which alpha~ is integral (the proof needs q >= 2)."""
    rows, ok = [], True
    for t in (5, 2, 3, 13):
        dl = delta(t)
        P, s0 = setup(t, dl, dl)
        found = sorted((int(v), tuple(str(x) for x in u))
                       for v, u, _ in generic_classes(t, 1, P, s0, 14, generic=False) if u[3] == 0)
        want = sorted([(8, ("1", "0", "0", "0")), (8, ("-1", "0", "0", "0")), (8, ("0", "1", "0", "0")),
                       (8, ("0", "-1", "0", "0"))])
        qmin = None
        for q in range(1, 30):
            cl = bclasses(t, q, P)
            if all(Fr(v).denominator == 1 for v in cl[2].values()):
                qmin = q
                break
        good = found == want and qmin is not None and qmin >= 2
        ok &= good
        rows.append({"t": t, "classes": found, "least_q_alpha_tilde_integral": qmin})
        print("    (R) t=%d principal free: integral classes with c_1 in Q Theta and chi <= 14 (q = 1): %s; least "
              "integer q with alpha~ integral: %s: %s" % (t, found, qmin, good), flush=True)
    return rows, ok


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "shape_search.json"
    res, allok, secs = {}, True, {}
    for key, fn in (("H", check_H), ("I", check_I), ("J", check_J), ("L", check_L), ("M", check_M),
                    ("N", check_N), ("O", check_O), ("P", check_P), ("R", check_R)):
        t0 = time.time()
        r, ok = fn()
        secs[key] = round(time.time() - t0, 1)
        res[key] = r
        allok &= ok
        print("(%s) %s" % (key, ok), flush=True)
    with open(out, "w", newline="\n") as f:
        json.dump({"results": res, "seconds": secs, "_all_checks_pass": allok}, f, indent=1, sort_keys=True,
                  default=str)
    print("ALL CHECKS PASS:", allok)


if __name__ == "__main__":
    main()
