"""cm_ranks.py -- the annihilators and ranks of Theorem 2.3 (rm-semiregularity note) for q in F totally positive.

For q in F totally positive put q_k = sigma_k(q) and exp(i sqrt(q_k) Theta_k) = C_k + i sqrt(q_k) S_k, so C_k, S_k are
polynomials in Theta_k with coefficients in F, and B(Theta, q)_R = U_1 U_2, U_k = <C_k, S_k>.  A class of B_R is
    c = M00 C_1C_2 + M01 C_1S_2 + M10 S_1C_2 + M11 S_1S_2,
and it is rational iff M00, M11 in Q, M01 in F and M10 = sigma(M01) (the Hermitian model of cm_fields.bclass).  The
proof of Theorem 2.3 should go through with lambda_+- = +- i sqrt(q_k) in block k.  Checks, in the frame of semireg.py
(Theta_k = Omega_k = sum_{j in block k} dz_j dzbar_j):
 (i)   n = 4, 6, 8: R + G', G' = span{dzbar_a dzbar_b + q_k^-1 d_a d_b : a < b in block k}, annihilates C_1C_2,
       C_1S_2, S_1C_2, S_1S_2 (exact, coefficients in F) and has dimension n^2/2 (rank modulo p); the common kernel
       modulo two primes has dimension n^2/2 and its H^1(T)-part dimension m(m+1).  So K(B) = R + G' and
       K(B) n H^1(T) = R.  (Ranks modulo p are lower bounds for ranks over C, so the modular kernels are upper bounds.)
 (ii)  ranks r_k(c), all k, n = 4 and 6, at random rational points and at rational points lambda u u^sigma of the
       rank-one cone (u in F C_1 + F S_1; also C_1C_2 and S_1S_2): the block formula of Theorem 2.3 (exact ranks of
       the N_{u,v} over Q(sqrt t), sympy) = the direct ranks modulo p = the prediction of Theorem 2.3(iv) with
       (number of nonzero c_ef, tensor rank) = (4, rank M).  The c_ef are computed numerically (floats), so the
       count 4 is a sanity check; for non-biquadratic K it follows from the transitivity of Galois on CM types.
 (iii) biquadratic control q = 3 + sqrt 5 = 2 phi^2: the rational plane spanned by the real part and sqrt(2)^-1 times
       the imaginary part of exp(i sqrt 2 L_phi), L_phi = phi Theta_1 + sigma(phi) Theta_2, is annihilated by
       R + G'; ranks (1,8,12,8,1) and dim K = 16 at a random rational point of it (n = 4).
Fields: t = 5, q = (5 - 2 sqrt 5)/4 (K = Q(zeta_5)); t = 5, q = 4 + sqrt 5; t = 2, q = 1 + sqrt(2)/2; t = 2,
q = 3 + sqrt 2; t = 13, q = (5 + sqrt 13)/2; t = 5, q = 3 + sqrt 5 (biquadratic).
Usage: python -I -B cm_ranks.py      Expected last line "ALL CHECKS PASS: True".
"""
import cmath
import importlib.util
import math
import os
import random
import sys
from fractions import Fraction as Fr

import sympy

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SR = load("semireg", os.path.normpath(os.path.join(HERE, "..", "..", "rm-semiregularity", "code", "semireg.py")))
TL = load("twist_lib", os.path.join(HERE, "twist_lib.py"))
Q2 = TL.Q2

FIELDS = [(5, Q2(5, Fr(5, 4), Fr(-1, 2))), (5, Q2(5, 4, 1)), (2, Q2(2, 1, Fr(1, 2))), (2, Q2(2, 3, 1)),
          (13, Q2(13, Fr(5, 2), Fr(1, 2))), (5, Q2(5, 3, 1))]


def is_sq(x):
    x = Fr(x)
    return x >= 0 and math.isqrt(x.numerator) ** 2 == x.numerator and math.isqrt(x.denominator) ** 2 == x.denominator


def kind(t, q):
    N = q.norm()
    return "biquadratic" if is_sq(N) else ("cyclic" if is_sq(N / t) else "non-Galois")


def cs(n, t, qk):
    """[C, S] for one block: dicts power -> coefficient in F, exp(i sqrt(qk) Theta) = C + i sqrt(qk) S."""
    m = n // 2
    C = {2 * j: Q2(t, Fr(1, math.factorial(2 * j))) * pw(t, -qk, j) for j in range(m // 2 + 1)}
    S = {2 * j + 1: Q2(t, Fr(1, math.factorial(2 * j + 1))) * pw(t, -qk, j) for j in range((m - 1) // 2 + 1)}
    return [C, S]


def pw(t, x, j):
    out = Q2(t, 1)
    for _ in range(j):
        out = out * x
    return out


def tensor_coef(n, t, q, M):
    """coefficients x_ab of Theta_1^a Theta_2^b of sum M[al][be] P1[al] P2[be], P1 = (C_1, S_1), P2 = (C_2, S_2)."""
    m = n // 2
    P1, P2 = cs(n, t, q), cs(n, t, q.conj())
    coef = {(a, b): Q2(t) for a in range(m + 1) for b in range(m + 1)}
    for al in range(2):
        for be in range(2):
            if M[al][be] == 0:
                continue
            for a, ca in P1[al].items():
                for b, cb in P2[be].items():
                    coef[(a, b)] = coef[(a, b)] + M[al][be] * ca * cb
    return coef


def hermitian(coef):
    return all(coef[(b, a)] == coef[(a, b)].conj() for (a, b) in coef)


def ext(n, coef):
    return SR.poly(n, {k: v for k, v in coef.items() if not v == 0})


def as_pairs(vec):
    return {k: (v.u, v.v) if isinstance(v, Q2) else v for k, v in vec.items()}


def kernel_basis(n, t, q):
    """R + G' (G' with coefficient q_k^-1 in block k)."""
    m = n // 2
    qi = [q.inv(), q.conj().inv()]
    out = [{x: Q2(t, c) for x, c in r.items()} for r in SR.rm_basis(n)]
    for k, blk in enumerate((range(m), range(m, n))):
        for a in blk:
            for b in blk:
                if a < b:
                    out.append({("O", a, b): Q2(t, 1), ("L", a, b): qi[k]})
    return out


def rank_combs(n, combs, md):
    cols = [{x: md.red((c.u, c.v)) for x, c in comb.items()} for comb in combs]
    return SR.rank_mod(cols, md.p)


def mods(t, k=2):
    return [SR.Mod(p, t) for p in SR.PRIMES if SR.sqrt_mod(t, p) is not None][:k]


def flat_ranks(n, t, coef):
    m = n // 2
    st = sympy.sqrt(t)
    co = {(a, b): (sympy.Rational(v.u.numerator, v.u.denominator) + sympy.Rational(v.v.numerator, v.v.denominator) * st)
          * math.factorial(a) * math.factorial(b) for (a, b), v in coef.items()}
    return SR.flattening_ranks(m, co)


def predicted(m, n4, tr):
    out = []
    for k in range(4 * m + 1):
        tot = 0
        for u in range(m + 1):
            v = k - u
            if v < 0 or v > m:
                continue
            ends = (u in (0, m)) + (v in (0, m))
            if ends == 0:
                rho = n4
            elif ends == 1:
                rho = 2
            elif (u, v) in ((0, 0), (m, m)):
                rho = 1
            else:
                rho = tr
            tot += math.comb(m, u) * math.comb(m, v) * rho
        out.append(tot)
    return out


def c_ef(t, q, M):
    """the coordinates of c in the basis l_ef = exp(i(e sqrt(q_1) Theta_1 + f sqrt(q_2) Theta_2)), numerically."""
    r1, r2 = cmath.sqrt(q.emb(0)), cmath.sqrt(q.emb(1))
    Mf = [[(M[a][b].emb(0) if isinstance(M[a][b], Q2) else float(M[a][b])) for b in range(2)] for a in range(2)]
    out = {}
    for e in (1, -1):
        for f in (1, -1):
            w1, w2 = [0.5, e / (2j * r1)], [0.5, f / (2j * r2)]
            out[(e, f)] = sum(Mf[a][b] * w1[a] * w2[b] for a in range(2) for b in range(2))
    return out


def rnd_q(rng):
    return Fr(rng.randint(-9, 9) or 1, rng.randint(1, 5))


def rnd_F(rng, t):
    return Q2(t, rnd_q(rng), rnd_q(rng))


def check_i(t, q):
    ok = True
    for n in (4, 6, 8):
        m = n // 2
        basis = [[[Q2(t, int(i == al and j == be)) for j in range(2)] for i in range(2)]
                 for al in range(2) for be in range(2)]
        classes = [ext(n, tensor_coef(n, t, q, M)) for M in basis]
        kb = kernel_basis(n, t, q)
        exact = all(SR.annihilates(n, comb, classes) for comb in kb)
        mdl = mods(t)
        ind = all(rank_combs(n, kb, md) == n * n // 2 for md in mdl)
        up = [SR.kerdim_mod(n, [as_pairs(c) for c in classes], md) for md in mdl]
        upT = [SR.kerdim_mod(n, [as_pairs(c) for c in classes], md, kinds={"T"}) for md in mdl]
        good = exact and ind and all(u == n * n // 2 for u in up) and all(u == m * (m + 1) for u in upT)
        ok &= good
        print("  (i) n=%d: R+G' annihilates B exactly: %s; independent (dim %d): %s; modular dim K(B) %s, "
              "H^1(T)-part %s: %s" % (n, exact, n * n // 2, ind, up, upT, good), flush=True)
    return ok


def check_ranks(t, q, rng):
    ok = True
    for n in (4, 6):
        m = n // 2
        md = mods(t, 1)[0]
        pts = []
        for _ in range(3 if n == 4 else 1):
            c01 = rnd_F(rng, t)
            pts.append(("random", [[Q2(t, rnd_q(rng)), c01], [c01.conj(), Q2(t, rnd_q(rng))]]))
        for _ in range(2 if n == 4 else 1):
            a, b, lam = rnd_F(rng, t), rnd_F(rng, t), rnd_q(rng)
            pts.append(("cone", [[Q2(t, lam) * a * a.conj(), Q2(t, lam) * a * b.conj()],
                                 [Q2(t, lam) * b * a.conj(), Q2(t, lam) * b * b.conj()]]))
        if n == 4:
            one, zero = Q2(t, 1), Q2(t)
            pts.append(("cone C1C2", [[one, zero], [zero, zero]]))
            pts.append(("cone S1S2", [[zero, zero], [zero, one]]))
        kb = kernel_basis(n, t, q)
        for name, M in pts:
            coef = tensor_coef(n, t, q, M)
            herm = hermitian(coef)
            det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
            tr = 1 if det == 0 else 2
            ce = c_ef(t, q, M)
            n4 = sum(1 for v in ce.values() if abs(v) > 1e-9)
            vec = ext(n, coef)
            inB = all(not SR.act_comb(comb, vec, n) for comb in kb)
            fl = flat_ranks(n, t, coef)
            di = SR.direct_ranks(n, as_pairs(vec), md)
            pr = predicted(m, n4, tr)
            kd = SR.kerdim_mod(n, [as_pairs(vec)], md)
            want_kd = n * (2 * n - 1) - pr[2]
            good = herm and inB and n4 == 4 and (name == "random") == (tr == 2) and fl == di == pr and kd == want_kd
            ok &= good
            print("  (ii) n=%d %-10s rational %s, K(B) kills it %s, nonzero c_ef %d, tensor rank %d: ranks %s "
                  "(block formula = direct mod p = predicted: %s), dim K(c) mod p %d: %s"
                  % (n, name, herm, inB, n4, tr, di[:n + 1], fl == di == pr, kd, good), flush=True)
    return ok


def check_biquadratic(rng):
    t, q, n, m = 5, Q2(5, 3, 1), 4, 2
    phi = Q2(t, Fr(1, 2), Fr(1, 2))
    assert q == Q2(t, 2) * phi * phi
    # exp(i sqrt 2 L) = Re + i sqrt 2 Im', L = phi Theta_1 + sigma(phi) Theta_2: Re = sum (-2)^j L^2j/(2j)!,
    # Im' = sum (-2)^j L^(2j+1)/(2j+1)!
    def Lpow(k):
        return {(a, k - a): Q2(t, math.comb(k, a)) * pw(t, phi, a) * pw(t, phi.conj(), k - a)
                for a in range(k + 1) if a <= m and k - a <= m}
    x, y = rnd_q(rng), rnd_q(rng)
    coef = {(a, b): Q2(t) for a in range(m + 1) for b in range(m + 1)}
    for k in range(2 * m + 1):
        f = Fr((-2) ** (k // 2), math.factorial(k)) * (x if k % 2 == 0 else y)
        for key, v in Lpow(k).items():
            coef[key] = coef[key] + Q2(t, f) * v
    herm = hermitian(coef)
    vec = ext(n, coef)
    kb = kernel_basis(n, t, q)
    inB = all(not SR.act_comb(comb, vec, n) for comb in kb)
    md = mods(t, 1)[0]
    di = SR.direct_ranks(n, as_pairs(vec), md)
    fl = flat_ranks(n, t, coef)
    kd = SR.kerdim_mod(n, [as_pairs(vec)], md)
    good = herm and inB and fl == di and di[:n + 1] == [1, 8, 12, 8, 1] and kd == 16
    print("  (iii) biquadratic q = 3 + sqrt 5 = 2 phi^2, point of the rational plane: rational %s, K(B) kills it %s, "
          "ranks %s (block formula agrees: %s), dim K(c) mod p %d: %s" % (herm, inB, di[:n + 1], fl == di, kd, good),
          flush=True)
    return good


def main():
    rng = random.Random(20261025)
    ok = True
    for t, q in FIELDS:
        kd = kind(t, q)
        print("t = %d, q = %s (sigma_1 q = %.4f, sigma_2 q = %.4f, N(q) = %s, K %s)"
              % (t, q, q.emb(0), q.emb(1), q.norm(), kd), flush=True)
        assert q.totally_positive()
        ok &= check_i(t, q)
        if kd != "biquadratic":
            ok &= check_ranks(t, q, rng)
    ok &= check_biquadratic(rng)
    print("ALL CHECKS PASS:", ok)


if __name__ == "__main__":
    main()
