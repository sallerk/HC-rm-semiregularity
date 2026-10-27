"""
cone6.py -- checks for ../NOTE.txt, Theorems 9.15, 9.16, Proposition 9.17 and Example 9.18 (the rank-one cone of
B in dimension 6, OPEN 11(l)) and Proposition 10.7 (sums of two twisted divisor sheaves, OPEN 11(f)).

Setting as in semireg.py: frame of NOTE Theorem 3.1 (Theta_k = sum over block k of dz_j dzbar_j after rescaling),
B = U_1 (x) U_2, U_k = <C_k, S_k>, exp(i sqrt q Theta_k) = C_k + i sqrt q S_k; for n = 6 (blocks of size 3)
C_k = 1 - (q/2) Theta_k^2 and S_k = Theta_k - (q/6) Theta_k^3.  A class of B is c = x alpha + y beta + z alpha~ +
w sqrt t beta~ = r C_1 C_2 + q v S_1 S_2 + g S_1 C_2 + gbar C_1 S_2 (r = x + z, v = z - x, g = y + w sqrt t); the
rank-one cone is q r v = N(g).

 (U) NOTE Theorem 9.15 (n = 6: on the cone, (M_a) + (M_b) force the all-degree mechanism).  (i) Ranks of
     contraction with c on HT^k = wedge^k W, all k: (1, 12, 48, 74, 48, 12, 1) at six cone points (block formula of
     (G)), and by direct modular ranks in every degree at C_1 C_2 and S_1 S_2 (q = 3).  (ii) The count: the only
     integers e_1 >= 12, e_3 >= 74 with 2 - 2 e_1 + 2 * 48 - e_3 = 0 are (12, 74); the profile is the coefficient
     list of (1 + 6x + 6x^2 + x^3)^2.  (iii) One block (a complex 3-torus): the ranks of contraction with a C + b S
     are (1, 6, 6, 1) (direct, several (a, b) and q).  (iv) The ideal ann(u_1) (x) wedge W_2 + wedge W_1 (x)
     ann(u_2) annihilates c = u_1 u_2: every kernel element of wedge^2 W_1 and wedge^3 W_1 on u_1 (modular
     basis), times every monomial of wedge W_2, kills c; likewise with the blocks exchanged.  (v) dim K(c) = 18 =
     n^2/2 (modular upper bound, equal to 66 - r_2) at three rational cone points for t = 57, q = 3.
 (V) NOTE Theorem 9.16 (the split model).  (i) Markman's identity e^Theta (1 - (N/2) Theta^2 + (N/3) Theta^3) =
     alpha_0 + beta_0 with q = N - 1 on a threefold (N symbolic): ch(I_{C_1 u ... u C_N}) = e^-Theta (alpha_0 +
     beta_0).  (ii) The class ch(E') = f(Theta_1) f(Theta_2), f = 1 - (N/2) Theta^2 + (N/3) Theta^3 (E' = I (x) I
     on X_0 x X_0, N = 2, 4): ranks (1, 12, 48, 74, 48, 12, 1) (block formula; direct modular ranks for N = 4) and
     Mukai square 0.  (iii) On one block, A = e_Theta - (1/q) Lambda (Lambda = sum_j i(dzbar_j^*) i(dz_j^*)) maps
     C -> 3 S and S -> -(3/q) C (a rotation of U, exact, q symbolic), A Omega = 0, A is skew for the Mukai pairing
     of the block, and [A, w] and [e_Theta, w] lie in W for the six generators w of W (exact, q = 3).
 (W) NOTE Proposition 9.17 and Example 9.18 (integral points of the cone).  (i) Symbolic: (r, q v, g) =
     lambda (N(a), N(b), b sigma(a)) for c = lambda u sigma(u), u = a C_1 + b S_1, and the cone equation;
     Riemann-Roch on the cone, int c e^{s_1 Theta_1 + s_2 Theta_2} = P f_1(s_1) f_2(s_2) with f_k(s) =
     a_k (s^3 - 3 q s) + b_k (3 s^2 - q); the Bogomolov discriminant on the cone, c_1^2 - 2 r ch_2 =
     sum_k (sigma_k(g)^2 + q r^2) Theta_k^2; and chi = 32 (r chi(F) - nu) with chi(F) = q^3 P v, nu = q^2 P N(g).
     (ii) Exact enumeration (Gauss-reduced, all classes below a height bound, height h = 32 P (q^3 (x^2 + z^2) +
     q^2 (y^2 + t w^2))) of the integral cone points on the principal free lattices (t, q) = (57, 3), (57, 1),
     (5, 1), (5, 3), (13, 1), (2, 1), (3, 1) and on two non-free lattices for (57, 3); for the smallest classes:
     integrality of every component, chi = 0, chi(F) = int c = q^3 P v, nu = q^2 P N(g) = r chi(F), a = q^2 P g^2
     integral on Lambda, and g = 0 when r = 0; the gcd of the ranks that occur below the bound.
 (X) NOTE Proposition 10.7 (n = 4).  c = ch O_D1(M_1) + ch O_D2(M_2) = e^{L_m1} - e^{L_{m1 - g1}} + e^{L_m2} -
     e^{L_{m2 - g2}}.  (i) Symbolic identities (t symbolic): the B-conditions with sum eps = 0, the formulas for
     sum eps h^2, sum eps N(h), sum eps N(h)^2, Hn = sum eps N(h) h, q = -Hn G/N(G), the coordinates and chi =
     W [(sum eps N(h))^2 + q Tr(G^2)]; the parametrization p_1 = g_2 s, p_2 = -g_1 s; the factorizations in the
     cases s rational and g_2 = lambda g_1.  (ii) The content of Theta_1 Theta_2 on the principal free lattice is
     1/d_F (t = 1 mod 4) or 2/d_F, for 9 fields.  (iii) The remaining case (g_1 = g_2 irrational, Tr(a) <= 11): for
     every real quadratic field the finite set of a, and the finitely many candidate scalings from the
     factorization of 3 (r_0^2 - 1) L^2; no solution; as a control, the configuration of NOTE 10.5(d) is recovered
     from its a.  (iv) Brute force (t = 5, 13, 2, boxes): every sum of two twisted divisor sheaves found in B lies
     in P.  (v) The classes of NOTE 10.5(c), (d): chi = 32, 89400928 and 6592.
Usage:  python3 -I cone6.py OUT.json      Expected last line "ALL CHECKS PASS: True".  About 20 seconds.
"""

import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint  # noqa: E402
import sympy  # noqa: E402

import semireg as SR  # noqa: E402
import shape_search as SS  # noqa: E402

PROFILE6 = [1, 12, 48, 74, 48, 12, 1]


def coefs(m, e, T1, T2):
    P_ = sympy.Poly(sympy.expand(e), T1, T2)
    return {(i, j): P_.coeff_monomial(T1 ** i * T2 ** j) * math.factorial(i) * math.factorial(j)
            for i in range(m + 1) for j in range(m + 1)}


def cs(m, q, X):
    return (sum((-q) ** (k // 2) * X ** k / sympy.factorial(k) for k in range(0, m + 1, 2)),
            sum((-q) ** (k // 2) * X ** k / sympy.factorial(k) for k in range(1, m + 1, 2)))


def frame_class(n, e, T1, T2):
    """the polynomial e(T1, T2) (rational coefficients) as an element of the exterior algebra (Omega frame)."""
    P_ = sympy.Poly(sympy.expand(e), T1, T2)
    co = {(i, j): Fr(int(sympy.numer(c)), int(sympy.denom(c))) for (i, j), c in P_.terms()}
    return SR.poly(n, co)


def block_class(nn, f):
    """sum_k f[k] Theta^k on a single block of nn coordinates (Theta = sum_j dz_j dzbar_j)."""
    vec = {}
    for S_ in range(1 << nn):
        k = SR.popc(S_)
        fk = f.get(k, 0)
        if fk == 0:
            continue
        cur = {0: 1}
        for l in range(nn):
            if S_ >> l & 1:
                cur = SR.wedge(nn + l, SR.wedge(l, cur))
        for mk, cc in cur.items():
            vec[mk] = vec.get(mk, 0) + fk * math.factorial(k) * cc
    return {k: v for k, v in vec.items() if v != 0}


# ----------------------------------------------------------------------------- (U)
def apply_mono(ops, vec):
    """apply the product of the operators ops (list, rightmost first applied last as written) to vec."""
    for kind, a in reversed(ops):
        vec = SR.wedge(a, vec) if kind == "w" else SR.contr(a, vec)
        if not vec:
            break
    return vec


def check_U():
    res, ok = {}, True
    T1, T2 = sympy.symbols("T1 T2")
    # (i) cone points, n = 6
    st = sympy.sqrt(57)
    pts = []
    for qq in (3, sympy.Rational(1, 2)):
        (C1, S1), (C2, S2) = cs(3, qq, T1), cs(3, qq, T2)
        pts += [("C1C2 q=%s" % qq, C1 * C2), ("S1S2 q=%s" % qq, S1 * S2)]
    (C1, S1), (C2, S2) = cs(3, 3, T1), cs(3, 3, T2)
    pts += [("(C1+15S1)(C2+15S2) q=3", (C1 + 15 * S1) * (C2 + 15 * S2)),
            ("(C1+(4+19w)S1)(C2+(4+19wbar)S2) q=3, t=57",
             (C1 + (4 + 19 * (1 + st) / 2) * S1) * (C2 + (4 + 19 * (1 - st) / 2) * S2))]
    block = {nm: SR.flattening_ranks(3, coefs(3, e, T1, T2))[:7] for nm, e in pts}
    good1 = all(v == PROFILE6 for v in block.values())
    md = next(SR.Mod(p, 57) for p in SR.PRIMES if SR.sqrt_mod(57, p) is not None)
    direct = {}
    for nm, e in (("C1C2 q=3", C1 * C2), ("S1S2 q=3", S1 * S2)):
        direct[nm] = SR.direct_ranks(6, frame_class(6, e, T1, T2), md)
    good1 &= all(v == PROFILE6 + [0] * 6 for v in direct.values())
    ok &= good1
    res["ranks"] = {"block": block, "direct": direct}
    print("    (U) n=6 cone points, ranks on HT^k (block formula) %s; direct modular, all degrees: %s: %s"
          % (block, direct, good1), flush=True)
    # (ii) the count
    sols = [(e1, 98 - 2 * e1) for e1 in range(12, 300) if 98 - 2 * e1 >= 74]
    x = sympy.symbols("x")
    sq = sympy.Poly(sympy.expand((1 + 6 * x + 6 * x ** 2 + x ** 3) ** 2), x).all_coeffs()[::-1]
    alt = sum((-1) ** k * r for k, r in enumerate(PROFILE6))
    good2 = sols == [(12, 74)] and [int(c) for c in sq] == PROFILE6 and alt == 0 \
        and [math.comb(12, k) for k in range(7)][2] == 66 and 66 - PROFILE6[2] == 18
    ok &= good2
    res["count"] = {"solutions_e1_e3": sols, "square": [int(c) for c in sq], "alternating": alt}
    print("    (U) count: e_1 >= 12, e_3 >= 74, 2 - 2 e_1 + 96 - e_3 = 0 -> %s; profile = (1 + 6x + 6x^2 + x^3)^2: "
          "%s" % (sols, good2), flush=True)
    # (iii) one block: ranks (1, 6, 6, 1) for a C + b S
    md3 = SR.Mod(SR.PRIMES[0], 1)
    blk = {}
    for q in (Fr(3), Fr(1, 2), Fr(7)):
        for a, b in ((1, 0), (0, 1), (1, 2), (3, -1), (2, 5)):
            f = {0: Fr(a), 1: Fr(b), 2: -q * a / 2, 3: -q * b / 6}
            blk[(str(q), a, b)] = SR.direct_ranks(3, block_class(3, f), md3)
    good3 = all(v == [1, 6, 6, 1, 0, 0, 0] for v in blk.values())
    ok &= good3
    res["one_block"] = {str(k): v for k, v in blk.items()}
    print("    (U) one block, ranks of contraction with a C + b S (q = 3, 1/2, 7; five (a, b)): all (1, 6, 6, 1): "
          "%s" % good3, flush=True)
    # (iv) the tensor ideal annihilates c = u_1 u_2 (n = 6, q = 3, u_1 = C_1 + 2 S_1, u_2 = 3 C_2 - S_2)
    q = Fr(3)
    n = 6
    u1 = SR.poly(n, {(0, 0): Fr(1), (1, 0): Fr(2), (2, 0): -q / 2, (3, 0): -q * 2 / 6})
    u2 = SR.poly(n, {(0, 0): Fr(3), (0, 1): Fr(-1), (0, 2): -q * 3 / 2, (0, 3): q / 6})
    c = frame_class(6, (C1 + 2 * S1) * (3 * C2 - S2), T1, T2)
    assert c == {k: v for k, v in SR.poly(n, {(i, j): a_ * b_ for (i, a_) in
                                                  ((0, 1), (1, 2), (2, -q / 2), (3, -q / 3))
                                                  for (j, b_) in ((0, 3), (1, -1), (2, -3 * q / 2), (3, q / 6))
                                                  }).items() if v != 0}
    p = md3.p
    gens = {1: [("w", n + j) for j in range(3)] + [("c", j) for j in range(3)],
            2: [("w", n + j) for j in range(3, 6)] + [("c", j) for j in range(3, 6)]}
    tested, killed = 0, True
    for blk_i, (ua, other) in ((1, (u1, 2)), (2, (u2, 1))):
        G = gens[blk_i]
        for deg in (2, 3):
            monos = list(itertools.combinations(G, deg))
            cols = []
            for mo in monos:
                v = apply_mono(list(mo), {k: md3.red(v) for k, v in ua.items()})
                cols.append({kk: vv % p for kk, vv in v.items() if vv % p})
            keys = sorted(set(kk for col in cols for kk in col))
            M = flint.nmod_mat([[col.get(kk, 0) for col in cols] for kk in keys] or [[0] * len(cols)], p)
            K, nul = M.nullspace()
            kern = [[int(K[i, j]) for i in range(len(monos))] for j in range(nul)]
            for kv in kern:
                for r in range(0, 7):
                    for mo2 in itertools.combinations(gens[other], r):
                        v2 = apply_mono(list(mo2), {k: md3.red(v) for k, v in c.items()})
                        tot = {}
                        for coef, mo in zip(kv, monos):
                            if coef:
                                for kk, vv in apply_mono(list(mo), v2).items():
                                    tot[kk] = (tot.get(kk, 0) + coef * vv) % p
                        tested += 1
                        if any(vv % p for vv in tot.values()):
                            killed = False
            res.setdefault("kernel_dims", {})["block%d_deg%d" % (blk_i, deg)] = len(kern)
    good4 = killed and res["kernel_dims"] == {"block1_deg2": 9, "block1_deg3": 19, "block2_deg2": 9,
                                              "block2_deg3": 19}
    ok &= good4
    print("    (U) ann(u_1) (x) wedge W_2 + wedge W_1 (x) ann(u_2) kills c = u_1 u_2: kernel dims %s, %d products "
          "tested: %s" % (res["kernel_dims"], tested, good4), flush=True)
    # (v) dim K(c) = 18 at three rational cone points, t = 57, q = 3
    al, be, alt_, bet = SR.b_basis(6, Fr(3))
    kd = {}
    # (C_1 + 15 S_1)(C_2 + 15 S_2) = (x, y, z, w) = (-37, 15, 38, 0); 57 S_1 S_2 = (-19/2, 0, 19/2, 0)
    # and the irrational class (C_1 + g S_1)(C_2 + gbar S_2), g = 4 + 19 omega = 27/2 + (19/2) sqrt 57:
    # r = 1, q v = N(g) = -4962, (x, z) = ((1 - v)/2, (1 + v)/2), (y, w) = (27/2, 19/2)
    v_ = Fr(-4962, 3)
    for nm, co in (("(C1+15S1)(C2+15S2)", (-37, 15, 38, 0)), ("57 S1S2", (Fr(-19, 2), 0, Fr(19, 2), 0)),
                   ("(C1+gS1)(C2+gbarS2), g = 4 + 19 omega", ((1 - v_) / 2, Fr(27, 2), (1 + v_) / 2, Fr(19, 2)))):
        keys = set(al) | set(be) | set(alt_) | set(bet)
        cls = {k: (co[0] * al.get(k, 0) + co[1] * be.get(k, 0) + co[2] * alt_.get(k, 0), co[3] * bet.get(k, 0))
               for k in keys}
        kd[nm] = SR.kerdim_mod(6, [cls], md)
    good5 = all(v == 18 for v in kd.values())
    ok &= good5
    res["dim_K"] = kd
    print("    (U) dim K(c) at rational cone points (t = 57, q = 3): %s (= 18 = n^2/2): %s" % (kd, good5),
          flush=True)
    return res, ok


# ----------------------------------------------------------------------------- (V)
def check_V():
    res, ok = {}, True
    T, N_, T1, T2 = sympy.symbols("T N T1 T2")

    def trunc(e, X, m):
        P_ = sympy.Poly(sympy.expand(e), X)
        return sum(c * X ** k for (k,), c in P_.terms() if k <= m)
    f = 1 - N_ / 2 * T ** 2 + N_ / 3 * T ** 3
    lhs = trunc(sympy.exp(T).series(T, 0, 4).removeO() * f, T, 3)
    qd = N_ - 1
    rhs = 1 + T - qd / 2 * T ** 2 - qd / 6 * T ** 3
    good1 = sympy.expand(lhs - rhs) == 0
    ok &= good1
    print("    (V) e^Theta (1 - (N/2) Theta^2 + (N/3) Theta^3) = alpha_0 + beta_0 with q = N - 1 (threefold): %s"
          % good1, flush=True)
    # (ii) ranks and Mukai square of ch(E') = f(T1) f(T2)
    md = SR.Mod(SR.PRIMES[0], 1)
    rk = {}
    for NN in (2, 4):
        fe = (1 - sympy.Rational(NN, 2) * T1 ** 2 + sympy.Rational(NN, 3) * T1 ** 3) * \
             (1 - sympy.Rational(NN, 2) * T2 ** 2 + sympy.Rational(NN, 3) * T2 ** 3)
        rk[NN] = SR.flattening_ranks(3, coefs(3, fe, T1, T2))[:7]
        if NN == 4:
            rk["direct N=4"] = SR.direct_ranks(6, frame_class(6, fe, T1, T2), md)[:7]
        # Mukai square, int T1^3 T2^3 = 36 in the frame
        P_ = sympy.Poly(sympy.expand(fe), T1, T2)
        tau = sum(cc * (-1) ** (i + j) * T1 ** i * T2 ** j for (i, j), cc in P_.terms())
        prod = sympy.Poly(sympy.expand(tau * fe), T1, T2).coeff_monomial(T1 ** 3 * T2 ** 3) * 36
        rk["mukai N=%d" % NN] = int(prod)
    good2 = rk[2] == PROFILE6 and rk[4] == PROFILE6 and rk["direct N=4"] == PROFILE6 and \
        rk["mukai N=2"] == 0 and rk["mukai N=4"] == 0
    ok &= good2
    res["model_ranks"] = {str(k): v for k, v in rk.items()}
    print("    (V) ch(E') = f(Theta_1) f(Theta_2), N = 2, 4: ranks %s, Mukai square 0: %s" % (rk, good2),
          flush=True)
    # (iii) the rotation on one block
    n = 3
    masks = list(range(1 << (2 * n)))
    idx = {m_: i for i, m_ in enumerate(masks)}
    q = sympy.symbols("q", positive=True)

    def opmat(fn):
        M = sympy.zeros(len(masks), len(masks))
        for m_ in masks:
            for k, v in fn({m_: 1}).items():
                M[idx[k], idx[m_]] += v
        return M
    Wg = [opmat(lambda v, j=j: SR.wedge(n + j, v)) for j in range(n)] + \
         [opmat(lambda v, j=j: SR.contr(j, v)) for j in range(n)]
    E = sympy.zeros(64, 64)
    Lm = sympy.zeros(64, 64)
    for j in range(n):
        E += opmat(lambda v, j=j: SR.wedge(j, SR.wedge(n + j, v)))
        Lm += opmat(lambda v, j=j: SR.contr(n + j, SR.contr(j, v)))
    A = E - Lm / q

    def vec(cl):
        xx = sympy.zeros(64, 1)
        for k, v in cl.items():
            xx[idx[k]] += v
        return xx
    Tp = [vec({0: 1})]
    for k in range(3):
        Tp.append(E * Tp[-1])          # Theta^k with Theta = E(1) = sum_j dz_j dzbar_j
    C, S = Tp[0] - q / 2 * Tp[2], Tp[1] - q / 6 * Tp[3]
    rot = sympy.simplify(A * C - 3 * S) == sympy.zeros(64, 1) and \
        sympy.simplify(A * S + 3 / q * C) == sympy.zeros(64, 1)
    Om = vec({(1 << n) - 1: 1})
    kill = A * Om == sympy.zeros(64, 1)
    # Mukai pairing of the block: <u, v> = top coefficient of tau(u) ^ v, tau = (-1)^(j(j-1)/2) on degree j
    top = (1 << (2 * n)) - 1
    Tau = sympy.diag(*[(-1) ** (SR.popc(m_) * (SR.popc(m_) - 1) // 2) for m_ in masks])

    def pairing_matrix():
        Pm = sympy.zeros(64, 64)
        for a_ in masks:
            for b_ in masks:
                if a_ & b_ or (a_ | b_) != top:
                    continue
                s_, xx = 0, b_
                while xx:
                    jj = (xx & -xx).bit_length() - 1
                    s_ += SR.popc(a_ >> (jj + 1))
                    xx &= xx - 1
                Pm[idx[a_], idx[b_]] = -1 if s_ & 1 else 1
        return Pm
    Pm = pairing_matrix()
    G = Tau.T * Pm          # <u, v> = u^T G v
    skew = sympy.simplify(A.T * G + G * A) == sympy.zeros(64, 64)
    Aq = A.subs(q, 3)
    Wv = sympy.Matrix.hstack(*[w.reshape(64 * 64, 1) for w in Wg])
    rW = Wv.rank()
    normA = all(Wv.row_join((Aq * w - w * Aq).reshape(64 * 64, 1)).rank() == rW for w in Wg)
    normE = all(Wv.row_join((E * w - w * E).reshape(64 * 64, 1)).rank() == rW for w in Wg)
    good3 = rot and kill and skew and normA and normE
    ok &= good3
    res["rotation"] = {"A C = 3 S, A S = -(3/q) C": rot, "A Omega = 0": kill, "skew": skew,
                       "[A, W] in W": normA, "[e_Theta, W] in W": normE}
    print("    (V) one block, A = e_Theta - (1/q) Lambda: %s: %s" % (res["rotation"], good3), flush=True)
    return res, ok


# ----------------------------------------------------------------------------- (W)
def lattice6(t, ds, A=None):
    n, D = 6, 12
    if A is None:
        A = [[Fr(int(i == j)) for i in range(D)] for j in range(D)]
    AT = [list(r) for r in zip(*A)]
    M1 = SS.gmm(AT, SS.gmm(SS.gmat(SS.twoform_gen(t, ds, SS.Fe(t, 1)), D), A))
    Ms = SS.gmm(AT, SS.gmm(SS.gmat(SS.twoform_gen(t, ds, SS.Fe(t, 0, 1)), D), A))
    L1, Ls = SS.gform(M1), SS.gform(Ms)
    Pw = {}
    for aa in range(n + 1):
        X = {0: Fr(1)}
        for _ in range(aa):
            X = SS.wmul(X, L1)
        for bb in range(n + 1 - aa):
            Pw[(aa, bb)] = X
            X = SS.wmul(X, Ls)
    top = (1 << D) - 1
    s0 = 1 if Pw[(n, 0)].get(top, 0) > 0 else -1
    Pf = abs(Pw[(n, 0)].get(top, 0)) / math.factorial(n)
    St = SS.gmm(SS.frmat(SS.qmat(M1).inv()), Ms)
    return Pw, s0, Pf, top, St


def bclasses6(t, q, Pw, n=6):
    cl = [{}, {}, {}, {}]
    for k in range(n + 1):
        c = Fr((-q) ** (k // 2), math.factorial(k))
        if k % 2 == 0:
            cl[0] = SS.add((1, cl[0]), (c, Pw[(k, 0)]))
            cl[2] = SS.add((1, cl[2]), (c / Fr(t) ** (k // 2), Pw[(0, k)]))
        else:
            cl[1] = SS.add((1, cl[1]), (c, Pw[(k, 0)]))
            cl[3] = SS.add((1, cl[3]), (c / Fr(t) ** (k // 2), Pw[(0, k)]))
    return cl


def cone_points(t, q, ev, od, Pf, H):
    """all integral cone points (x, y, z, w) != 0 with height <= H (exact: Gauss-reduced enumeration)."""
    he = [32 * Pf * q ** 3, 32 * Pf * q ** 3]
    ho = [32 * Pf * q ** 2, 32 * Pf * q ** 2 * t]
    E_ = SS.enum2(ev, he, H) + [(Fr(0), (Fr(0), Fr(0)))]
    O_ = SS.enum2(od, ho, H) + [(Fr(0), (Fr(0), Fr(0)))]
    byN = {}
    for hv, (y, w) in O_:
        byN.setdefault(y * y - t * w * w, []).append((hv, y, w))
    out = []
    for hv, (x, z) in E_:
        for hv2, y, w in byN.get(q * (z * z - x * x), []):
            if hv + hv2 == 0 or hv + hv2 > H:
                continue
            out.append((hv + hv2, (x, y, z, w)))
    out.sort()
    return out


def canon(co):
    """representative modulo c -> -c and tau: (x, y, z, w) -> (x, -y, z, -w)."""
    x, y, z, w = co
    reps = [(x, y, z, w), (-x, -y, -z, -w), (x, -y, z, -w), (-x, y, -z, w)]
    return max(reps, key=lambda r: (r[0] + r[2], r[1], r[3], r[2]))


def check_W():
    res, ok = {}, True
    # (i) symbolic identities
    a0, a1, b0, b1, lam, t, q, s1, s2 = sympy.symbols("a0 a1 b0 b1 lam t q s1 s2")
    rt = sympy.sqrt(t)
    a, abar, b, bbar = a0 + a1 * rt, a0 - a1 * rt, b0 + b1 * rt, b0 - b1 * rt
    r_, qv, g, gbar = lam * a * abar, lam * b * bbar, lam * b * abar, lam * bbar * a
    cone = sympy.expand(r_ * qv - g * gbar) == 0
    T1, T2 = sympy.symbols("T1 T2")
    A1, B1, A2, B2 = sympy.symbols("A1 B1 A2 B2")
    (C1, S1), (C2, S2) = cs(3, q, T1), cs(3, q, T2)
    c = (A1 * C1 + B1 * S1) * (A2 * C2 + B2 * S2)

    def integ(e):     # int over the frame: T1^3 T2^3 -> 36 = (3!)^2 P with P = 1
        P_ = sympy.Poly(sympy.expand(e), T1, T2)
        return P_.coeff_monomial(T1 ** 3 * T2 ** 3) * 36

    ex = sum((s1 * T1 + s2 * T2) ** k / sympy.factorial(k) for k in range(7))
    fA = lambda A_, B_, s_: A_ * (s_ ** 3 - 3 * q * s_) + B_ * (3 * s_ ** 2 - q)
    rr = sympy.expand(integ(c * ex) - fA(A1, B1, s1) * fA(A2, B2, s2)) == 0
    # Bogomolov on the cone: c = r C1C2 + q v S1S2 + g S1C2 + gb C1S2 with q r v = g gb
    rr_, vv, g1, g2 = sympy.symbols("r v g1 g2")
    cc = rr_ * C1 * C2 + q * vv * S1 * S2 + g1 * S1 * C2 + g2 * C1 * S2
    P_ = sympy.Poly(sympy.expand(cc), T1, T2)
    ch1 = P_.coeff_monomial(T1) * T1 + P_.coeff_monomial(T2) * T2
    ch2 = sum(P_.coeff_monomial(T1 ** i * T2 ** (2 - i)) * T1 ** i * T2 ** (2 - i) for i in range(3))
    disc = sympy.expand(ch1 ** 2 - 2 * rr_ * ch2)
    want = (g1 ** 2 + q * rr_ ** 2) * T1 ** 2 + (g2 ** 2 + q * rr_ ** 2) * T2 ** 2
    bog = sympy.expand((disc - want).subs(g2, q * rr_ * vv / g1)) == 0
    # chi = 32 (r chi(F) - nu), chi(F) = q^3 v, nu = q^2 g1 g2 (P = 1)
    tau = sum(cc_ * (-1) ** (i + j) * T1 ** i * T2 ** j for (i, j), cc_ in P_.terms())
    chi = integ(tau * cc)
    chiF = integ(cc)
    chiok = sympy.expand(chi - 32 * (rr_ * q ** 3 * vv - q ** 2 * g1 * g2)) == 0 and \
        sympy.expand(chiF - q ** 3 * vv) == 0
    good1 = cone and rr and bog and chiok
    ok &= good1
    res["identities"] = {"cone": cone, "riemann_roch": rr, "bogomolov": bog, "chi": chiok}
    print("    (W) identities: rational cone points lambda u sigma(u) lie on the cone %s; int c e^{s.Theta} = "
          "P f_1(s_1) f_2(s_2) %s; c_1^2 - 2 r ch_2 = sum (g_k^2 + q r^2) Theta_k^2 on the cone %s; chi = "
          "32 (r chi(F) - nu), chi(F) = q^3 P v %s: %s" % (cone, rr, bog, chiok, good1), flush=True)
    # (ii) enumeration
    rng = random.Random(9152)
    cases = [(57, Fr(3), "principal", None, 3000000), (57, Fr(1), "principal", None, 300000),
             (5, Fr(1), "principal", None, 2000), (5, Fr(3), "principal", None, 20000),
             (13, Fr(1), "principal", None, 12000), (2, Fr(1), "principal", None, 2000),
             (3, Fr(1), "principal", None, 5000)]
    for N in (2, 3):
        cases.append((57, Fr(3), "N=%d" % N, SS.glattice(rng, N, 12), 2600000))
    tab = []
    for t_, q_, kind, Alat, H in cases:
        t0 = time.time()
        Pw, s0, Pf, top, St = lattice6(t_, [SS.delta(t_)] * 3, Alat)
        cl = bclasses6(t_, q_, Pw)
        ev, od = SS.lattice2([cl[0], cl[2]]), SS.lattice2([cl[1], cl[3]])
        # the primitive integral class on the line r = g = 0 (x = -z): make sure the bound reaches it
        b1, b2 = (ev[0][0], ev[1][0]), (ev[0][1], ev[1][1])
        k1, k2 = b2[0] + b2[1], -(b1[0] + b1[1])
        gk = math.gcd(int(k1 * k1.denominator * k2.denominator), int(k2 * k1.denominator * k2.denominator))
        k1, k2 = k1 * k1.denominator * k2.denominator / gk, k2 * k1.denominator * k2.denominator / gk
        xl = k1 * b1[0] + k2 * b2[0]
        hline = 32 * Pf * q_ ** 3 * 2 * xl * xl
        H = max(H, int(hline) + 1)
        pts = cone_points(t_, q_, ev, od, Pf, H)

        def muk(A_, B_):
            tA = {kk: (-v if (SS.popc(kk) // 2) % 2 else v) for kk, v in A_.items()}
            return s0 * SS.wmul(tA, B_).get(top, 0)
        reps, seen = [], set()
        for hv, co in pts:
            cn = canon(co)
            if cn not in seen:
                seen.add(cn)
                reps.append((hv, cn))
        good, rows = True, []
        ranks = sorted(set(abs(co[0] + co[2]) for _, co in pts))
        gr = 0
        for rv in ranks:
            gr = math.gcd(gr, int(rv)) if rv.denominator == 1 else gr
        chosen = reps[:4]
        nz = [e for e in reps if e[1][0] + e[1][2] != 0]
        if nz:
            mr = min(abs(e[1][0] + e[1][2]) for e in nz)
            chosen += [e for e in nz if abs(e[1][0] + e[1][2]) == mr][:1]
        irr = [e for e in reps if e[1][3] != 0]
        chosen += irr[:1]
        chosen = [e for i_, e in enumerate(chosen) if e not in chosen[:i_]]
        for hv, co in chosen:
            x, y, z, w = co
            cls = SS.lin(cl, [x, y, z, w])
            integral = all(Fr(v).denominator == 1 for v in cls.values())
            r, v = x + z, z - x
            chiF = q_ ** 3 * Pf * v
            nu = q_ ** 2 * Pf * (y * y - t_ * w * w)
            intc = s0 * cls.get(top, 0)
            Gm = [[y * (i == j) + w * St[i][j] for j in range(12)] for i in range(12)]
            aint = SS.integral([[q_ ** 2 * Pf * xx for xx in rr2] for rr2 in SS.gmm(Gm, Gm)])
            gok = r != 0 or (y == 0 and w == 0)
            row_ok = integral and muk(cls, cls) == 0 and intc == chiF and nu.denominator == 1 and \
                nu == r * chiF and aint and gok
            good &= row_ok
            rows.append({"height": str(hv), "x,y,z,w": [str(u) for u in co], "rank": str(r), "v": str(v),
                         "g": "%s + %s sqrt %d" % (y, w, t_), "chi(F)": str(chiF), "ok": row_ok})
        good &= len(reps) > 0
        ok &= good
        tab.append({"t": t_, "q": str(q_), "lattice": kind, "P": str(Pf), "height_bound": H, "classes": len(pts),
                    "up_to_symmetry": len(reps), "gcd_of_ranks": gr, "listed": rows, "ok": good})
        print("    (W) t=%d q=%s %s (P = %s): %d integral cone points with height <= %d (%d up to sign and tau), "
              "gcd of ranks %d (%.1fs): %s" % (t_, q_, kind, Pf, len(pts), H, len(reps), gr, time.time() - t0,
                                               good), flush=True)
        for rw in rows:
            print("        h=%s (x,y,z,w)=(%s) r=%s v=%s g=%s chi(F)=%s" % (rw["height"], ", ".join(rw["x,y,z,w"]),
                  rw["rank"], rw["v"], rw["g"], rw["chi(F)"]), flush=True)
    res["enumeration"] = tab
    return res, ok


# ----------------------------------------------------------------------------- (X)
def check_X():
    res, ok = {}, True
    # (i) symbolic identities, F = Q(sqrt t), elements u + v sqrt t
    t = sympy.symbols("t", positive=True)
    rt = sympy.sqrt(t)
    sy = sympy.symbols("m1a m1b m2a m2b g1a g1b g2a g2b sa sb lam qq")
    m1a, m1b, m2a, m2b, g1a, g1b, g2a, g2b, sa, sb, lam, qq = sy
    m1, m2, g1, g2 = m1a + m1b * rt, m2a + m2b * rt, g1a + g1b * rt, g2a + g2b * rt
    cj = lambda e: sympy.expand(e).subs(rt, -rt)
    Nm = lambda e: sympy.expand(e * cj(e))
    Tr = lambda e: sympy.expand(e + cj(e))
    nodes = [(1, m1), (-1, m1 - g1), (1, m2), (-1, m2 - g2)]
    p1, p2 = 2 * m1 - g1, 2 * m2 - g2
    G = g1 + g2
    i1 = sympy.expand(sum(e * h ** 2 for e, h in nodes) - (g1 * p1 + g2 * p2)) == 0
    sN = sum(e * Nm(h) for e, h in nodes)
    i2 = sympy.expand(sN - (Tr(p1 * cj(g1)) + Tr(p2 * cj(g2))) / 2) == 0
    i3 = sympy.expand(sum(e * Nm(h) ** 2 for e, h in nodes)
                      - (Tr(p1 * cj(g1)) * (Nm(p1) + Nm(g1)) + Tr(p2 * cj(g2)) * (Nm(p2) + Nm(g2))) / 4) == 0
    Hn = sum(e * Nm(h) * h for e, h in nodes)
    i4 = sympy.expand(Hn - (p1 * Tr(p1 * cj(g1)) + g1 * (Nm(p1) + Nm(g1))
                            + p2 * Tr(p2 * cj(g2)) + g2 * (Nm(p2) + Nm(g2))) / 4) == 0
    # the five B-conditions (10.5(a)) with sum eps = 0, at the two real embeddings
    x_, y_ = sympy.symbols("x y")
    pts = [(e, h, cj(h)) for e, h in nodes]
    conds = [sympy.expand(sum(e * (a_ ** 2 + qq) * b_ ** j for e, a_, b_ in pts)) for j in range(3)] + \
            [sympy.expand(sum(e * a_ ** i * (b_ ** 2 + qq) for e, a_, b_ in pts)) for i in range(2)]
    # with sum eps = 0 they read: sum eps h^2 = 0, Hn + q Gbar = 0, sum eps N(h)^2 + q conj(sum eps h^2) = 0,
    # conj(sum eps h^2) = 0, conj(Hn) + q G = 0 (a = h, b = hbar at the first real embedding; a^2 b = N(h) h)
    seh2 = sympy.expand(sum(e * h ** 2 for e, h in nodes))
    want5 = [seh2, sympy.expand(Hn + qq * cj(G)),
             sympy.expand(sum(e * Nm(h) ** 2 for e, h in nodes) + qq * cj(seh2)),
             cj(seh2), sympy.expand(cj(Hn) + qq * G)]
    i5 = all(sympy.expand(c_ - w_) == 0 for c_, w_ in zip(conds, want5))
    # coordinates (10.5(a)) and chi (Proposition 7.1, principal: W = 4): x = -sN/(2q), z = sN/(2q), y + w sqrt t = G
    xx = sum(e * (1 - Nm(h) / qq) for e, h in nodes) / 2
    zz = sum(e * (1 + Nm(h) / qq) for e, h in nodes) / 2
    i6 = sympy.simplify(xx + sN / (2 * qq)) == 0 and sympy.simplify(zz - sN / (2 * qq)) == 0
    chi = 8 * qq ** 2 * (xx ** 2 + zz ** 2) + 4 * qq * Tr(G * G)
    i7 = sympy.simplify(chi - 4 * (sN ** 2 + qq * Tr(G * G))) == 0
    # parametrization p1 = g2 s, p2 = -g1 s
    s = sa + sb * rt
    sub = {m1a: 0, m1b: 0, m2a: 0, m2b: 0}
    P1, P2 = g2 * s, -g1 * s
    Aq, Bq = Tr(g2 * cj(g1) * s), Tr(g1 * cj(g2) * s)
    E2 = sympy.expand(Tr(P1 * cj(g1)) * (Nm(P1) + Nm(g1)) + Tr(P2 * cj(g2)) * (Nm(P2) + Nm(g2)))
    i8 = sympy.expand(E2 - (Aq * (Nm(g2) * Nm(s) + Nm(g1)) - Bq * (Nm(g1) * Nm(s) + Nm(g2)))) == 0
    H4 = sympy.expand(P1 * Tr(P1 * cj(g1)) + g1 * (Nm(P1) + Nm(g1)) + P2 * Tr(P2 * cj(g2)) + g2 * (Nm(P2) + Nm(g2)))
    i9 = sympy.expand(H4 - (s * (g2 * Aq + g1 * Bq) + Nm(s) * (g1 * Nm(g2) + g2 * Nm(g1)) + g1 * Nm(g1)
                            + g2 * Nm(g2))) == 0
    sNs = sympy.expand((Tr(P1 * cj(g1)) + Tr(P2 * cj(g2))) / 2)
    i10 = sympy.expand(sNs - Tr(s * (g2 * cj(g1) - g1 * cj(g2))) / 2) == 0
    # s rational: (II) = s tau (N(g2) - N(g1)) (s^2 - 1), tau = Tr(g2 gbar1)
    tau_ = Tr(g2 * cj(g1))
    i11 = sympy.expand(E2.subs(sb, 0) - sa * tau_ * (Nm(g2) - Nm(g1)) * (sa ** 2 - 1)) == 0
    # g2 = lam g1: (II) and 4 Hn
    E2p = sympy.expand(E2.subs({g2a: lam * g1a, g2b: lam * g1b}))
    i12 = sympy.expand(E2p - lam * Nm(g1) * Tr(s) * (Nm(s) - 1) * Nm(g1) * (lam ** 2 - 1)) == 0
    H4p = sympy.expand(H4.subs({g2a: lam * g1a, g2b: lam * g1b}))
    i13 = sympy.expand(H4p - g1 * Nm(g1) * (1 + lam) * (lam * (s * s + 2 * Nm(s)) + 1 - lam + lam ** 2)) == 0
    good1 = all([i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13])
    ok &= good1
    res["identities"] = [i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13]
    print("    (X) symbolic identities (B-conditions, sum eps N, Hn, coordinates, chi = 4 [(sum eps N)^2 + "
          "q Tr(G^2)], parametrization, the cases s rational and g_2 = lambda g_1): %s: %s"
          % (res["identities"], good1), flush=True)
    # (ii) content of Theta_1 Theta_2 on the principal free lattice, n = 4
    cont = {}
    for tt in (2, 3, 5, 6, 7, 13, 17, 21, 57):
        P, s0 = SS.setup(tt, SS.delta(tt), SS.delta(tt))
        T12 = SS.add((Fr(1, 4), P[(2, 0)]), (-Fr(1, 4 * tt), P[(0, 2)]))
        dF = tt if tt % 4 == 1 else 4 * tt
        cont[tt] = (str(SS.content(T12)), str(Fr(1, dF) if tt % 4 == 1 else Fr(2, dF)))
    good2 = all(a_ == b_ for a_, b_ in cont.values())
    ok &= good2
    res["content_T12"] = cont
    print("    (X) content of Theta_1 Theta_2 on O_F^4 (principal): %s (= 1/d_F or 2/d_F): %s" % (cont, good2),
          flush=True)
    # (iii) the case g_1 = g_2 = g irrational, chi = 4 Tr(a) <= 44, a = 4 q g^2, over all real quadratic fields
    found, nas, cands = [], 0, 0
    for tt in range(2, 400):
        if any(tt % (pp * pp) == 0 for pp in range(2, int(math.isqrt(tt)) + 1)):
            continue
        sols, na = dioph(tt, 11)
        nas += na
        cands += len(sols[0])
        found += sols[1]
    # positive control: the configuration of 10.5(d) (g_1 = g_2 = 8 - 3 sqrt 5, q = 25631) is recovered
    g5, q5 = (Fr(8), Fr(-3)), Fr(25631)
    g5sq = (g5[0] ** 2 + 5 * g5[1] ** 2, 2 * g5[0] * g5[1])
    a5 = (4 * q5 * g5sq[0], 4 * q5 * g5sq[1])
    ks5, sols5 = dioph_a(5, a5)
    p1 = (Fr(-2008), Fr(981))                                   # 2 m_1 - g
    s5 = ((p1[0] * g5[0] - 5 * p1[1] * g5[1]) / (g5[0] ** 2 - 5 * g5[1] ** 2),
          (p1[1] * g5[0] - p1[0] * g5[1]) / (g5[0] ** 2 - 5 * g5[1] ** 2))
    control = any(s_ in (s5, (-s5[0], -s5[1])) for _, s_ in sols5)
    good3 = not found and nas > 0 and control
    ok &= good3
    res["g1=g2 case"] = {"a's": nas, "candidate scalings": cands, "solutions": found,
                         "control_10.5(d)": {"k": ks5, "s": str(s5), "recovered": control}}
    print("    (X) g_1 = g_2 irrational, Tr(a) <= 11, all squarefree t < 400: %d elements a, %d candidate "
          "scalings, solutions %s; control: the 10.5(d) configuration is recovered from its a (k = %s, s = %s): "
          "%s: %s" % (nas, cands, found, ks5, s5, control, good3), flush=True)
    # (iv) brute force
    bf = {}
    for tt, GB, MB in ((5, 5, 15), (13, 4, 10), (2, 4, 10)):
        bf[tt] = brute(tt, GB, MB)
    good4 = all(v["off_P"] == 0 and v["in_B"] > 0 for v in bf.values())
    ok &= good4
    res["brute"] = bf
    print("    (X) brute force (g totally positive with Tr(g) <= 2 GB, m in a box): %s: %s" % (bf, good4),
          flush=True)
    # (v) the classes of 10.5(c), (d) (t = 5): chi
    r5 = sympy.sqrt(5)
    ex = {}
    for nm, g_, mm1, qv_ in (("10.5(d) O_D(M) + O_D'(M')", 8 - 3 * r5, -1000 + 489 * r5, 25631),
                             ("10.5(c) a = phi", sympy.Integer(1), (1 - r5) / 2, 1)):
        mm2 = g_ - mm1 if nm.startswith("10.5(c)") else 1008 - 492 * r5
        nd = [(1, mm1), (-1, mm1 - g_), (1, mm2), (-1, mm2 - g_)]
        N5 = lambda e: sympy.expand(e * sympy.expand(e).subs(r5, -r5))
        sN5 = sum(e * N5(h) for e, h in nd)
        Hn5 = sympy.expand(sum(e * N5(h) * h for e, h in nd))
        GG = 2 * g_
        qv = sympy.simplify(-Hn5 * GG / N5(GG))
        chi5 = 4 * (sN5 ** 2 + qv * sympy.expand(GG ** 2 + sympy.expand(GG ** 2).subs(r5, -r5)))
        ex[nm] = (str(qv), str(sympy.simplify(chi5)))
    chi_sig = 4 * 1 * ((18 + 10 * r5) ** 2 + (18 - 10 * r5) ** 2)
    ex["10.5(d) 10 O_D(L_phi^2) - 3 O_D'(L_2) (signed)"] = ("1", str(sympy.expand(chi_sig)))
    good5 = ex["10.5(d) O_D(M) + O_D'(M')"] == ("25631", "89400928") and ex["10.5(c) a = phi"] == ("1", "32") \
        and ex["10.5(d) 10 O_D(L_phi^2) - 3 O_D'(L_2) (signed)"][1] == "6592"
    ok &= good5
    res["examples"] = ex
    print("    (X) (q, chi) of the classes of 10.5(c), (d): %s (> 44, resp. > 28 on P): %s" % (ex, good5),
          flush=True)
    return res, ok


def ratsqrt(x):
    x = Fr(x)
    if x < 0:
        return None
    a, b = math.isqrt(x.numerator), math.isqrt(x.denominator)
    return Fr(a, b) if a * a == x.numerator and b * b == x.denominator else None


def dioph_a(t, a):
    """the case g_1 = g_2 = g for one totally positive irrational a = 4 q g^2 in O_F (pair (u, v): u + v sqrt t)
    with N(a) a square n^2.  Then g = mu g_a (g_a = a + n, mu > 0 rational with mu g_a in O_F), and s = u + v sqrt t
    solves s^2 + 2 N(s) + 1 = -K rho, rho = gbar_a^2/N(g_a) = r0 + r1 sqrt t, K = 1/(mu^4 n tau^2), tau =
    Tr(a) + 2 n.  Rational u needs D = 1 + 2 K r0 + (4 r0^2 - 3) K^2 to be a square, i.e. (M n tau^2 + r0)^2 +
    3 (r0^2 - 1) = square (M = mu^4): finitely many M by factorization.  Returns the candidate k (mu = k/c,
    g_a = c * primitive) and the solutions (mu, s) with u, v rational."""
    om = (Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else (Fr(0), Fr(1))

    def mul(x, y):
        return (x[0] * y[0] + t * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    def N(x):
        return x[0] * x[0] - t * x[1] * x[1]
    n = ratsqrt(N(a))
    ga = (a[0] + n, a[1])
    tau = 2 * a[0] + 2 * n
    Nga = N(ga)
    assert Nga == n * tau and ga[0] > 0 and Nga > 0
    rb = mul((ga[0], -ga[1]), (ga[0], -ga[1]))
    r0, r1 = rb[0] / Nga, rb[1] / Nga
    bb = ga[1] / om[1]
    aa = ga[0] - bb * om[0]
    L0 = aa.denominator * bb.denominator
    c = Fr(math.gcd(int(aa * L0), int(bb * L0)), L0)          # g_a = c * (primitive element of O_F)
    B = 3 * (r0 * r0 - 1)
    base = n * tau * tau / c ** 4                               # M n tau^2 = k^4 base, mu = k/c
    L = 1
    for xq in (base, r0, B):
        L = L * xq.denominator // math.gcd(L, xq.denominator)
    L2B = L * L * B
    assert L2B.denominator == 1 and L2B > 0
    ks, sols = [], []
    for d in sympy.divisors(int(L2B)):
        e = int(L2B) // d
        if e <= d or (e - d) % 2:
            continue
        Z = Fr(e - d, 2) / L
        k4 = (Z - r0) / base
        if k4 <= 0 or k4.denominator != 1:
            continue
        kk = round(float(k4) ** 0.25)
        for k_ in (kk - 1, kk, kk + 1):
            if k_ > 0 and k_ ** 4 == k4:
                ks.append(k_)
                mu = Fr(k_) / c
                K = 1 / (mu ** 4 * n * tau * tau)
                sD = ratsqrt(1 + 2 * K * r0 + (4 * r0 * r0 - 3) * K * K)
                uu = ratsqrt((sD - 1 - K * r0) / 6) if sD is not None else None
                if uu:
                    for su in (uu, -uu):
                        sols.append((mu, (su, -K * r1 / (2 * su))))
    return ks, sols


def dioph(t, TR):
    """all totally positive irrational a in O_F with Tr(a) <= TR and N(a) a square, through dioph_a."""
    om = (Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else (Fr(0), Fr(1))

    def inOF(x):
        b = x[1] / om[1]
        a = x[0] - b * om[0]
        return a.denominator == 1 and b.denominator == 1
    cands, found, na = [], [], 0
    for A2 in range(1, TR + 1):
        u = Fr(A2, 2)
        vmax = int(u / math.sqrt(t)) + 2
        for k in range(-4 * vmax, 4 * vmax + 1):
            v = Fr(k, 2)
            a = (u, v)
            if v == 0 or not inOF(a) or u * u - t * v * v <= 0 or ratsqrt(u * u - t * v * v) is None:
                continue
            na += 1
            ks, sols = dioph_a(t, a)
            cands += [(t, str(a), k_) for k_ in ks]
            found += [(t, str(a), str(mu), str(s_)) for mu, s_ in sols]
    return (cands, found), na


def brute(t, GB, MB):
    om = (Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else (Fr(0), Fr(1))

    def el(a, b):
        return (a + b * om[0], b * om[1])

    def mul(x, y):
        return (x[0] * y[0] + t * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    def N(x):
        return x[0] * x[0] - t * x[1] * x[1]

    def inOF(x):
        b = x[1] / om[1]
        a = x[0] - b * om[0]
        return a.denominator == 1 and b.denominator == 1
    gs = [el(a, b) for a in range(-3 * GB, 3 * GB + 1) for b in range(0, 2 * GB + 1)]
    gs = [g for g in gs if g[0] > 0 and N(g) > 0 and 2 * g[0] <= 2 * GB]
    ms = [el(a, b) for a in range(-MB, MB + 1) for b in range(-MB, MB + 1)]
    inB, offP, least = 0, 0, None
    for g1, g2 in itertools.combinations_with_replacement(gs, 2):
        G = (g1[0] + g2[0], g1[1] + g2[1])
        Ng2 = N(g2)
        for m1 in ms:
            p1 = (2 * m1[0] - g1[0], 2 * m1[1] - g1[1])
            num = mul(g1, p1)
            num = mul(num, (g2[0], -g2[1]))
            p2 = (-num[0] / Ng2, -num[1] / Ng2)
            m2 = ((p2[0] + g2[0]) / 2, (p2[1] + g2[1]) / 2)
            if not inOF(m2):
                continue
            nodes = [(1, m1), (-1, (m1[0] - g1[0], m1[1] - g1[1])), (1, m2), (-1, (m2[0] - g2[0], m2[1] - g2[1]))]
            if sum(e * N(h) ** 2 for e, h in nodes) != 0:
                continue
            Hn = (sum(e * N(h) * h[0] for e, h in nodes), sum(e * N(h) * h[1] for e, h in nodes))
            HG = mul(Hn, G)
            if HG[1] != 0 or -HG[0] / N(G) <= 0:
                continue
            q = -HG[0] / N(G)
            sN = sum(e * N(h) for e, h in nodes)
            inB += 1
            if not (sN == 0 and G[1] == 0):
                offP += 1
                chi = 4 * (sN ** 2 + q * 2 * (G[0] ** 2 + t * G[1] ** 2))
                least = chi if least is None else min(least, chi)
    return {"pairs_g": len(gs), "box_m": MB, "in_B": inB, "off_P": offP, "least_chi_off_P": str(least)}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "cone6.json"
    res, allok, secs = {}, True, {}
    for key, fn in (("U", check_U), ("V", check_V), ("W", check_W), ("X", check_X)):
        t0 = time.time()
        r, ok = fn()
        secs[key] = round(time.time() - t0, 1)
        res[key] = r
        allok &= ok
        print("(%s) %s  [%.1fs]" % (key, ok, secs[key]), flush=True)
    with open(out, "w", newline="\n") as f:
        json.dump({"results": res, "seconds": secs, "_all_checks_pass": allok}, f, indent=1, sort_keys=True,
                  default=str)
    print("ALL CHECKS PASS:", allok)


if __name__ == "__main__":
    main()
