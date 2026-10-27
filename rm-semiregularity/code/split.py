"""
split.py -- checks for ../NOTE.txt, Proposition 9.19, Lemma 9.20, Theorem 9.21, Corollary 9.22 and Proposition 9.23
(OPEN 11(l): does Markman's split model X_0 x X_0, E' = I' [x] I occur inside an RM family?).

Frames.  As in semireg.py: H^1(X, C) has the basis dz_j (index j) and dzbar_j (index n + j), Theta-type 2-forms are
sums of dzbar_l ^ dz_l, and HT^2 = wedge^2 W, W = <dzbar_j, d_j>, acts by wedge and contraction.  so(V_C) =
wedge^2 V_C, V_C = H^1 + H_1 (4n labels: wedge with dz_j, dzbar_j, contraction with their duals), acts on H^*(X, C)
by e_i ^ e_j -> (c_i c_j - c_j c_i)/2.  The "eigen-frame" is the frame of NOTE Theorem 3.1 (the two eigenspaces of
F are the coordinate blocks; rational coefficients).  The "product frame" is that of X_0 x X_0 (coordinates 0, 1, 2
and 3, 4, 5 on the two factors, Theta' = sum_{j<3} dzbar_j dz_j, Theta'' the same for j >= 3, M = (1/2) sum_{j<3}
(dzbar_j dz_{j+3} + dzbar_{j+3} dz_j)); there O = Z + Z omega acts through the regular representation, the
eigen-projections are phi_k = pi_1 + omega_k pi_2 (omega_k = sigma_k(omega)), phi_k^* theta_0 = Theta' +
2 omega_k M + omega_k^2 Theta'', and Theta_k = sigma_k(gamma) phi_k^* theta_0 for the polarization
theta_0 (x) Tr(gamma .) (gamma totally positive).  Classes there have coefficients in F = Q(sqrt t), stored as
pairs (rational part, sqrt t part); ranks over F are computed exactly as half the rank over Q of the realified
matrix [[A, t B], [B, A]] (python-flint, fmpz).

 (Y) NOTE Lemma 9.20, Theorem 9.21(a), (b) (eigen-frame, exact over Q).  (i) One block (a complex 3-torus, V of
     dimension 12): for u = a C + b S (q = 3, 1/2; five (a, b)) the stabilizer of u in so(V) has dimension 35
     (= dim sl_6), the stabilizer of the line C u is the same, and v -> v.u is injective on V.  (ii) n = 6,
     rational cone points C_1 C_2, S_1 S_2, (C_1 + 15 S_1)(C_2 + 15 S_2), (C_1 + 2 S_1)(3 C_2 - S_2) (q = 3) and
     Markman's class f(Theta_1) f(Theta_2), f = 1 - (N/2) T^2 + (N/3) T^3 (N = 2, 4): stab = line stab = 70,
     contained in so(V_1) + so(V_2) with 35 + 35, annihilator in V zero.  (iii) The rest of B (NOTE 9.21(d)): the
     pointwise stabilizer Lie Spin(V)_B of B (the common kernel of the four basis classes) has dimension
     70 = 35 + 35; five rational points of B off the cone (Mukai square != 0) and off P, P~ have stab = line stab =
     70 = 35 + 35, so stab(c) = Lie Spin(V)_B there; alpha, beta, alpha~, beta~ have 143 (= dim sl_12) >= 70.
     Comparison with n = 4: cone points 42 = 21 + 21 (21 = dim spin(7)), generic 30 = 15 + 15 (the reason why
     dim K(c) = 10 > 8 on the n = 4 cone, NOTE 3.1(iv)).
 (Z) NOTE Propositions 9.19 and 9.23 at the diagonal member X_0 (x) O = X_0 x X_0 (product frame, exact over F),
     for (t, omega, gamma, q, N) = (5, phi, 1, 2, 3), (5, phi, 1, 1, 3), (5, phi, phi^2, 3, 2),
     (2, sqrt 2, 1, 1, 2), (2, 1 + sqrt 2, 1, 1, 2), (57, (1 + sqrt 57)/2, 1, 3, 4) (N = d + 1).  (i) The regular
     representation is self-adjoint for Tr(gamma x y), (1, omega_k) are its eigenvectors, Theta_1 + Theta_2 =
     Tr(gamma) Theta' + 2 Tr(gamma omega) M + Tr(gamma omega^2) Theta'', and the four basis classes of B are
     rational.  (ii) dim K(B) = 18, K(B) n H^1(T) = R (12), K(B) n (H^2(O) + H^0(wedge^2 T)) = G (6).
     (iii) dim K(ch E') = 18 = pi_1^* K(ch I) + pi_2^* K(ch I) (9 + 9, exact kernels on X_0), K(ch E') n H^1(T) =
     12, K(ch E') n (H^2(O) + H^0(wedge^2 T)) = 0.  (iv) R n K(ch E') = K(B) n K(ch E') = R_0 (6), the diagonal
     deformations sum_{a,b} c_ab (dzbar_a d_b + dzbar_{a+3} d_{b+3}), c symmetric (lower bound by the explicit
     basis, exact).  (v) Markman's 9-dimensional diagonal D = {pi_1^* a + pi_2^* a^* : a in K(alpha_0 + beta_0)}
     ([M2] Corollary 8.5.2) lies in K((alpha_0 + beta_0)' (alpha_0 + beta_0)''), e^{-Theta' - Theta''} times that
     class is ch(E'), D_E' n K(B) has dimension 6 (computed as D n K(e^{Theta' + Theta''} B)), and R_0 lies in D.
     (vi) ch(E') is not in the span of the 16 classes Theta_1^i Theta_2^j, so it is in B for no q.
 (AA) NOTE Corollary 9.22.  (i) At t = 5, omega = phi, N = 3 and gamma = 1, 2, 1/3 with q = (N - 1)/gamma^2
     (product frame): the Galois-twisted class c_tw = phi_1^* ch(I) phi_2^* ch(I) has rational coefficients and does
     not depend on gamma, e^{Theta/gamma} c_tw = (C_1 + S_1/gamma)(C_2 + S_2/gamma), a class of B with coordinates
     r = 1, q v = gamma^-2, g = gamma^-1 (q r v = N(g): a cone class), while e^{Theta/gamma} c_tw = (C_1 + S_1)
     (C_2 + S_2) holds only for gamma = 1; c_tw != ch(E').  The stabilizers of c_tw, of C_1 C_2, of ch(E') and of
     the cone classes (C_1 + S_1/gamma)(C_2 + S_2/gamma) (gamma = 2, 1/3) have dimension 70; those of c_tw, C_1 C_2
     and the cone classes meet so(V'_1), so(V'_2) (the product blocks) in 0 and so(V_1), so(V_2) (the eigen-blocks,
     coefficients in F) in 35 each, and the reverse holds for ch(E').  (ii) On the lattice H_1(X_0, Z) (x) O_F with
     E = J (theta_0 principal, gamma = 1; code of cone6.py, exact integral cohomology) for t = 2, 3, 5, 13, 57 and
     N = 2, ..., 5: e^{-Theta} (C_1 + S_1)(C_2 + S_2) (q = N - 1) is integral, has rank 1, c_1 = 0, Mukai square 0,
     and int = d_F^3 (2N)^2 = d_F^3 int ch(E'); ch(E') = f(Theta') f(Theta'') is integral there as well.
     (iii) The proof of integrality in NOTE 9.22(b), on the integral basis of H^*(X, Z), X = X_0 (x) O, given by
     the monomials in pr_1^* x_j, pr_1^* y_j, pr_2^* x_j, pr_2^* y_j (x_j, y_j a symplectic basis of H^1(X_0, Z),
     theta_0 = sum x_j y_j), with phi_k^* x = pr_1^* x + omega_k pr_2^* x: for eleven Z-bases (1, omega) of orders
     O (the maximal orders of Q(sqrt t), t = 2, 3, 5, 10, 13, 21, 57, the orders Z[sqrt 5], Z[3 phi], Z[1 + sqrt 13]
     of index 2, 3, 2, and the basis (1, 1 + sqrt 2) of Z[sqrt 2]) and N = 1, ..., 9: (phi_k^* theta_0)^4 = 0,
     phi_1^* ch(I) has coefficients in O, c_tw is rational and integral, of rank 1 with c_1 = 0 and Mukai square 0,
     and int c_tw = disc(O)^3 (2N)^2; ch(E') is integral with int 4 N^2 and Mukai square 0.  For four of the bases,
     N = 3 and gamma = 1, 2, 1/3: c_tw = e^{-Theta/gamma} (C_1 + S_1/gamma)(C_2 + S_2/gamma) (Theta_k =
     gamma phi_k^* theta_0, q = (N - 1)/gamma^2), computed as a product in the exterior algebra, and the formula
     with C_k + S_k fails for gamma != 1.
Usage:  python3 -I split.py OUT.json      Expected last line "ALL CHECKS PASS: True".  About 35 seconds.
"""

import itertools
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint  # noqa: E402
import sympy  # noqa: E402

import semireg as SR  # noqa: E402
import shape_search as SS  # noqa: E402
import cone6 as C6  # noqa: E402

ZERO, ONE = (Fr(0), Fr(1))
ONEF = (Fr(1), Fr(0))


# ----------------------------------------------------------------------------- exact ranks
def rank_Q(cols):
    """rank over Q of the matrix whose columns are the dicts cols (row key -> rational)."""
    keys = sorted(set(k for c in cols for k in c), key=repr)
    if not keys:
        return 0
    den = 1
    for c in cols:
        for v in c.values():
            dv = Fr(v).denominator
            den = den * dv // math.gcd(den, dv)
    idx = {k: i for i, k in enumerate(keys)}
    M = flint.fmpz_mat(len(keys), len(cols))
    for j, c in enumerate(cols):
        for k, v in c.items():
            M[idx[k], j] = int(Fr(v) * den)
    return M.rank()


def rank_F(cols, t):
    """rank over F = Q(sqrt t) of columns with entries (a, b) = a + b sqrt t (realification)."""
    out = []
    for c in cols:
        re, rt = {}, {}
        for k, (a, b) in c.items():
            if a:
                re[(k, 0)] = a
                rt[(k, 1)] = a
            if b:
                re[(k, 1)] = b
                rt[(k, 0)] = t * b
        out += [re, rt]
    r = rank_Q(out)
    assert r % 2 == 0
    return r // 2


def kernel_Q(n, cls, basis):
    """exact rational kernel basis of contraction with the rational class cls on the given HT^2 basis elements."""
    cols = [SR.act(x, cls, n) for x in basis]
    keys = sorted(set(k for c in cols for k in c))
    den = 1
    for c in cols:
        for v in c.values():
            den = den * Fr(v).denominator // math.gcd(den, Fr(v).denominator)
    idx = {k: i for i, k in enumerate(keys)}
    M = flint.fmpz_mat(len(keys), len(cols))
    for j, c in enumerate(cols):
        for k, v in c.items():
            M[idx[k], j] = int(Fr(v) * den)
    X, nul = M.nullspace()
    return [{basis[i]: Fr(int(X[i, j])) for i in range(len(basis)) if X[i, j] != 0} for j in range(nul)]


# ----------------------------------------------------------------------------- so(V) acting on H^*
def op(lab, vec, n):
    return SR.wedge(lab, vec) if lab < 2 * n else SR.contr(lab - 2 * n, vec)


def biv(i, j, vec, n):
    out = {}
    for k, v in op(i, op(j, vec, n), n).items():
        out[k] = out.get(k, 0) + v / 2
    for k, v in op(j, op(i, vec, n), n).items():
        out[k] = out.get(k, 0) - v / 2
    return {k: v for k, v in out.items() if v != 0}


def block_labels(n, coords):
    out = set()
    for j in coords:
        out |= {j, n + j, 2 * n + j, 3 * n + j}
    return out


def stab_rational(c, n, blocks):
    """dimensions of the stabilizer of the rational class c in so(V_C), of the line, of its intersections with the
    so(V_k) of coordinate blocks, and of the annihilator of c in V_C."""
    B = list(itertools.combinations(range(4 * n), 2))
    cols = [biv(i, j, c, n) for i, j in B]
    r = rank_Q(cols)
    res = {"dim_so": len(B), "stab": len(B) - r, "line_stab": len(B) + 1 - rank_Q(cols + [c])}
    for nm, lab in blocks.items():
        sel = [cols[k] for k, (i, j) in enumerate(B) if i in lab and j in lab]
        res["stab n so(%s)" % nm] = len(sel) - rank_Q(sel)
    res["ann_V"] = 4 * n - rank_Q([op(l, c, n) for l in range(4 * n)])
    return res


def stab_multi(classes, n, blocks):
    """dimension of the common stabilizer in so(V_C) of the rational classes (e.g. the pointwise stabilizer of B),
    and of its intersections with the so(V_k) of coordinate blocks."""
    B = list(itertools.combinations(range(4 * n), 2))
    cols = []
    for i, j in B:
        col = {}
        for ci, c in enumerate(classes):
            for k, v in biv(i, j, c, n).items():
                col[(ci, k)] = v
        cols.append(col)
    res = {"stab": len(B) - rank_Q(cols)}
    for nm, lab in blocks.items():
        sel = [cols[k] for k, (i, j) in enumerate(B) if i in lab and j in lab]
        res["stab n so(%s)" % nm] = len(sel) - rank_Q(sel)
    return res


# ----------------------------------------------------------------------------- classes with coefficients in F
class Fld:
    def __init__(self, t):
        self.t = t

    def mul(self, x, y):
        return (x[0] * y[0] + self.t * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    @staticmethod
    def conj(x):
        return (x[0], -x[1])

    def inv(self, x):
        nx = x[0] * x[0] - self.t * x[1] * x[1]
        return (x[0] / nx, -x[1] / nx)

    @staticmethod
    def tr(x):
        return 2 * x[0]


class FClasses:
    """classes A + sqrt t B as pairs of bitmask dicts; product frame of X_0 x X_0 (n = 6)."""

    def __init__(self, t, n=6):
        self.t, self.n, self.K = t, n, Fld(t)

    def add(self, *terms):
        A, B = {}, {}
        for (a, b), (X, Y) in terms:
            for k, v in X.items():
                if a:
                    A[k] = A.get(k, 0) + a * v
                if b:
                    B[k] = B.get(k, 0) + b * v
            for k, v in Y.items():
                if a:
                    B[k] = B.get(k, 0) + a * v
                if b:
                    A[k] = A.get(k, 0) + self.t * b * v
        return ({k: v for k, v in A.items() if v}, {k: v for k, v in B.items() if v})

    def mul(self, P, Q):
        A = SS.add((1, SS.wmul(P[0], Q[0])), (self.t, SS.wmul(P[1], Q[1])))
        B = SS.add((1, SS.wmul(P[0], Q[1])), (1, SS.wmul(P[1], Q[0])))
        return (A, B)

    def poly(self, coefs, X):
        out, P = ({}, {}), ({0: Fr(1)}, {})
        for ck in coefs:
            if not isinstance(ck, tuple):
                ck = (Fr(ck), Fr(0))
            out = self.add((ONEF, out), (ck, P))
            P = self.mul(P, X)
        return out

    def act(self, x, P):
        out = {k: (v, Fr(0)) for k, v in SR.act(x, P[0], self.n).items()}
        for k, v in SR.act(x, P[1], self.n).items():
            out[k] = (out.get(k, (Fr(0), Fr(0)))[0], v)
        return out

    def ann_cols(self, combs, classes):
        """columns (stacked over the classes) of contraction with the HT^2 elements combs (label -> F pair)."""
        cols = []
        for comb in combs:
            col = {}
            for ci, P in enumerate(classes):
                for x, cx in comb.items():
                    for k, val in self.act(x, P).items():
                        u = self.K.mul(cx, val)
                        o = col.get((ci, k), (Fr(0), Fr(0)))
                        col[(ci, k)] = (o[0] + u[0], o[1] + u[1])
            cols.append({k: v for k, v in col.items() if v != (0, 0)})
        return cols

    def kdim(self, combs, classes):
        return len(combs) - rank_F(self.ann_cols(combs, classes), self.t)


def rat(D):
    return (D, {})


def tf(a, b, n=6):
    """dzbar_b ^ dz_a."""
    return SR.wedge(n + b, SR.wedge(a, {0: Fr(1)}))


def product_frame():
    Tp = SS.add(*[(1, tf(j, j)) for j in range(3)])
    Tpp = SS.add(*[(1, tf(j + 3, j + 3)) for j in range(3)])
    M = SS.add(*[(Fr(1, 2), tf(j, j + 3)) for j in range(3)], *[(Fr(1, 2), tf(j + 3, j)) for j in range(3)])
    return Tp, Tpp, M


def unit_combs(kinds=None, n=6):
    return [{x: ONEF} for x in SR.ht2(n) if kinds is None or x[0] in kinds]


def omega_of(t):
    return (Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else (Fr(0), Fr(1))


def diagonal_member(t, om, c, q, N):
    """the data of NOTE 9.19, 9.22, 9.23 in the product frame."""
    FC, K = FClasses(t), Fld(t)
    Tp, Tpp, M = product_frame()

    def pull(lam):            # phi_lam^* theta_0 = Theta' + 2 lam M + lam^2 Theta''
        return FC.add((ONEF, rat(Tp)), ((2 * lam[0], 2 * lam[1]), rat(M)), (K.mul(lam, lam), rat(Tpp)))
    th1, th2 = pull(om), pull(K.conj(om))
    T1, T2 = FC.add((c, th1)), FC.add((K.conj(c), th2))
    q = Fr(q)
    C1, S1 = FC.poly([1, 0, -q / 2], T1), FC.poly([0, 1, 0, -q / 6], T1)
    C2, S2 = FC.poly([1, 0, -q / 2], T2), FC.poly([0, 1, 0, -q / 6], T2)
    Bc = [FC.mul(C1, C2), FC.mul(S1, S2), FC.add((ONEF, FC.mul(S1, C2)), (ONEF, FC.mul(C1, S2))),
          FC.add(((Fr(0), Fr(1)), FC.mul(S1, C2)), ((Fr(0), Fr(-1)), FC.mul(C1, S2)))]
    f = [1, 0, -Fr(N, 2), Fr(N, 3)]
    chE = FC.mul(FC.poly(f, rat(Tp)), FC.poly(f, rat(Tpp)))
    return dict(FC=FC, K=K, Tp=Tp, Tpp=Tpp, M=M, th1=th1, th2=th2, T1=T1, T2=T2, C1=C1, S1=S1, C2=C2, S2=S2,
                B=Bc, chE=chE, f=f)


# ----------------------------------------------------------------------------- (Y)
def check_Y():
    res, ok = {}, True
    T1, T2 = sympy.symbols("T1 T2")
    # (i) one block
    one = {}
    for q in (Fr(3), Fr(1, 2)):
        for a, b in ((1, 0), (0, 1), (1, 2), (3, -1), (2, 5)):
            u = C6.block_class(3, {0: Fr(a), 1: Fr(b), 2: -q * a / 2, 3: -q * b / 6})
            one[(str(q), a, b)] = stab_rational(u, 3, {})
    good1 = all(v["stab"] == 35 and v["line_stab"] == 35 and v["ann_V"] == 0 for v in one.values())
    ok &= good1
    res["one_block"] = {str(k): v for k, v in one.items()}
    print("    (Y) one block (V of dimension 12), u = a C + b S: stab = line stab = 35, ann_V = 0 at %d classes: %s"
          % (len(one), good1), flush=True)
    # (ii) n = 6 cone points and Markman's class
    n = 6
    blocks = {"V_1": block_labels(n, range(3)), "V_2": block_labels(n, range(3, 6))}
    q = 3
    (C1, S1), (C2, S2) = C6.cs(3, q, T1), C6.cs(3, q, T2)
    fN = lambda T, N: 1 - sympy.Rational(N, 2) * T ** 2 + sympy.Rational(N, 3) * T ** 3
    pts = [("C1 C2", C1 * C2), ("S1 S2", S1 * S2), ("(C1+15S1)(C2+15S2)", (C1 + 15 * S1) * (C2 + 15 * S2)),
           ("(C1+2S1)(3C2-S2)", (C1 + 2 * S1) * (3 * C2 - S2)),
           ("f(T1)f(T2) N=2", fN(T1, 2) * fN(T2, 2)), ("f(T1)f(T2) N=4", fN(T1, 4) * fN(T2, 4))]
    six = {}
    for nm, e in pts:
        six[nm] = stab_rational(C6.frame_class(6, e, T1, T2), n, blocks)
    want = {"dim_so": 276, "stab": 70, "line_stab": 70, "stab n so(V_1)": 35, "stab n so(V_2)": 35, "ann_V": 0}
    good2 = all(v == want for v in six.values())
    ok &= good2
    res["n6_cone"] = six
    print("    (Y) n = 6 cone points and f(Theta_1) f(Theta_2): %s: %s" % (six, good2), flush=True)
    # (iii) comparison
    cmpr = {}
    al, be, alt, bet = SR.b_basis(6, Fr(3))
    cmpr["n=6 B pointwise (Lie Spin(V)_B)"] = stab_multi([al, be, alt, bet], 6, blocks)
    top6 = (1 << 12) - 1
    generic = {"2 alpha - beta + 3 alpha~ + beta~": (2, -1, 3, 1),
               "alpha + 2 beta - alpha~ + 3 beta~": (1, 2, -1, 3),
               "5 alpha - 3 beta + 2 alpha~ + 7 beta~": (5, -3, 2, 7),
               "-2 alpha + beta + 4 alpha~ + beta~": (-2, 1, 4, 1),
               "alpha + beta~": (1, 0, 0, 1)}
    for nm, co in generic.items():
        cl = SS.add(*[(x, b) for x, b in zip(co, (al, be, alt, bet))])
        r = stab_rational(cl, 6, blocks)
        r["Mukai square != 0 (off the cone)"] = SS.wmul(SS.tau(cl), cl).get(top6, 0) != 0
        cmpr["n=6 generic " + nm] = r
    for nm, cl in (("alpha", al), ("beta", be), ("alpha~", alt), ("beta~", bet)):
        cmpr["n=6 " + nm] = stab_rational(cl, 6, blocks)
    blocks4 = {"V_1": block_labels(4, range(2)), "V_2": block_labels(4, range(2, 4))}
    (c1, s1), (c2, s2) = C6.cs(2, 1, T1), C6.cs(2, 1, T2)
    for nm, e in (("n=4 C1 C2", c1 * c2), ("n=4 (C1+2S1)(C2+2S2)", (c1 + 2 * s1) * (c2 + 2 * s2))):
        cmpr[nm] = stab_rational(C6.frame_class(4, e, T1, T2), 4, blocks4)
    al4, be4, alt4, bet4 = SR.b_basis(4, Fr(1))
    cmpr["n=4 generic"] = stab_rational(SS.add((2, al4), (-1, be4), (3, alt4), (5, bet4)), 4, blocks4)
    good3 = (cmpr["n=6 B pointwise (Lie Spin(V)_B)"] == {"stab": 70, "stab n so(V_1)": 35, "stab n so(V_2)": 35}
             and all(cmpr["n=6 generic " + nm] == dict(want, **{"Mukai square != 0 (off the cone)": True})
                     for nm in generic)
             and all(cmpr["n=6 " + nm]["stab"] == 143 for nm in ("alpha", "beta", "alpha~", "beta~"))
             and all(cmpr[k]["stab"] == 42 and cmpr[k]["stab n so(V_1)"] == 21 and cmpr[k]["stab n so(V_2)"] == 21
                     for k in ("n=4 C1 C2", "n=4 (C1+2S1)(C2+2S2)"))
             and cmpr["n=4 generic"]["stab"] == 30 and cmpr["n=4 generic"]["stab n so(V_1)"] == 15)
    ok &= good3
    res["comparison"] = cmpr
    print("    (Y) all of B (pointwise stabilizer, five points off the cone, alpha, beta, alpha~, beta~) and n = 4 "
          "(cone and generic): %s: %s" % (cmpr, good3), flush=True)
    return res, ok


# ----------------------------------------------------------------------------- (Z)
CONFIGS = [(5, (Fr(1, 2), Fr(1, 2)), (Fr(1), Fr(0)), 2, 3),
           (5, (Fr(1, 2), Fr(1, 2)), (Fr(1), Fr(0)), 1, 3),
           (5, (Fr(1, 2), Fr(1, 2)), (Fr(3, 2), Fr(1, 2)), 3, 2),
           (2, (Fr(0), Fr(1)), (Fr(1), Fr(0)), 1, 2),
           (2, (Fr(1), Fr(1)), (Fr(1), Fr(0)), 1, 2),
           (57, (Fr(1, 2), Fr(1, 2)), (Fr(1), Fr(0)), 3, 4)]


def delta_basis():
    out = []
    for a in range(3):
        for b in range(a, 3):
            el = {}
            for x, y in {(a, b), (b, a)}:
                el[("T", x, y)] = ONEF
                el[("T", x + 3, y + 3)] = ONEF
            out.append(el)
    return out


def check_Z():
    res, ok = [], True
    for t, om, c, q, N in CONFIGS:
        t0 = time.time()
        D_ = diagonal_member(t, om, c, q, N)
        FC, K = D_["FC"], D_["K"]
        row, good = {"t": t, "omega": str(om), "c": str(c), "q": q, "N": N}, True
        # (i) regular representation, eigenvectors, Theta_1 + Theta_2, rationality of B
        w2 = K.mul(om, om)                                   # omega^2 = u omega + v
        u_ = w2[1] / om[1]
        v_ = w2[0] - u_ * om[0]
        A = [[Fr(0), v_], [Fr(1), u_]]                       # multiplication by omega in the basis (1, omega)
        basis = [ONEF, om]
        Sc = [[K.tr(K.mul(c, K.mul(basis[i], basis[j]))) for j in range(2)] for i in range(2)]
        AtS = [[sum(A[k][i] * Sc[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        SA = [[sum(Sc[i][k] * A[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        selfadj = AtS == SA
        eig = True
        for omk in (om, K.conj(om)):
            vec = [ONEF, omk]                                # (1, omega_k) A = omega_k (1, omega_k)
            for j in range(2):
                lhs = (Fr(0), Fr(0))
                for i in range(2):
                    p = K.mul(vec[i], (A[i][j], Fr(0)))
                    lhs = (lhs[0] + p[0], lhs[1] + p[1])
                eig &= lhs == K.mul(omk, vec[j])
        Th = FC.add((ONEF, D_["T1"]), (ONEF, D_["T2"]))
        Th_want = FC.add(((K.tr(c), Fr(0)), rat(D_["Tp"])), ((2 * K.tr(K.mul(c, om)), Fr(0)), rat(D_["M"])),
                         ((K.tr(K.mul(c, w2)), Fr(0)), rat(D_["Tpp"])))
        dTh = FC.add((ONEF, Th), ((Fr(-1), Fr(0)), Th_want))
        ratB = all(not P[1] for P in D_["B"])
        g1 = selfadj and eig and not dTh[0] and not dTh[1] and ratB
        row["(i) self-adjoint, eigenvectors, Theta, B rational"] = g1
        good &= g1
        # (ii) K(B)
        B = D_["B"]
        dims = {"K(B)": FC.kdim(unit_combs(), B), "R": FC.kdim(unit_combs("T"), B),
                "G": FC.kdim(unit_combs("OL"), B)}
        g2 = dims == {"K(B)": 18, "R": 12, "G": 6}
        # (iii) K(ch E')
        chE = D_["chE"]
        dims.update({"K(chE')": FC.kdim(unit_combs(), [chE]), "K(chE') n T": FC.kdim(unit_combs("T"), [chE]),
                     "K(chE') n OL": FC.kdim(unit_combs("OL"), [chE])})
        chI = C6.block_class(3, {0: Fr(1), 2: -Fr(N, 2), 3: Fr(N, 3)})
        KI = kernel_Q(3, chI, SR.ht2(3))
        emb = [{(kd, a + s, b + s): (x, Fr(0)) for (kd, a, b), x in kv.items()} for s in (0, 3) for kv in KI]
        dims["dim K(ch I)"] = len(KI)
        dims["pi^* K(ch I) in K(chE')"] = FC.kdim(emb, [chE])
        dims["rank pi_1^* K(ch I) + pi_2^* K(ch I)"] = rank_F([dict(e) for e in emb], t)
        g3 = (dims["K(chE')"] == 18 and dims["K(chE') n T"] == 12 and dims["K(chE') n OL"] == 0
              and dims["dim K(ch I)"] == 9 and dims["pi^* K(ch I) in K(chE')"] == 18
              and dims["rank pi_1^* K(ch I) + pi_2^* K(ch I)"] == 18)
        # (iv) intersections
        dims["R n K(chE')"] = FC.kdim(unit_combs("T"), B + [chE])
        dims["K(B) n K(chE')"] = FC.kdim(unit_combs(), B + [chE])
        Dl = delta_basis()
        dims["Delta in K(B) n K(chE')"] = FC.kdim(Dl, B + [chE])
        dims["rank Delta"] = rank_F([dict(e) for e in Dl], t)
        g4 = (dims["R n K(chE')"] == 6 and dims["K(B) n K(chE')"] == 6 and dims["Delta in K(B) n K(chE')"] == 6
              and dims["rank Delta"] == 6)
        # (v) Markman's diagonal D
        d = N - 1
        chF = C6.block_class(3, {0: Fr(1), 1: Fr(1), 2: -Fr(d, 2), 3: -Fr(d, 6)})
        KF = kernel_Q(3, chF, SR.ht2(3))
        sgn = {"O": -1, "T": 1, "L": -1}
        Dm = []
        for kv in KF:
            el = {}
            for (kd, a, b), x in kv.items():
                el[(kd, a, b)] = (x, Fr(0))
                el[(kd, a + 3, b + 3)] = (sgn[kd] * x, Fr(0))
            Dm.append(el)
        lin = [1, 1, -Fr(d, 2), -Fr(d, 6)]
        chFF = FC.mul(FC.poly(lin, rat(D_["Tp"])), FC.poly(lin, rat(D_["Tpp"])))
        L = SS.add((1, D_["Tp"]), (1, D_["Tpp"]))
        eL = FC.poly([Fr(1, math.factorial(k)) for k in range(7)], rat(L))
        emL = FC.poly([Fr((-1) ** k, math.factorial(k)) for k in range(7)], rat(L))
        tw = FC.add((ONEF, FC.mul(emL, chFF)), ((Fr(-1), Fr(0)), chE))
        dims["dim D"] = len(Dm)
        dims["D in K(ch F[x]F)"] = FC.kdim(Dm, [chFF])
        dims["D_E' n K(B)"] = FC.kdim(Dm, [FC.mul(eL, b) for b in B])
        dims["rank D + Delta"] = rank_F([dict(e) for e in Dm + Dl], t)
        g5 = (not tw[0] and not tw[1] and dims["dim D"] == 9 and dims["D in K(ch F[x]F)"] == 9
              and dims["D_E' n K(B)"] == 6 and dims["rank D + Delta"] == 9)
        # (vi) ch(E') not in the span of Theta_1^i Theta_2^j
        mons = []
        P1 = ({0: Fr(1)}, {})
        for i in range(4):
            P2 = ({0: Fr(1)}, {})
            for j in range(4):
                mons.append(FC.mul(P1, P2))
                P2 = FC.mul(P2, D_["T2"])
            P1 = FC.mul(P1, D_["T1"])

        def as_col(P):
            col = {k: (v, Fr(0)) for k, v in P[0].items()}
            for k, v in P[1].items():
                col[k] = (col.get(k, (Fr(0), Fr(0)))[0], v)
            return col
        rm = rank_F([as_col(P) for P in mons], t)
        rme = rank_F([as_col(P) for P in mons + [chE]], t)
        dims["rank monomials, + ch E'"] = (rm, rme)
        g6 = rm == 16 and rme == 17
        good &= g2 and g3 and g4 and g5 and g6
        row["dims"] = dims
        row["ok"] = good
        ok &= good
        res.append(row)
        print("    (Z) t=%d omega=%s c=%s q=%s N=%d: %s; %s (%.1fs): %s"
              % (t, om, c, q, N, g1, dims, time.time() - t0, good), flush=True)
    return res, ok


# ----------------------------------------------------------------------------- (AA)
def biv_cols_F(P, n=6):
    cols = {}
    for i, j in itertools.combinations(range(4 * n), 2):
        col = {k: (v, Fr(0)) for k, v in biv(i, j, P[0], n).items()}
        for k, v in biv(i, j, P[1], n).items():
            col[k] = (col.get(k, (Fr(0), Fr(0)))[0], v)
        cols[(i, j)] = col
    return cols


def sub_cols(cols, vecs, K):
    out = []
    for v, w in itertools.combinations(vecs, 2):
        col = {}
        for i, vi in v.items():
            for j, wj in w.items():
                if i == j:
                    continue
                cf = K.mul(vi, wj)
                key = (i, j)
                if i > j:
                    key, cf = (j, i), (-cf[0], -cf[1])
                for k, x in cols[key].items():
                    u = K.mul(cf, x)
                    o = col.get(k, (Fr(0), Fr(0)))
                    col[k] = (o[0] + u[0], o[1] + u[1])
        out.append({k: x for k, x in col.items() if x != (0, 0)})
    return out


def eigen_vectors(om, K, k, n=6):
    """basis of V_k (eigen-block k) in the product frame: du_j, dubar_j and their duals."""
    omk = om if k == 1 else K.conj(om)
    omo = K.conj(om) if k == 1 else om
    dinv = K.inv((omo[0] - omk[0], omo[1] - omk[1]))
    a, b = K.mul(omo, dinv), K.mul((Fr(-1), Fr(0)), dinv)
    vecs = []
    for j in range(3):
        vecs += [{j: ONEF, j + 3: omk}, {n + j: ONEF, n + j + 3: omk},
                 {2 * n + j: a, 2 * n + j + 3: b}, {3 * n + j: a, 3 * n + j + 3: b}]
    return vecs


def product_vectors(coords, n=6):
    return [{l: ONEF} for j in coords for l in (j, n + j, 2 * n + j, 3 * n + j)]


def stab_F(P, om, K, t):
    cols = biv_cols_F(P)
    res = {"stab": 276 - rank_F(list(cols.values()), t)}
    for nm, vecs in (("so(V'_1) product", product_vectors(range(3))),
                     ("so(V'_2) product", product_vectors(range(3, 6))),
                     ("so(V_1) eigen", eigen_vectors(om, K, 1)), ("so(V_2) eigen", eigen_vectors(om, K, 2))):
        sc = sub_cols(cols, vecs, K)
        res["stab n " + nm] = len(sc) - rank_F(sc, t)
    return res


def fzero(P):
    return not P[0] and not P[1]


MINUS = (Fr(-1), Fr(0))


def check_AA():
    res, ok = {}, True
    t, om, _, _, N = CONFIGS[0]
    d = N - 1
    want_tw = {"stab": 70, "stab n so(V'_1) product": 0, "stab n so(V'_2) product": 0, "stab n so(V_1) eigen": 35,
               "stab n so(V_2) eigen": 35}
    want_E = {"stab": 70, "stab n so(V'_1) product": 35, "stab n so(V'_2) product": 35, "stab n so(V_1) eigen": 0,
              "stab n so(V_2) eigen": 0}
    ctw0, good1 = None, True
    for gam in (Fr(1), Fr(2), Fr(1, 3)):
        q = d / gam ** 2                                     # Theta_k = gamma phi_k^* theta_0, q = d/gamma^2
        D_ = diagonal_member(t, om, (gam, Fr(0)), q, N)
        FC, Bc = D_["FC"], D_["B"]
        ctw = FC.mul(FC.poly(D_["f"], D_["th1"]), FC.poly(D_["f"], D_["th2"]))
        if ctw0 is None:
            ctw0 = ctw
        Th = FC.add((ONEF, D_["T1"]), (ONEF, D_["T2"]))
        lhs = FC.mul(FC.poly([Fr(1, math.factorial(k)) / gam ** k for k in range(7)], Th), ctw)  # e^{Theta/gamma}
        ig = (1 / gam, Fr(0))
        cone = FC.mul(FC.add((ONEF, D_["C1"]), (ig, D_["S1"])), FC.add((ONEF, D_["C2"]), (ig, D_["S2"])))
        old = FC.mul(FC.add((ONEF, D_["C1"]), (ONEF, D_["S1"])), FC.add((ONEF, D_["C2"]), (ONEF, D_["S2"])))
        r_, qv, g_ = Fr(1), gam ** -2, 1 / gam              # coordinates of cone in B
        coords = FC.add((ONEF, Bc[0]), ((qv, Fr(0)), Bc[1]), ((g_, Fr(0)), Bc[2]))
        facts = {"c_tw rational": not ctw[1], "c_tw independent of gamma": ctw == ctw0,
                 "c_tw != ch E'": ctw != D_["chE"],
                 "e^{Theta/gamma} c_tw = (C1 + S1/gamma)(C2 + S2/gamma)": fzero(FC.add((ONEF, lhs), (MINUS, cone))),
                 "(C1 + S1/gamma)(C2 + S2/gamma) = C1 C2 + gamma^-2 S1 S2 + gamma^-1 (S1 C2 + C1 S2)":
                     fzero(FC.add((ONEF, cone), (MINUS, coords))),
                 "q r v = N(g)": q * r_ * (qv / q) == g_ * g_}
        old_holds = fzero(FC.add((ONEF, lhs), (MINUS, old)))
        facts["e^{Theta/gamma} c_tw = (C1 + S1)(C2 + S2) (only for gamma = 1)"] = old_holds
        g1 = all(v for k, v in facts.items() if "only for" not in k) and old_holds == (gam == 1)
        st = {}
        if gam == 1:
            for nm, P in (("c_tw", ctw), ("C1 C2", Bc[0]), ("ch E'", D_["chE"])):
                st[nm] = stab_F(P, om, D_["K"], t)
            g1 &= st["c_tw"] == want_tw and st["C1 C2"] == want_tw and st["ch E'"] == want_E
        else:
            st["(C1 + S1/gamma)(C2 + S2/gamma)"] = stab_F(cone, om, D_["K"], t)
            g1 &= st["(C1 + S1/gamma)(C2 + S2/gamma)"] == want_tw
        good1 &= g1
        res["diagonal member t=5, gamma=%s" % gam] = {"q": str(q), "facts": facts, "stabilizers": st, "ok": g1}
        print("    (AA) t = 5, N = 3, gamma = %s, q = %s: %s; stabilizers %s: %s" % (gam, q, facts, st, g1),
              flush=True)
    ok &= good1
    # (ii) integrality on H_1(X_0, Z) (x) O_F with E = J
    rows, good2 = [], True
    for tt in (2, 3, 5, 13, 57):
        Pw, s0, Pf, top, St = C6.lattice6(tt, [SS.Fe(tt, 1)] * 3)
        dF = tt if tt % 4 == 1 else 4 * tt
        eneg = {}
        for k in range(7):
            eneg = SS.add((1, eneg), (Fr((-1) ** k, math.factorial(k)), Pw[(k, 0)]))
        Tp = {(1 << (4 * b)) | (1 << (4 * b + 2)): Fr(1) for b in range(3)}
        Tpp = {(1 << (4 * b + 1)) | (1 << (4 * b + 3)): Fr(1) for b in range(3)}
        for NN in (2, 3, 4, 5):
            qq = Fr(NN - 1)
            cl = C6.bclasses6(tt, qq, Pw)
            v = 1 / qq
            ctw_l = SS.wmul(eneg, SS.lin(cl, [(1 - v) / 2, 1, (1 + v) / 2, 0]))
            fT = lambda T: SS.add((1, {0: Fr(1)}), (-Fr(NN, 2), SS.wmul(T, T)),
                                  (Fr(NN, 3), SS.wmul(T, SS.wmul(T, T))))
            chE_l = SS.wmul(fT(Tp), fT(Tpp))
            integ = all(Fr(x).denominator == 1 for x in ctw_l.values())
            integE = all(Fr(x).denominator == 1 for x in chE_l.values())
            rank1 = ctw_l.get(0, 0) == 1
            c1zero = not any(SS.popc(k) == 2 for k in ctw_l)
            muk = s0 * SS.wmul(SS.tau(ctw_l), ctw_l).get(top, 0)
            intc = s0 * ctw_l.get(top, 0)
            intE = s0 * chE_l.get(top, 0)
            r_ok = integ and integE and rank1 and c1zero and muk == 0 and intE == 4 * NN * NN and \
                intc == dF ** 3 * 4 * NN * NN and Pf == dF ** 3
            good2 &= r_ok
            rows.append({"t": tt, "N": NN, "Pf": str(Pf), "integral": integ, "chE' integral": integE,
                         "rank 1, c_1 = 0": rank1 and c1zero, "Mukai square": str(muk), "int c_tw": str(intc),
                         "int ch E'": str(intE), "ok": r_ok})
    ok &= good2
    res["integrality"] = rows
    print("    (AA) e^-Theta (C_1 + S_1)(C_2 + S_2) on H_1(X_0, Z) (x) O_F (E = J), t = 2, 3, 5, 13, 57, N = 2..5: "
          "integral, rank 1, c_1 = 0, Mukai square 0, int = d_F^3 (2N)^2: %s" % good2, flush=True)
    for r in rows:
        if r["N"] == 2:
            print("        t=%d Pf=%s int c_tw=%s int ch E'=%s" % (r["t"], r["Pf"], r["int c_tw"], r["int ch E'"]),
                  flush=True)
    # (iii) the proof of integrality: c_tw = phi_1^* ch(I) phi_2^* ch(I) on H^*(X_0 (x) O, Z), many orders and N
    rows3, good3 = [], True
    for idx, (tt, w, nm) in enumerate(ORDERS_INT):
        r3 = check_order_integral(tt, w, nm, gammas=idx in GAMMA_ORDERS)
        good3 &= r3["ok"]
        rows3.append(r3)
        short = {k: v for k, v in r3.items() if k not in ("N rows", "order")}
        print("        %s (t = %d, disc %s): %s" % (nm, tt, r3["disc(O)"], short), flush=True)
    ok &= good3
    res["integrality via phi_k^* (proof of 9.22(b))"] = rows3
    print("    (AA) c_tw = phi_1^* ch(I) phi_2^* ch(I) on H^*(X_0 (x) O, Z), %d Z-bases of orders, N = 1..9: "
          "phi_1^* ch(I) in H^* (x) O, c_tw rational and integral, rank 1, c_1 = 0, Mukai square 0, "
          "int = disc(O)^3 (2N)^2; gamma = 1, 2, 1/3 formula: %s" % (len(ORDERS_INT), good3), flush=True)
    return res, ok


ORDERS_INT = [(2, (Fr(0), Fr(1)), "Z[sqrt 2]"), (2, (Fr(1), Fr(1)), "Z[sqrt 2], basis (1, 1 + sqrt 2)"),
              (3, (Fr(0), Fr(1)), "Z[sqrt 3]"), (5, (Fr(1, 2), Fr(1, 2)), "Z[phi]"),
              (5, (Fr(0), Fr(1)), "Z[sqrt 5]"), (5, (Fr(3, 2), Fr(3, 2)), "Z[3 phi]"),
              (10, (Fr(0), Fr(1)), "Z[sqrt 10]"),
              (13, (Fr(1, 2), Fr(1, 2)), "Z[(1 + sqrt 13)/2]"), (13, (Fr(1), Fr(1)), "Z[1 + sqrt 13]"),
              (21, (Fr(1, 2), Fr(1, 2)), "Z[(1 + sqrt 21)/2]"), (57, (Fr(1, 2), Fr(1, 2)), "Z[(1 + sqrt 57)/2]")]
GAMMA_ORDERS = (0, 3, 5, 10)


def check_order_integral(t, om, name, gammas=False):
    """X = X_0 (x) O = X_0 x X_0, O = Z + Z omega; H^1(X, Z) has the basis pr_1^* x_j (bit j), pr_1^* y_j
    (bit 3 + j), pr_2^* x_j (bit 6 + j), pr_2^* y_j (bit 9 + j), theta_0 = sum_j x_j y_j principal; monomials are an
    integral basis of H^*(X, Z).  phi_k^* x = pr_1^* x + omega_k pr_2^* x on H^1(X_0) (phi_k = pi_1 +
    omega_k pi_2)."""
    FC, K = FClasses(t), Fld(t)
    one, top = ({0: Fr(1)}, {}), (1 << 12) - 1

    def gen(b):
        return ({1 << b: Fr(1)}, {})

    def pull_theta(w):
        th = ({}, {})
        for j in range(3):
            th = FC.add((ONEF, th), (ONEF, FC.mul(FC.add((ONEF, gen(j)), (w, gen(6 + j))),
                                                  FC.add((ONEF, gen(3 + j)), (w, gen(9 + j))))))
        return th
    a = [pull_theta(om), pull_theta(K.conj(om))]
    pw = []
    for ak in a:
        P_ = [one]
        for _ in range(4):
            P_.append(FC.mul(P_[-1], ak))
        pw.append(P_)
    a4zero = all(fzero(P_[4]) for P_ in pw)
    prod = {(i, j): FC.mul(pw[0][i], pw[1][j]) for i in (0, 2, 3) for j in (0, 2, 3)}
    tp = SS.add(*[(1, SS.wmul({1 << j: Fr(1)}, {1 << (3 + j): Fr(1)})) for j in range(3)])
    tpp = SS.add(*[(1, SS.wmul({1 << (6 + j): Fr(1)}, {1 << (9 + j): Fr(1)})) for j in range(3)])
    cube = lambda T: SS.wmul(T, SS.wmul(T, T))
    s0 = SS.wmul(cube(tp), cube(tpp)).get(top, 0) / 36            # int (theta'^3/6)(theta''^3/6) = 1
    disc = 4 * om[1] ** 2 * t                                      # (omega_1 - omega_2)^2
    row, good = {"order": name, "disc(O)": str(disc), "(phi_k^* theta_0)^4 = 0": a4zero, "N rows": []}, a4zero
    good &= s0 in (1, -1)
    for N in range(1, 10):
        f = {0: Fr(1), 2: -Fr(N, 2), 3: Fr(N, 3)}
        z1 = FC.add(*[((f[i], Fr(0)), pw[0][i]) for i in f])     # phi_1^* ch(I)
        inO = True                                                  # u + v sqrt t = a + b omega, a, b in Z
        for k in set(z1[0]) | set(z1[1]):
            u, v = z1[0].get(k, Fr(0)), z1[1].get(k, Fr(0))
            b = v / om[1]
            inO &= b.denominator == 1 and (u - b * om[0]).denominator == 1
        ctw = FC.add(*[((f[i] * f[j], Fr(0)), prod[(i, j)]) for i in f for j in f])
        rational = not ctw[1]
        c = ctw[0]
        integral = rational and all(x.denominator == 1 for x in c.values())
        rank1 = c.get(0, 0) == 1 and not any(SS.popc(k) == 2 for k in c)
        muk = s0 * SS.wmul(SS.tau(c), c).get(top, 0)
        intc = s0 * c.get(top, 0)
        fT = lambda T: SS.add((1, {0: Fr(1)}), (-Fr(N, 2), SS.wmul(T, T)), (Fr(N, 3), cube(T)))
        chE = SS.wmul(fT(tp), fT(tpp))
        intE = s0 * chE.get(top, 0)
        chE_ok = all(x.denominator == 1 for x in chE.values()) and intE == 4 * N * N and \
            SS.wmul(SS.tau(chE), chE).get(top, 0) == 0
        r_ok = inO and integral and rank1 and muk == 0 and intc == disc ** 3 * 4 * N * N and chE_ok
        good &= r_ok
        row["N rows"].append({"N": N, "phi_1^* ch(I) in H^* (x) O": inO, "rational": rational, "integral": integral,
                              "rank 1, c_1 = 0": rank1, "Mukai square": str(muk), "int c_tw": str(intc),
                              "ch E' integral, int 4N^2, Mukai square 0": chE_ok, "ok": r_ok})
    row["N = 1..9 ok"] = all(r["ok"] for r in row["N rows"])
    if gammas:                                                     # NOTE 9.22(b) for rational gamma, N = 3
        N, gm = 3, {}
        fl = [Fr(1), Fr(0), -Fr(N, 2), Fr(N, 3)]
        ctw = FC.mul(FC.poly(fl, a[0]), FC.poly(fl, a[1]))
        for gam in (Fr(1), Fr(2), Fr(1, 3)):
            q = (N - 1) / gam ** 2
            new, old = one, one
            for k in range(2):
                Tk = FC.add(((gam, Fr(0)), a[k]))                 # Theta_k = gamma phi_k^* theta_0
                Ek = FC.poly([Fr((-1) ** m, math.factorial(m)) for m in range(7)], a[k])    # e^{-Theta_k/gamma}
                Ck, Sk = FC.poly([1, 0, -q / 2], Tk), FC.poly([0, 1, 0, -q / 6], Tk)
                new = FC.mul(new, FC.mul(Ek, FC.add((ONEF, Ck), ((1 / gam, Fr(0)), Sk))))
                old = FC.mul(old, FC.mul(Ek, FC.add((ONEF, Ck), (ONEF, Sk))))
            hn, ho = fzero(FC.add((ONEF, ctw), (MINUS, new))), fzero(FC.add((ONEF, ctw), (MINUS, old)))
            gm[str(gam)] = {"c_tw = e^{-Theta/gamma} (C1 + S1/gamma)(C2 + S2/gamma)": hn,
                            "c_tw = e^{-Theta/gamma} (C1 + S1)(C2 + S2)": ho}
            good &= hn and ho == (gam == 1)
        row["gamma (N = 3)"] = gm
    row["ok"] = good
    return row


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "split.json"
    res, allok, secs = {}, True, {}
    for key, fn in (("Y", check_Y), ("Z", check_Z), ("AA", check_AA)):
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
