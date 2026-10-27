"""check_product.py -- exact checks of the product argument of NOTE.md Section 3 (Theorem 3.4).

For a free lattice Lambda = O_F^4 with theta = delta g1 (e12 + g2 e34) (shape_search.twoform) the decomposition
V = V_A + V_B (V_A = <e1, e2>_F, V_B = <e3, e4>_F) is theta-orthogonal and Lambda = Lambda_A + Lambda_B.  For the
surface lattices: W_A = |int theta_A,1 theta_A,2| = |N(g1)|, J_A = {g : Tr(g theta_A) integral} = g1^-1 O_F; for B the
same with g1 g2.  H_{A,Z} = {c00 in Z, c10 in J_A, c01 = sigma(c10), c11 in (1/W_A) Z} (in the e-coordinates of
check_sym2: e_ab <-> theta_A,1^a theta_A,2^b).  Checks, for 9 free types per field and 40 random X in L_Herm:
 (P1) W = 4 W_A W_B;
 (P2) S(X) lies in H_{A,Z} (x) H_{B,Z}: the matrix S_AB = F_A^-1 S F_B^-T is integral;
 (P3) Pi(X) = (W^2/4) det(eps) = (W^2/4) det(S_AB) det(b(f^B_J, f^A_K)), and
      |(W^2/4) det(b(f^B, f^A))| = 4 d_F (W_A N(J_A)) (W_B N(J_B)), with W_A N(J_A) = W_B N(J_B) = 1 here,
      N(J) = covol(J)/covol(O_F).
Usage: python -I -B check_product.py      Expected last line "ALL CHECKS PASS: True".
"""
import importlib.util
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
TL = CS.TL                       # one copy of twist_lib, so that Q2 instances are compatible
SS, Q2, Herm = TL.SS, TL.Q2, TL.Herm

ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print("(%s) %s" % (name, bool(cond)), flush=True)


def q2(t, g):
    return Q2(t, Fr(g.u), Fr(g.v))


def surface_basis(t, g):
    """columns f_J (e-coordinates e00, e01, e10, e11) of H_{Z} for theta_A = delta g e12 on O_F^2."""
    gq = q2(t, g)
    WA = abs(gq.norm())
    om = Q2(t, Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else Q2(t, 0, 1)
    z0 = Q2(t, 0)
    js = [gq.inv(), om * gq.inv()]                     # Z-basis of J = g^-1 O_F
    cols = [(Q2(t, 1), z0, z0, z0), (z0, z0, z0, Q2(t, Fr(1) / WA))] + [(z0, j.conj(), j, z0) for j in js]
    NJ = Fr(1) / WA                                    # covol(g^-1 O_F)/covol(O_F) = 1/|N(g)|
    return [[cols[J][I] for J in range(4)] for I in range(4)], WA, NJ


def main():
    rng = random.Random(20261026)
    p1 = p2 = p3 = True
    for t in (5, 2, 3, 13):
        gs = SS.reps(t, 9)[:3]
        b = CS.bmat(t)
        for g1 in gs:
            for g2 in gs:
                L = TL.free_lattice(t, g1, g2)
                FA, WA, NJA = surface_basis(t, g1)
                FB, WB, NJB = surface_basis(t, g1 * g2)
                p1 &= L.W == 4 * WA * WB
                FAi, FBi = CS.inv4(FA), CS.inv4(FB)
                G = CS.mm(CS.tr_(FB), CS.mm(b, FA))          # G[J][K] = b(f^B_J, f^A_K)
                dG = CS.det4(G)
                p3 &= dG.is_rational() and abs(L.W ** 2 / 4 * dG.u) == 4 * CS.dF(t) * (WA * NJA) * (WB * NJB)
                p3 &= WA * NJA == 1 and WB * NJB == 1
                for _ in range(40):
                    m, a, k = (CS.comb(rng, B) for B in (L.L4, L.J2, L.J6))
                    X = Herm(t, rng.randint(-3, 3), Q2(t, a[0], a[1]), Q2(t, m[0], m[1]), m[2], Q2(t, k[0], k[1]),
                             Fr(rng.randint(-5, 5), L.W))
                    assert L.integral(X)
                    S = CS.to_S(X)
                    SAB = CS.mm(FAi, CS.mm(S, CS.tr_(FBi)))
                    p2 &= all(e.is_rational() and e.u.denominator == 1 for row in SAB for e in row)
                    Pi = L.W ** 2 / 4 * CS.det4(CS.mm(S, b)).u
                    p3 &= Pi == L.W ** 2 / 4 * (CS.det4(SAB) * dG).u
                print("    t=%2d %-36s W=%-6s W_A=%-4s W_B=%-5s det b(f^B,f^A)=%s" %
                      (t, L.label, L.W, WA, WB, dG), flush=True)
    check("P1", p1)
    check("P2", p2)
    check("P3", p3)
    print("ALL CHECKS PASS:", ok)


if __name__ == "__main__":
    main()
