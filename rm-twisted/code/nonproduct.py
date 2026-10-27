"""nonproduct.py -- the content of Pi on lattices built to be non-products (data for NOTE.md Section 3.5).

Lambda = O^4 for the order O = Z + f O_F (conductor f), with theta = delta (e12 + e34 + omega e13).  In the coordinates
x' = (x1, x2 + omega x3, x3, x4) theta is the principal form delta (e'12 + e'34), and Lambda becomes the sublattice
Lambda' = {(x1, x2 + omega x3, x3, x4) : x_i in O} of O_F^4 (twist_lib.y_lattice with g1 = g2 = 1).  Since omega is
not in O, the obvious decompositions fail (whether some other theta-orthogonal decomposition splits Lambda is not
decided here).  Also: the same with omega e13 replaced by omega e13 + omega e24, and superlattices
Lambda' + Z gamma / f.  For each lattice: the multiplier ring index (via the action of sqrt t), W, and the content of
Pi on L_Herm divided by 4 d_F (exact, all of L_Herm).
Usage: python -I -B nonproduct.py"""
import importlib.util
import os
import random
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("check_sym2", os.path.join(HERE, "check_sym2.py"))
CS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CS)
TL = CS.TL
SS = TL.SS


def omega(t):
    return (Fr(1, 2), Fr(1, 2)) if t % 4 == 1 else (Fr(0), Fr(1))


def fe_coords(t, u, v):
    """coordinates of u + v sqrt t in the Z-basis (1, omega) of O_F."""
    if t % 4 == 1:            # omega = (1 + sqrt t)/2: u + v sqrt t = (u - v) + 2 v omega
        return (u - v, 2 * v)
    return (u, v)


def mul(t, a, b):
    return (a[0] * b[0] + t * a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def lattice_matrix(t, f, extra_e24=False):
    """Z-basis (columns, 8 coordinates in the basis (e_a, omega e_a) of O_F^4) of Lambda'."""
    w = omega(t)
    obasis = [(Fr(1), Fr(0)), (f * w[0], f * w[1])]           # Z-basis of O = Z + f O_F
    cols = []
    for a in range(4):                                        # x_a runs over O
        for o in obasis:
            vec = [(Fr(0), Fr(0))] * 4
            vec[a] = o
            if a == 2:                                        # x3 contributes omega x3 to the second coordinate
                vec[1] = mul(t, w, o)
            if extra_e24 and a == 3:                          # x4 contributes omega x4 to the first coordinate
                vec[0] = mul(t, w, o)
            cols.append([c for (u, v) in vec for c in fe_coords(t, u, v)])
    return [[cols[j][i] for j in range(8)] for i in range(8)]


def main():
    rng = random.Random(20261027)
    for t in (5, 2, 3, 13):
        one = SS.Fe(t, 1)
        for f in (2, 3):
            for e24 in (False, True):
                A = lattice_matrix(t, f, e24)
                L = TL.y_lattice(t, one, one, A, "O^4, f=%d%s" % (f, ", +e24" if e24 else ""))
                cont = CS.content_Pi(L)
                print("t=%2d %-24s W=%-8s content/(4 d_F) = %s" % (t, L.label, L.W, cont / (4 * CS.dF(t))),
                      flush=True)
            # a superlattice Lambda' + Z gamma / f of the first one
            A = lattice_matrix(t, f)
            gam = [Fr(rng.randrange(f), f) for _ in range(8)]
            B = [[sum(A[i][k] * x for k, x in enumerate(col)) for i in range(8)] for col in
                 [[int(k == j) for k in range(8)] for j in range(8)]]
            gens = [[A[i][j] for i in range(8)] for j in range(8)] + [[sum(A[i][k] * gam[k] for k in range(8))
                                                                        for i in range(8)]]
            D = 1
            for v in gens:
                for x in v:
                    D = D * x.denominator // __import__("math").gcd(D, x.denominator)
            H = SS.flint.fmpz_mat([[int(x * D) for x in v] for v in gens]).hnf()
            rows = [[Fr(int(H[r, c]), D) for c in range(8)] for r in range(H.nrows()) if any(H[r, c] != 0
                                                                                             for c in range(8))]
            A2 = [[rows[c][r] for c in range(8)] for r in range(8)]
            L = TL.y_lattice(t, one, one, A2, "O^4 + Z gamma/%d" % f)
            cont = CS.content_Pi(L)
            print("t=%2d %-24s W=%-8s content/(4 d_F) = %s" % (t, L.label, L.W, cont / (4 * CS.dF(t))), flush=True)


if __name__ == "__main__":
    main()
