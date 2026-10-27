"""multiplier.py -- content of Pi versus the multiplier order, on random lattices (data for NOTE.md Section 3.5).

For Lambda_Y = O_F^4 + sum Z gamma_i / N (shape_search.random_lattice) inside the free type (g1, g2): the conductor f
of the multiplier order O(Lambda_Y) = Z + f O_F (least f with f omega Lambda_Y in Lambda_Y, omega the generator of
O_F), and the content of Pi on L_Herm(Lambda_Y) (exact).  Tests the guess content in 4 disc(O) Z, disc(O) = d_F f^2.
Usage: python -I -B multiplier.py NLAT"""
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


def conductor(t, SY, N):
    n = len(SY)
    if t % 4 == 1:
        Om = [[(Fr(int(i == j)) + SY[i][j]) / 2 for j in range(n)] for i in range(n)]
    else:
        Om = SY
    for f in range(1, N ** 4 + 1):
        if all((f * x).denominator == 1 for row in Om for x in row):
            return f
    return None


def main():
    nlat = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    rng = random.Random(20261028)
    bad = 0
    for t in (5, 2, 3, 13):
        gs = SS.reps(t, 9)[:3]
        for N in (2, 3, 4):
            for i in range(nlat):
                g1, g2 = rng.choice(gs), rng.choice(gs)
                A = SS.random_lattice(rng, N)
                M1, Ms, SY, P, s0, T12, W = SS.ylattice(t, g1, g2, A)
                f = conductor(t, SY, N)
                L = TL.Lat(t, P, s0, "rand N=%d #%d (%s,%s)" % (N, i, g1, g2))
                cont = CS.content_Pi(L)
                ratio = cont / (4 * CS.dF(t) * f * f)
                bad += ratio.denominator != 1
                print("t=%2d %-36s W=%-9s f=%d  content/(4 d_F)=%s  content/(4 disc O)=%s" %
                      (t, L.label, L.W, f, cont / (4 * CS.dF(t)), ratio), flush=True)
    print("lattices with content not in 4 disc(O) Z:", bad)


if __name__ == "__main__":
    main()
