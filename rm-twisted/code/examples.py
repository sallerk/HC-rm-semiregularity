"""examples.py -- the four sharp examples of Theorem 3.6 (NOTE.md), checked exactly.

Each is a generic integral twisted secant class c' = e^{-L_h} c, c = x alpha + y beta + z alpha~ + w sqrt t beta~ in
B(Theta, q), on a free (product) lattice of twist_scan2.py, with Theta normalized to the primitive integral class of
its ray as there.  Checked: c' integral (degree lattices and the slow path), chi(c') = the stated value, c' off P u P~
and off the cone, the projections chi(c_P), chi(c_P~), and Pi = chi(c_P) chi(c_P~) = (W^2/4) det eps(c') in 4 d_F Z;
the value of chi equals the bound s(d_F) = min{a + b : a, b >= 1, a + b even, 4 d_F | ab} of Theorem 3.6.
Usage: python -I -B examples.py      Expected last line "ALL CHECKS PASS: True"."""
import importlib.util
import os
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("check_sym2", os.path.join(HERE, "check_sym2.py"))
CS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CS)
TL = CS.TL
SS, Q2 = TL.SS, TL.Q2

EXAMPLES = [  # t, (g1, g2), q, h, (x, y, z, w), chi
    (5, ((Fr(-3, 2), Fr(-3, 2)), (Fr(-3, 2), Fr(-3, 2))), Fr(1, 9), (0, Fr(1, 3)), (Fr(-3, 2), 0, Fr(-5, 2), Fr(-1, 3)),
     12),
    (2, ((0, -2), (0, -1)), Fr(1, 8), (Fr(1, 2), Fr(1, 4)), (0, 1, 0, Fr(-1, 2)), 12),
    (3, ((0, -1), (1, 0)), Fr(1, 12), (Fr(1, 2), 0), (0, -1, 0, Fr(-2, 3)), 14),
    (13, ((Fr(-1, 2), Fr(-1, 2)), (1, 0)), Fr(1), (0, 0), (0, Fr(1, 6), 0, Fr(-1, 6)), 28),
]


def s_bound(dF):
    best = None
    for s in range(2, 200, 2):
        if any((a * (s - a)) % (4 * dF) == 0 for a in range(1, s)):
            return s
    return best


def main():
    ok = True
    for t, ((a1, b1), (a2, b2)), q, (h0, h1), (x, y, z, w), chi in EXAMPLES:
        g1, g2 = SS.Fe(t, a1, b1), SS.Fe(t, a2, b2)
        dl = SS.delta(t)
        P, s0 = SS.setup(t, dl * g1, dl * g1 * g2)
        cth = SS.content(P[(1, 0)])
        P = {(a, b): {k: v / cth ** (a + b) for k, v in P[(a, b)].items()} for (a, b) in P}
        L = TL.Lat(t, P, s0, "free(%s,%s)" % (g1, g2))
        x, y, z, w = (Fr(u) for u in (x, y, z, w))
        X = TL.twist(TL.secant(t, q, x, y, z, w), Q2(t, Fr(h0), Fr(h1)))
        fast = L.integral(X)
        slow = all(v.denominator == 1 for v in L.form(X).values())
        c = L.chi(X)
        W = L.W
        r, v, Ng = x + z, z - x, y * y - t * w * w
        n = r * q * v - Ng
        kind = TL.kind(t, q, x, y, z, w)
        cP, cPt = 2 * W * q * (q * x * x + y * y), 2 * W * q * (q * z * z + t * w * w)
        Pi = W * W / 4 * CS.det4(CS.mm(CS.to_S(X), CS.bmat(t))).u
        dF = CS.dF(t)
        good = (fast and slow and c == chi and kind == "generic" and cP + cPt == c and cP * cPt == Pi
                and Pi % (4 * dF) == 0 and s_bound(dF) == chi and (c / 2 - W * q * n) == cP)
        ok &= good
        print("t=%2d %-30s W=%-4s q=%-5s h=(%s,%s) u=(%s,%s,%s,%s): integral %s/%s, chi=%s, kind %s, "
              "(chi_P, chi_P~)=(%s,%s), Pi=%s=4 d_F*%s, s(d_F)=%s: %s"
              % (t, L.label, W, q, h0, h1, x, y, z, w, fast, slow, c, kind, cP, cPt, Pi, Pi / (4 * dF), s_bound(dF),
                 good), flush=True)
    print("ALL CHECKS PASS:", ok)


if __name__ == "__main__":
    main()
