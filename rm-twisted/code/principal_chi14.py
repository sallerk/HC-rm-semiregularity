"""principal_chi14.py -- all nonzero integral classes of B with chi <= 14 on the principal free lattice, q = 1
(LaTeX note on the semiregularity step, remark after Corollary 5.7; referee R1, point 5(iii)).

Lattice: O_F^4 with E = delta (e_01 + e_23), delta a generator of the inverse different (Theta = Tr E principal),
built by twist_lib.free_lattice(t, 1, 1); the enumeration is shape_search.generic_classes(..., generic=False),
exact over the two rank-2 lattices of even and odd parts.  For each class: (x, y, z, w), r = x + z, g = y + w sqrt t
and whether g is totally nonnegative (exact: y >= 0 and y^2 >= t w^2), totally nonpositive, or indefinite.
Expected: for t = 5 the classes +-alpha, +-beta (chi = 8) and the four classes +-(beta +- sqrt5 beta~)/2
(x = z = 0, y = +-1/2, w = +-1/2, chi = 12, g indefinite); for t = 2, 3, 13 only +-alpha, +-beta.
Usage: python -I -B principal_chi14.py      Expected last line "ALL CHECKS PASS: True".
"""
import importlib.util
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8")
TLP = __file__.replace("principal_chi14.py", "twist_lib.py")
spec = importlib.util.spec_from_file_location("twist_lib", TLP)
TL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TL)
SS = TL.SS


def sign_type(t, y, w):
    if y >= 0 and y * y >= t * w * w:
        return "g totally >= 0"
    if y <= 0 and y * y >= t * w * w:
        return "g totally <= 0"
    return "g indefinite"


H = Fr(1, 2)
BASE = {(8, (1, 0, 0, 0)), (8, (-1, 0, 0, 0)), (8, (0, 1, 0, 0)), (8, (0, -1, 0, 0))}
EXTRA5 = {(12, (0, a, 0, b)) for a in (H, -H) for b in (H, -H)}
ok = True
for t in (5, 2, 3, 13):
    one = SS.Fe(t, 1)
    L = TL.free_lattice(t, one, one)
    out = SS.generic_classes(t, 1, L.P, L.s0, 14, generic=False)
    found = {(int(chi), tuple(Fr(u) for u in co)) for chi, co, _ in out}
    want = BASE | (EXTRA5 if t == 5 else set())
    good = found == want and len(out) == len(found)
    ok &= good
    print("t = %d, q = 1, principal free lattice: %d classes with chi <= 14; as expected: %s" % (t, len(out), good))
    for chi, (x, y, z, w) in sorted(found):
        print("  chi=%d  (x, y, z, w) = (%s, %s, %s, %s)  r = %s  g = %s + %s sqrt%d  (%s)"
              % (chi, x, y, z, w, x + z, y, w, t, sign_type(t, y, w)))
    if t == 5:
        ind = all(sign_type(t, co[1], co[3]) == "g indefinite" for chi, co in found if chi == 12)
        ok &= ind
        print("  the chi = 12 classes have r = 0 and g indefinite (c_1 = L_g not nef): %s" % ind)
print("ALL CHECKS PASS: %s" % ok)
