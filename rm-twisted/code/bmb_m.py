"""bmb_m.py -- the invariant m of [BMB] Proposition 19.36 against our Pi (NOTE.md Section 3.3, Remark after Theorem 3.4).

[BMB]'s proof of Proposition 19.36(ii) (v5.1.0, pp. 235-236) writes a class of R = span(Theta_1^a Theta_2^b) as
v = sum V_ab Theta_1^a Theta_2^b / (a! b!), so V_ab = a! b! x_ab, puts J = antidiag(1, -2, 1) and m = sigma_2(J V J Vbar),
Vbar = sigma(V) = V^T (Hermitian).  Here c = x alpha + y beta + z alpha~ + w sqrt(t) beta~ (q rational, the matrix of
Section 2), twisted by e^{-L_h} (h = h0 + h1 sqrt t).  Exact symbolic checks (t, q, h0, h1, varpi = int Theta_1^2
Theta_2^2 symbolic):
    chi = (varpi/4) tr(J V J V^T);   with m = (varpi/4)^2 sigma_2(J V J V^T):   4 Pi = chi^2 - 4 m,
    Pi = chi(c_P) chi(c_P~);   m = ((chi(c_P) - chi(c_P~))/2)^2;   det(J V J V^T) = 0.
Usage: python -I -B bmb_m.py      Expected last line "ALL CHECKS PASS: True".
"""
import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
x, y, z, w, q, t, h0, h1, vp = sp.symbols('x y z w q t h0 h1 varpi')
s = sp.sqrt(t)


def sigm(M):
    return M.applyfunc(lambda e: sp.expand(sp.expand(e).subs(sp.sqrt(t), -sp.sqrt(t))))


r, v, g = x + z, z - x, y + w * s
gb = y - w * s
X = sp.Matrix([[r, gb, -r * q / 2], [g, q * v, -q * g / 2], [-r * q / 2, -q * gb / 2, r * q ** 2 / 4]])
h = h0 + h1 * s
Eh = sp.Matrix([[1, 0, 0], [-h, 1, 0], [h ** 2 / 2, -h, 1]])
Xt = (Eh * X * sigm(Eh).T).applyfunc(sp.expand)
herm = sp.simplify(sigm(Xt) - Xt.T) == sp.zeros(3, 3)
D = sp.diag(1, 1, 2)
V = D * Xt * D
J = sp.Matrix([[0, 0, 1], [0, -2, 0], [1, 0, 0]])
M = (J * V * J * V.T).applyfunc(sp.expand)
trM = sp.expand(M.trace())
s2 = sp.expand((trM ** 2 - (M * M).trace()) / 2)
chi = 2 * vp * (q ** 2 * (x ** 2 + z ** 2) + q * (y ** 2 + t * w ** 2))
chiP = 2 * vp * q * (q * x ** 2 + y ** 2)
chiPt = 2 * vp * q * (q * z ** 2 + t * w ** 2)
c1 = sp.expand(vp / 4 * trM - chi) == 0
m = vp ** 2 / 16 * s2
c2 = sp.expand(4 * chiP * chiPt - (chi ** 2 - 4 * m)) == 0
c3 = sp.expand(m - ((chiP - chiPt) / 2) ** 2) == 0
c4 = sp.expand(M.det()) == 0
print("twisted matrix Hermitian:", herm)
print("chi = (varpi/4) tr(J V J V^T):", c1)
print("4 Pi = chi^2 - 4 m with m = (varpi/4)^2 sigma_2(J V J V^T), Pi = chi(c_P) chi(c_P~):", c2)
print("m = ((chi(c_P) - chi(c_P~))/2)^2:", c3)
print("det(J V J V^T) = 0:", c4)
print("ALL CHECKS PASS:", herm and c1 and c2 and c3 and c4)
