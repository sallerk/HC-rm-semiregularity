"""bmb_m_all.py -- the identities of bmb_m.py on all of R, not only on twisted secant spaces (blind check R3, point 8).

A class of R is c = sum x_ij Theta_1^i Theta_2^j (i, j <= 2) with the 3x3 matrix X = (x_ij) over F = Q(sqrt t) Hermitian for
sigma (x_ji = sigma(x_ij)); all nine rational coordinates, t and varpi are symbols.  sqrt(t) is a symbol S, and every
expression is reduced modulo S^2 - t, so "== 0" is an exact polynomial identity.  As in the proof of Proposition 6.4 of
the note, eps(c) = sum x_ij J_i (x) J_j with J_0 = [[0,1],[0,0]], J_1 = [[-1,0],[0,1]], J_2 = [[0,0],[-2,0]] (the operators
J_{Theta^i} in the basis 1, a of F[a]/(a^2)), and Pi = (varpi^2/4) det eps(c).  chi(c) is Lemma 6.2(a).  [BMB]'s V = D X D,
D = diag(1, 1, 2) (V_ab = a! b! x_ab), J = antidiag(1, -2, 1), Vbar = sigma(V) entrywise, m = (varpi/4)^2 sigma_2(J V J Vbar).
Exact symbolic checks:
    J_i J_j + J_j J_i = -2 beta_0(e_i, e_j)  (the anticommutation rule used in the proof of Proposition 6.4);
    tr eps(c) = 0 and (varpi/4) tr eps(c)^2 = chi(c)  (Lemma 6.3);  Pi is rational (sigma-invariant);
    chi = (varpi/4) tr(J V J Vbar);   4 Pi = chi^2 - 4 m   (on all of R).
Usage: python -I -B bmb_m_all.py      Expected last line "ALL CHECKS PASS: True".
"""
import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
t, vp, S = sp.symbols('t varpi S')
a00, a11, a22, p1, p2, u1, u2, v1, v2 = sp.symbols('a00 a11 a22 p1 p2 u1 u2 v1 v2')


def red(e):
    """Normal form modulo S^2 - t."""
    return sp.expand(sp.rem(sp.expand(e), S ** 2 - t, S))


def sig(e):
    return red(sp.expand(e).subs(S, -S))


def is0(e):
    return red(e) == 0


x01, x02, x12 = p1 + p2 * S, u1 + u2 * S, v1 + v2 * S
X = sp.Matrix([[a00, x01, x02], [sig(x01), a11, x12], [sig(x02), sig(x12), a22]])
herm = all(is0(sig(X[i, j]) - X[j, i]) for i in range(3) for j in range(3))

Jop = [sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([[-1, 0], [0, 1]]), sp.Matrix([[0, 0], [-2, 0]])]
beta0 = sp.Matrix([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
anti = all((Jop[i] * Jop[j] + Jop[j] * Jop[i] + 2 * beta0[i, j] * sp.eye(2)) == sp.zeros(2, 2)
           for i in range(3) for j in range(3))

eps = sp.zeros(4, 4)
for i in range(3):
    for j in range(3):
        eps += X[i, j] * sp.kronecker_product(Jop[i], Jop[j])
eps = eps.applyfunc(red)
tr0 = is0(eps.trace())
chi = red(vp * (2 * a00 * a22 + a11 ** 2 + 2 * X[2, 0] * sig(X[2, 0])
                - 2 * (sig(X[1, 0]) * X[2, 1] + X[1, 0] * sig(X[2, 1]))))
c_eps = is0(vp / 4 * (eps * eps).trace() - chi)
Pi = red(vp ** 2 / 4 * eps.det(method='berkowitz'))
c_rat = is0(Pi - sig(Pi))

D = sp.diag(1, 1, 2)
V = D * X * D
Jb = sp.Matrix([[0, 0, 1], [0, -2, 0], [1, 0, 0]])
M = (Jb * V * Jb * V.applyfunc(sig)).applyfunc(red)
trM = red(M.trace())
s2 = red((trM ** 2 - (M * M).trace()) / 2)
m = red(vp ** 2 / 16 * s2)
c_chi = is0(vp / 4 * trM - chi)
c_pi = is0(4 * Pi - (chi ** 2 - 4 * m))
print("X Hermitian:", herm)
print("J_i J_j + J_j J_i = -2 beta_0:", anti)
print("tr eps = 0:", tr0)
print("(varpi/4) tr eps^2 = chi (Lemma 6.2(a)):", c_eps)
print("Pi rational:", c_rat)
print("chi = (varpi/4) tr(J V J Vbar) on all of R:", c_chi)
print("4 Pi = chi^2 - 4 m on all of R:", c_pi)
print("ALL CHECKS PASS:", herm and anti and tr0 and c_eps and c_rat and c_chi and c_pi)
