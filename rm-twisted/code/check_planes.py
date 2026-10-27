"""check_planes.py -- which eigenspace of K carries which plane (NOTE.md Proposition 3.1; referee point 2(a)).

For random twisted secant classes X = e^{-L_h} c (c = x alpha + y beta + z alpha~ + w sqrt t beta~), with
J_p = p om (p = (1, -h, (h^2 + q)/2) read as the symmetric matrix [[p0, p1], [p1, 2 p2]]), K = J_p (x) J_sigma(p)
(K^2 = q^2) and eps = S(X) b, checks exactly:
 (K1) eps K = K eps;
 (K2) eps^2 (K + q) = m_P~ (K + q) and eps^2 (K - q) = m_P (K - q), m_P = 2 chi(c_P)/W, m_P~ = 2 chi(c_P~)/W,
      i.e. eps^2 = 2 chi(c_P~)/W on H_+ = {K = q} and 2 chi(c_P)/W on H_- = {K = -q};
 (K3) for c in P (z = w = 0): eps (K + q) = 0 (eps vanishes on H_+); for c in P~ (x = y = 0): eps (K - q) = 0.
Usage: python -I -B check_planes.py      Expected last line "ALL CHECKS PASS: True"."""
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
TL, Q2 = CS.TL, CS.Q2
OM = CS.OM


def kmat(t, q, h):
    p = [Q2(t, 1), -h, (h * h + Q2(t, q)) * Fr(1, 2)]
    pm = [[p[0], p[1]], [p[1], p[2] * 2]]
    ps = [[e.conj() for e in row] for row in pm]
    J1 = [[sum((pm[a][k] * OM[k][b] for k in range(2)), Q2(t, 0)) for b in range(2)] for a in range(2)]
    J2 = [[sum((ps[a][k] * OM[k][b] for k in range(2)), Q2(t, 0)) for b in range(2)] for a in range(2)]
    return [[J1[a][a2] * J2[b][b2] for (a2, b2) in CS.IDX] for (a, b) in CS.IDX]


def add_scalar(A, s):
    return [[A[i][j] + (s if i == j else 0) for j in range(4)] for i in range(4)]


def scal(A, s):
    return [[A[i][j] * s for j in range(4)] for i in range(4)]


def main():
    rng = random.Random(20261024)
    k1 = k2 = k3 = True
    for t in (5, 2, 3, 13):
        b = CS.bmat(t)
        for trial in range(30):
            q = Fr(rng.randint(1, 9), rng.randint(1, 6))
            x, y, z, w = (CS.rfr(rng) for _ in range(4))
            if trial % 3 == 1:
                z = w = Fr(0)
            elif trial % 3 == 2:
                x = y = Fr(0)
            h = Q2(t, CS.rfr(rng, 4, 3), CS.rfr(rng, 4, 3))
            X = TL.twist(TL.secant(t, q, x, y, z, w), h)
            eps = CS.mm(CS.to_S(X), b)
            K = kmat(t, q, h)
            KK = CS.mm(K, K)
            assert all(KK[i][j] == (q * q if i == j else 0) for i in range(4) for j in range(4))
            k1 &= CS.mm(eps, K) == CS.mm(K, eps)
            mP, mPt = 4 * q * (q * x * x + y * y), 4 * q * (q * z * z + t * w * w)   # 2 chi/W for the two planes
            e2 = CS.mm(eps, eps)
            Kp, Km = add_scalar(K, q), add_scalar(K, -q)
            k2 &= CS.mm(e2, Kp) == scal(Kp, mPt) and CS.mm(e2, Km) == scal(Km, mP)
            zero = [[Q2(t, 0)] * 4 for _ in range(4)]
            if trial % 3 == 1:
                k3 &= CS.mm(eps, Kp) == zero
            if trial % 3 == 2:
                k3 &= CS.mm(eps, Km) == zero
    ok = all((k1, k2, k3))
    print("K1 eps commutes with K:", k1)
    print("K2 eps^2 = 2 chi(c_P~)/W on H_+, 2 chi(c_P)/W on H_-:", k2)
    print("K3 eps vanishes on H_+ for c in P, on H_- for c in P~:", k3)
    print("ALL CHECKS PASS:", ok)


if __name__ == "__main__":
    main()
