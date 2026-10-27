"""p88_check.py -- exact checks for the Phase 88 statements (WORK-2026-10-26.md).

 (A) Degree four: for every integral class Q_m + n Theta_1 Theta_2 (Q_m = s1(m) Theta_1^2 + s2(m) Theta_2^2),
     W (4 N(m) - n^2) lies in 4 d_F Z.  The quadratic form W (4 N(m) - n^2) on the lattice L4 of such classes has
     content gcd{Q(v_i), Q(v_i + v_j) - Q(v_i) - Q(v_j)} over a Z-basis v_1, v_2, v_3 (exact, all of L4); printed
     content / (4 d_F).
 (B) q in F totally positive: every nonzero integral class of B(Theta, q) has chi >= 8.  Odd part: on a Z-basis of
     the odd integral classes and on sums of pairs, a = (W/4) s(q) g^2 lies in O(Lambda) (g = c_1 coefficient),
     and chi_odd = 4 Tr(a); even part: exact enumeration of the integral classes r E1^q + n Theta_1 Theta_2 with
     chi <= 14 (r = rank, n in (1/D3) Z), least chi printed (must be >= 8).
Lattices: free types, the non-split lattices of Lemma 6.7 and superlattices, random O_F^4 + sum Z gamma_i / N,
for t = 2, 5, 3, 13, 7, 17.
Usage: python -I -B p88_check.py OUT.txt     Expected last line "ALL CHECKS PASS: True".
"""
import importlib.util
import math
import os
import random
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("p88_shape", os.path.join(HERE, "p88_shape.py"))
PS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(PS)
TL, SS, Q2, NP = PS.TL, PS.SS, PS.Q2, PS.NP


def frgcd(a, b):
    a, b = Fr(a), Fr(b)
    return Fr(math.gcd(a.numerator * b.denominator, b.numerator * a.denominator), a.denominator * b.denominator)


def content_deg4(L):
    t, W = L.t, L.W
    vs = [(Q2(t, L.L4[0][j], L.L4[1][j]), L.L4[2][j]) for j in range(3)]

    def Q(v):
        m, n = v
        return W * (4 * m.norm() - n * n)

    g = Fr(0)
    for i in range(3):
        g = frgcd(g, Q(vs[i]))
        for j in range(i + 1, 3):
            s = (vs[i][0] + vs[j][0], vs[i][1] + vs[j][1])
            g = frgcd(g, Q(s) - Q(vs[i]) - Q(vs[j]))
    return g


def odd_lattice(L, q):
    """Z-basis of {c in F : sigma(c) in J2 and -c q / 2 in J6} (the odd integral classes of B(Theta, q))."""
    t = L.t
    rows = []
    for Binv, mapf in ((L.J2inv, lambda c: c.conj()), (L.J6inv, lambda c: c * q * Fr(-1, 2))):
        imgs = [mapf(Q2(t, 1)), mapf(Q2(t, 0, 1))]               # images of the Q-basis 1, sqrt t
        for i in range(2):
            rows.append([Binv[i][0] * im.u + Binv[i][1] * im.v for im in imgs])
    classes = [{r: rows[r][k] for r in range(4) if rows[r][k] != 0} for k in range(2)]
    B = SS.lattice2(classes)
    return [Q2(t, B[0][j], B[1][j]) for j in range(2)]


def even_min(L, q, bound):
    """least chi of a nonzero integral class r E1^q + n T0 with chi <= bound (None if none)."""
    t, W = L.t, L.W
    Nq = q.norm()
    D3 = 1
    for j in range(3):
        D3 = D3 * L.L4[2][j].denominator // math.gcd(D3, L.L4[2][j].denominator)
    rmax = int(math.isqrt(int(bound / (W * Nq)) + 1)) + 1
    kmax = int(D3 * math.sqrt(float(bound / W))) + 1
    best = None
    for r in range(-rmax, rmax + 1):
        if W * Nq * r * r > bound or (r * Nq * W / 4).denominator != 1:
            continue
        m = q * Fr(-r, 2)
        for k in range(-kmax, kmax + 1):
            n = Fr(k, D3)
            if r == 0 and k == 0:
                continue
            chi = W * (Nq * r * r + n * n)
            if chi > bound or not L.in_L4(m, n):
                continue
            if best is None or chi < best:
                best = chi
    return best


def odd_min(L, q, bound):
    """least chi = W Tr(q c^2) of a nonzero odd integral class (c = c01) with chi <= bound (exact enumeration)."""
    W = L.W
    b1, b2 = odd_lattice(L, q)

    def chi(c):
        return W * (q * c * c).tr()

    G11, G22 = chi(b1), chi(b2)
    G12 = (chi(b1 + b2) - G11 - G22) / 2
    det = G11 * G22 - G12 * G12
    assert G11 > 0 and det > 0
    best = None
    k2max = math.isqrt(int(bound * G11 / det) + 1) + 1
    for k2 in range(-k2max, k2max + 1):
        rem = (bound - det * k2 * k2 / G11) / G11
        if rem < 0:
            continue
        cen = float(-G12 * k2 / G11)
        rad = math.sqrt(float(rem))
        for k1 in range(math.floor(cen - rad) - 1, math.ceil(cen + rad) + 2):
            if k1 == 0 and k2 == 0:
                continue
            x = chi(b1 * k1 + b2 * k2)
            if x <= bound and (best is None or x < best):
                best = x
    return best


def qs_for(t):
    out = []
    for d in (1, 2, 4):
        for a in range(1, 9):
            for b in range(-4, 5):
                q = Q2(t, Fr(a, d), Fr(b, d))
                if q.totally_positive() and q not in out:
                    out.append(q)
    rng = random.Random(t)
    rng.shuffle(out)
    return out[:10] + [Q2(t, 1)]


G2POS = {2: Q2(2, 2, 1), 5: Q2(5, Fr(5, 2), Fr(1, 2)), 3: Q2(3, 3, 1), 13: Q2(13, 4, 1), 7: Q2(7, 3, 1),
         17: Q2(17, Fr(5, 2), Fr(1, 2))}


def lattice_list(t, rng):
    one = Q2(t, 1)
    g2s = [one, G2POS[t]]
    out = [(one, g2, None, "free(1,%s)" % g2) for g2 in g2s]
    fs = {2: (2, 4), 5: (2, 3), 3: (2, 3), 13: (2,), 7: (2, 7), 17: (2,)}[t]
    for f in fs:
        for e24 in (False, True):
            A = NP.lattice_matrix(t, f, e24)
            out.append((one, one, A, "O^4 f=%d%s" % (f, " +e24" if e24 else "")))
            gam = [Fr(rng.randrange(f), f) for _ in range(8)]
            out.append((one, g2s[1], PS.superlattice(A, gam), "O^4 f=%d%s + Zgam/%d (1,%s)" % (
                f, " +e24" if e24 else "", f, g2s[1])))
    Ns = {2: (2, 4, 8), 5: (2, 3, 4, 5), 3: (2, 3, 4, 6), 13: (2, 13), 7: (2, 7, 14), 17: (2, 17)}[t]
    nr = {2: 12, 5: 6, 3: 6, 13: 5, 7: 5, 17: 5}[t]
    for N in Ns:
        for _ in range(nr):
            A = SS.random_lattice(rng, N)
            g2 = rng.choice(g2s)
            out.append((one, g2, A, "rand N=%d (1,%s)" % (N, g2)))
    return out


def main():
    out = open(sys.argv[1], "w", encoding="utf-8") if len(sys.argv) > 1 else None

    def say(msg):
        print(msg, flush=True)
        if out:
            out.write(msg + "\n")
            out.flush()

    ok = True
    rng = random.Random(20261026)
    nlat = nq = 0
    least, least_all = {}, {}
    for t in (2, 5, 3, 13, 7, 17):
        dF = TL.dfield(t)
        qs = qs_for(t)
        for g1, g2, A, label in lattice_list(t, rng):
            lat = PS.Lattice(t, PS.fe(t, g1), PS.fe(t, g2), A, label)
            L = lat.L
            cont = content_deg4(L)
            goodA = (cont / (4 * dF)).denominator == 1
            okB, mins = True, []
            for q in qs:
                # odd part
                bas = odd_lattice(L, q)
                for g0 in bas + [bas[0] + bas[1], bas[0] - bas[1]]:
                    g = g0.conj()                                   # c_1 = L_g, g = sigma(c01)
                    a = q.conj() * g * g * (L.W / 4)
                    chi_odd = TL.chi_herm(TL.Herm(t, 0, g, Q2(t, 0), 0, g0 * q * Fr(-1, 2), 0), L.W)
                    okB &= lat.preserves(a) and a.totally_positive() and chi_odd == 4 * a.tr() and a.tr() >= 2
                    okB &= L.integral(TL.Herm(t, 0, g, Q2(t, 0), 0, g0 * q * Fr(-1, 2), 0))
                # even part, and the least chi overall (a mixed class has chi >= both minima)
                mn = even_min(L, q, 14)
                mo = odd_min(L, q, 14)
                mins.append(mn)
                tot = min((x for x in (mn, mo) if x is not None), default=None)
                biq = PS.is_sq(q.norm()) is not None
                okB &= tot is None or tot >= (8 if biq else 12)
                key = (t, "biquadratic" if biq else "not biquadratic")
                if tot is not None and (least_all.get(key) is None or tot < least_all[key]):
                    least_all[key] = tot
                nq += 1
            nlat += 1
            ok &= goodA and okB
            em = min((m for m in mins if m is not None), default=None)
            least[t] = em if least.get(t) is None or (em is not None and em < least[t]) else least[t]
            say("t=%2d %-34s W=%-8s (A) content/(4 d_F)=%-8s %s  (B) least even chi<=14 over %d q: %s; odd a in "
                "O(Lambda): %s" % (t, label, L.W, cont / (4 * dF), goodA, len(qs), em, okB))
    say("lattices %d, (lattice, q) pairs %d, least even chi <= 14 by t: %s" % (nlat, nq, least))
    say("least chi <= 14 of a nonzero integral class of B(Theta, q), by t and kind of K: %s" % least_all)
    say("ALL CHECKS PASS: %s" % ok)


if __name__ == "__main__":
    main()
