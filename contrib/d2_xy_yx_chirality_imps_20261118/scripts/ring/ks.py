#!/usr/bin/env python3
"""Koo-Saleur analysis on rings for H = sum[a XY + b YX + ZZ] - h sum Z (2026-11-24, Claude Opus 5.5).
Uses ring_ed.py (momentum/parity ED, written by Codex gpt-6.2-astra, self-checked against dense ED).

Energy density h_j = a X_jY_{j+1} + b Y_jX_{j+1} + Z_jZ_{j+1} - (h/2)(Z_j + Z_{j+1}), centred at j + 1/2.
Ht_n = (L/2pi) sum_j e^{i n 2pi (j+1/2)/L} h_j  ~  v (L_n + Lbar_{-n})   (one CFT; Milsted-Vidal 2017).
With T|m> = e^{2 pi i m/L}|m>,  Ht_n maps momentum sector m to m - n.
For each state s in sector m0 - n of the ground parity:
   x_s = (E_s - E0) L / 2pi  (scaled energy, = v * Delta_s),   M_s = <s|Ht_n|0>,
   for n = -2, +2 the stress-tensor candidate has  c_s = 2 |M_s|^2 / v_s^2 with v_s = x_s / 2.
So c_s summed over the T-like states gives c, one velocity per tower.
usage: python ks.py L=18 delta=0.1875 h=1.783 nev=12 ns=1,2
"""
import sys, time
import numpy as np
import ring_ed as R

JZ = 1.0
_Binit = R.Bonds.__init__
def _binit(self, orb, p, delta, h):
    _Binit(self, orb, p, delta, h)
    if JZ != 1.0:   # rescale the ZZ part of the diagonal: diag = Jz (L - 2 #walls) - h (L - 2 N)
        L = orb.L; x = orb.rep[self.ids]
        self.diag = JZ * (L - 2 * R.pop(x ^ R.translate(x, L))) - h * (L - 2 * orb.down[self.ids])
R.Bonds.__init__ = _binit


def expand(orb, ids, m, v):
    L = orb.L
    ids_x, s = np.divmod(orb.address, L)
    coeff = np.zeros(orb.rep.size, dtype=complex)
    coeff[ids] = v / np.sqrt(orb.period[ids])
    return coeff[ids_x] * np.exp((-2j * np.pi / L) * ((m * s) % L))


def apply_Hn(psi, L, delta, h, n):
    a, b = (1 + delta) / 2, (1 - delta) / 2
    x = np.arange(psi.size, dtype=np.int64)
    out = np.zeros_like(psi)
    weight = np.array([1., -delta, delta, -1.])     # (H[flip, x]) / i, index 2 b_j + b_{j+1} (as in ring_ed)
    for j in range(L):
        k = (j + 1) % L
        bj = (x >> j) & 1; bk = (x >> k) & 1
        zj = 1 - 2 * bj; zk = 1 - 2 * bk
        ph = np.exp(1j * n * 2 * np.pi * (j + 0.5) / L)
        hd = JZ * zj * zk - (h / 2) * (zj + zk)
        term = hd * psi
        flip = x ^ ((1 << j) | (1 << k))
        term_off = np.zeros_like(psi)
        term_off[flip] = 1j * weight[2 * bj + bk] * psi     # amplitude from x to flip(x)
        out += ph * (term + term_off)
    return out * (L / (2 * np.pi))


def main(kw):
    global JZ
    JZ = float(kw.get("Jz", 1.0))
    L = int(kw.get("L", 12)); delta = float(kw.get("delta", 0.1875)); h = float(kw.get("h", 1.783))
    nev = int(kw.get("nev", 12)); ns = [int(x) for x in kw.get("ns", "1,2").split(",")]
    t0 = time.time()
    orb = R.Orbits(L)
    spec = {}
    for p in (0, 1):
        bonds = R.Bonds(orb, p, delta, h)
        for m in range(L):
            H, ids = bonds.block(m)
            spec[(m, p)] = R.lowest(H, nev)
    E0, (m0, p0) = min((spec[k][0], k) for k in spec)
    print("# L=%d delta=%g h=%g Jz=%g  E0=%.12f  E0/L=%.12f  ground (m,parity)=(%d,%d)  t=%.0fs" % (L, delta, h, JZ, E0, E0 / L, m0, p0, time.time() - t0), flush=True)
    for p in (0, 1):
        for m in range(L):
            dm = (m - m0 + L // 2) % L - L // 2
            xs = (spec[(m, p)] - E0) * L / (2 * np.pi)
            print("SPEC p=%d dm=%+d  x=(E-E0)L/2pi: %s" % (p, dm, " ".join("%.4f" % x for x in xs[:8])))
    bonds = R.Bonds(orb, p0, delta, h)
    H, ids = bonds.block(m0)
    e, V = R.lowest(H, 1, vector=True)
    psi0 = expand(orb, ids, m0, V[:, 0])
    for n in sorted(set(ns + [-x for x in ns])):
        phi = apply_Hn(psi0, L, delta, h, n)
        mt = (m0 - n) % L
        Hs, idt = bonds.block(mt)
        es, Vs = R.lowest(Hs, nev, vector=True)
        tot = np.vdot(phi, phi).real
        print("KS n=%+d  ||Ht_n|0>||^2 = %.6f  -> sector dm=%+d" % (n, tot, -n))
        acc = 0
        for i in range(len(es)):
            ps = expand(orb, idt, mt, Vs[:, i])
            M = np.vdot(ps, phi); acc += abs(M) ** 2
            x = (es[i] - E0) * L / (2 * np.pi)
            cs = 2 * abs(M) ** 2 / (x / 2) ** 2 if (abs(n) == 2 and x > 0) else np.nan
            vs = x / abs(n) if x > 0 else np.nan
            print("   s=%d x=%.4f |M|^2=%.6f  v_s=x/|n|=%.4f  c_s(n=2 only)=%.4f" % (i, x, abs(M) ** 2, vs, cs))
        print("   captured fraction of ||Ht_n|0>||^2 by these %d states: %.4f" % (len(es), acc / tot), flush=True)
    print("# done t=%.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main(dict(a.split("=", 1) for a in sys.argv[1:]))
