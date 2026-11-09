"""DMRG of one parity-class point on an open chain. Usage: python3 dmrg_pt.py L chis out.jsonl label t phi Delta V mu
(fermion units) -- or label spin:Jx,Jy,D,Jz,hz. Records E, S(l), central-window c fits (parity resolved and averaged)."""
import sys, json, time
import numpy as np
from model import bond, spin_from_fermion
from dmrg import dmrg
def cfit(S, L, lo, hi):
    S = np.asarray(S); l = np.arange(1, L); x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l >= lo * L) & (l <= hi * L)
    par = [float(6 * np.polyfit(x[m & (l % 2 == p)], S[m & (l % 2 == p)], 1)[0]) for p in (0, 1)]
    Sa = 0.5 * (S[:-1] + S[1:]); la = l[:-1] + 0.5; xa = np.log((2 * L / np.pi) * np.sin(np.pi * la / L)); ma = (la >= lo * L) & (la <= hi * L)
    return par + [float(6 * np.polyfit(xa[ma], Sa[ma], 1)[0])]
if __name__ == "__main__":
    L = int(sys.argv[1]); chis = tuple(int(c) for c in sys.argv[2].split(",")); fout = sys.argv[3]; label = sys.argv[4]
    if sys.argv[5].startswith("spin:"):
        p = tuple(float(x) for x in sys.argv[5][5:].split(","))
    else:
        t, phi, Dl, V, mu = (float(x) for x in sys.argv[5:10]); p = spin_from_fermion(t, phi, Dl, V, mu)
    h = bond(*p)
    if np.max(np.abs(h.imag)) < 1e-14: h = h.real
    t0 = time.time()
    from model import SZ
    hz = p[4]; ons = [None] * L
    ons[0] = (0.5 * hz * SZ).real if np.isrealobj(h) else 0.5 * hz * SZ
    ons[L - 1] = ons[0]                          # restore the full field on the two edge sites
    r = dmrg(h, L, chis=chis, onsite=ons)
    rec = dict(label=label, L=L, p=[float(x) for x in p], chis=list(chis), E=r["E"], trunc=r["hist"][-1][2],
               dE_last=r["hist"][-1][1] - r["hist"][-2][1], c_25_75=cfit(r["S"], L, 0.25, 0.75), c_12_88=cfit(r["S"], L, 0.125, 0.875),
               S_mid=float(r["S"][L // 2 - 1]), S=[round(float(s), 6) for s in r["S"]], sec=round(time.time() - t0, 1))
    with open(fout, "a") as f: f.write(json.dumps(rec) + "\n")
    print(label, L, rec["c_25_75"], rec["c_12_88"], rec["trunc"], rec["dE_last"], rec["sec"])
