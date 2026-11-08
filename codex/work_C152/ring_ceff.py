"""Ring ED c_eff for the parity class. c_eff(L1,L2) = 3 [S(L2) - S(L1)] / log(L2/L1), S = half-ring entropy of the
global ground state (lowest over the two Pi sectors). Usage as a module: ceff(params, Ls)."""
import sys, json, time
import numpy as np
from model import bond, ring_data

def ceff(p, Ls=(12, 16), k=2):
    h = bond(*p)
    rec = {}
    for L in Ls:
        r = ring_data(h, L, k=k)
        rec[L] = r
    S = [rec[L]['S'] for L in Ls]
    cs = [3 * (S[i + 1] - S[i]) / np.log(Ls[i + 1] / Ls[i]) for i in range(len(Ls) - 1)]
    gaps = [(rec[L]['levels'][1] - rec[L]['levels'][0]) * L for L in Ls]
    return dict(p=list(p), Ls=list(Ls), S=S, c=cs, gapL=gaps)

if __name__ == "__main__":
    # calibration: XX (c=1), TFIM critical (c=1/2), Heisenberg (c=1), Neel gapped, Ising+field paramagnet
    cal = {'XX': (1, 1, 0, 0, 0), 'TFIM': (1, 0, 0, 0, 2.0), 'Heis': (1, 1, 0, 1, 0), 'XXZ_D2': (1, 1, 0, 2, 0),
           'XY_DM_chiral': (1, 0.5, 0.6, 0, 0.3), 'XYZ_Jx=Jz': (1, 0.3, 0, 1, 0), 'Kitaev_topo': (1, 0, 0, 0, 0.5)}
    for name, p in cal.items():
        t = time.time(); r = ceff(p, (8, 12, 16)); print(name, p, [round(x, 3) for x in r['c']], [round(x, 2) for x in r['gapL']], round(time.time() - t, 1), flush=True)
