"""Checks of dmrg.py: (1) vs ED on an open L = 10 chain (XZ + 0.3 ZX + 0.4 YY); (2) vs exact free fermions
on open chains (XZ + 0.5 ZX, critical, c = 1; XZ + 0.3 ZX + 0.2(XX+ZZ) + 0.1 By); (3) c from S(l)."""
import numpy as np, time
from xzlib import bond, term, PAULI, I2
from dmrg import dmrg
import ff_open

def ed_open(h, L):
    N = 1 << L
    H = np.zeros((N, N), complex)
    for j in range(L - 1):
        H += np.kron(np.kron(np.eye(1 << j), h), np.eye(1 << (L - j - 2)))
    return np.linalg.eigvalsh(H)[0]

def cfit(S, L, lo=0.2, hi=0.8):
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l > lo * L) & (l < hi * L)
    out = []
    for par in (0, 1):
        mm = m & (l % 2 == par)
        out.append(6 * np.polyfit(x[mm], S[mm], 1)[0])
    return out

h = term("xz") + 0.3 * term("zx") + 0.4 * term("yy")
r = dmrg(h, 10, chis=(16, 32, 32))
print("ED open L=10:", ed_open(h, 10), " DMRG:", r["E"])
for name, cf in [("XZ+0.5ZX", dict(xz=1, zx=0.5)), ("XZ+0.3ZX+0.2(XX+ZZ)+0.1By", dict(xz=1, zx=0.3, xx=0.2, zz=0.2, By=0.1))]:
    h = sum(c * term(k if k != "By" else "By") for k, c in cf.items())
    for L in (40, 80):
        t0 = time.time()
        r = dmrg(h, L, chis=(16, 32, 64, 96, 96))
        E, S = ff_open.ground(L, cf)
        print(f"{name} L={L}: E_ff={E:.10f} E_dmrg={r['E']:.10f} dE={r['E']-E:.2e} max|dS|={np.abs(S-r['S']).max():.2e} "
              f"c_fit(ff)={np.round(cfit(S,L),3).tolist()} c_fit(dmrg)={np.round(cfit(r['S'],L),3).tolist()} ({time.time()-t0:.0f}s)", flush=True)
