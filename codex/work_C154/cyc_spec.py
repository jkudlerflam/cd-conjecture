"""C154: momentum-resolved ring spectra on the cyclic stratum h = a XZ + b YX + g ZY (Claude Opus agent, 2026-11-07).
Uses C155's RingTable (est_c155.py, unmodified). For each point and ring length prints, per momentum m, the lowest
levels as L (E - E0), plus E0/L and the ring entanglement fit c_fit (Calabrese-Cardy chord fit, l = 2..L/2).
Usage: python3 cyc_spec.py L1,L2,... a,b,g [a,b,g ...]"""
import sys, json, time
import numpy as np
import est_c155 as E

def bond_abg(a, b, g, B=(0, 0, 0)):
    J = np.zeros((3, 3)); J[0, 2] = a; J[1, 0] = b; J[2, 1] = g
    return E.bond_from_JB(J, np.asarray(B, float))

Ls = [int(x) for x in sys.argv[1].split(",")]
nev = 4
for spec in sys.argv[2:]:
    a, b, g = [float(x) for x in spec.split(",")]
    h = bond_abg(a, b, g)
    for L in Ls:
        t0 = time.time()
        d = E.ring_data(h, L, nev=nev)
        S = E.entropies(d['psi'], L)
        ls = np.arange(2, L // 2 + 1)
        X = np.log(L / np.pi * np.sin(np.pi * ls / L)) / 3
        cfit = np.polyfit(X, S[ls], 1)[0]
        lev = d['lev']; E0 = d['E0']
        allw = np.sort(lev.ravel())
        print(f"(a,b,g)=({a},{b},{g}) L={L}: E0/L={E0/L:.8f} m0={d['m0']} c_fit(ring)={cfit:.3f} S(L/2)={S[L//2]:.4f} "
              f"lowest L(E-E0)={np.round(L*(allw[:8]-E0),3).tolist()} ({time.time()-t0:.0f}s)", flush=True)
        for m in range(L):
            print(f"   m={m:2d}: " + " ".join(f"{L*(x-E0):7.3f}" for x in lev[m]), flush=True)
        rec = dict(a=a, b=b, g=g, L=L, E0=E0, m0=int(d['m0']), cfit=float(cfit), S=S.tolist(), lev=(L * (lev - E0)).tolist())
        with open("cyc_spec.jsonl", "a") as f:
            f.write(json.dumps(rec) + "\n")
