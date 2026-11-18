"""XXC (label-order-conserving) chains: d = n1+n2. Check that the ring spectrum equals the union of
twisted XXZ sector spectra with twists exp(2 pi i (k/nA - k'/nB)), and compare ground energies.
Claude Fable 5.2, work_C201."""
import numpy as np, itertools
from potts_xxz_inclusion import xxz_sector

def xxc_H(n1, n2, L, Delta, hz):
    d = n1+n2
    typ = [1]*n1 + [-1]*n2          # +1 type A ("up"), -1 type B ("down")
    states = list(itertools.product(range(d), repeat=L))
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states),)*2)
    for s in states:
        i = idx[s]
        for j in range(L):
            a, b = s[j], s[(j+1) % L]
            za, zb = typ[a], typ[b]
            H[i, i] += 0.5*Delta*(za*zb-1) + 0.5*hz*(za+zb)
            if za != zb:
                s2 = list(s); s2[j], s2[(j+1) % L] = b, a
                H[idx[tuple(s2)], i] += 1.0
    return H

def xxz_field(L, Delta, hz, phi, nup):
    H = xxz_sector(L, Delta, phi, nup)
    return H + hz*(2*nup-L)*np.eye(H.shape[0])

def run(n1, n2, L, Delta, hz):
    ev = np.linalg.eigvalsh(xxc_H(n1, n2, L, Delta, hz))
    cand = []
    for nup in range(L+1):
        nd = L-nup
        if (n1 == 0 and nup > 0) or (n2 == 0 and nd > 0): continue
        ks = range(nup) if (n1 > 1 and nup > 0) else [0]
        kps = range(nd) if (n2 > 1 and nd > 0) else [0]
        for k in ks:
            for kp in kps:
                phi = 2*np.pi*((k/nup if nup else 0) - (kp/nd if nd else 0))
                cand.append(np.linalg.eigvalsh(xxz_field(L, Delta, hz, phi, nup)))
    allx = np.concatenate(cand)
    worst = max(np.min(np.abs(allx-e)) for e in np.unique(np.round(ev, 9)))
    worst2 = max(np.min(np.abs(ev-e)) for e in np.unique(np.round(allx, 9)))
    print(f"(n1,n2)=({n1},{n2}) L={L} Delta={Delta} hz={hz}: E_XXC={ev[0]:.10f} min twisted XXZ={allx.min():.10f} "
          f"XXC-in-XXZ mismatch={worst:.1e} XXZ-in-XXC mismatch={worst2:.1e}")

for (n1, n2) in ((1, 2), (2, 2), (1, 3)):
    for L in (3, 4, 5):
        for Delta, hz in ((0.5, 0.0), (-0.7, 0.3), (1.0, 0.0)):
            if (n1+n2)**L > 2000: continue
            run(n1, n2, L, Delta, hz)
