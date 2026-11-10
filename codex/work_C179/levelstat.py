"""Level-spacing ratio <r> of the effective model (3.2) in resolved (N, K) sectors. GOE 0.531, GUE 0.600, Poisson 0.386.
Sectors with 2K = N pi (mod 2 pi) are skipped (there the k -> pi - k symmetry U acts inside the sector at phi = pi/2)."""
import numpy as np, sys, time
from ringk import reps, build_k, rstat
L, N = int(sys.argv[1]), int(sys.argv[2])
ms = [int(x) for x in sys.argv[3].split(',')]
data = reps(L, N)
print(f"L={L} N={N}: {len(data[0])} configurations, {len(data[3])} orbits", flush=True)
for label, w in (("phi=pi/2 |w|=1", 1j), ("phi=1.3 |w|=1", np.exp(1.3j)), ("phi=pi/2 |w|=2", 2j), ("w=0 (N_A-N_B conserved, not resolved)", 0.0)):
    rs = []
    for m in ms:
        K = 2 * np.pi * m / L
        if abs(np.exp(2j * K) - np.exp(1j * N * np.pi)) < 1e-9:
            continue
        H, R, P = build_k(L, N, K, w, data=data)
        e = np.linalg.eigvalsh(H.toarray())
        r, nr, ndeg = rstat(e)
        rs.append(r)
        print(f"  {label:40s} m={m:2d} dim={len(e):6d} <r>={r:.4f} (n={nr}, exact degeneracies {ndeg})", flush=True)
    print(f"  {label:40s} mean over sectors <r> = {np.mean(rs):.4f}", flush=True)
