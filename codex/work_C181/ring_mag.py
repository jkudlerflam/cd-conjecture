"""C181: valley 'magnetization curve' on rings. For each phi, the lowest level in every momentum sector of model (3.2)
on an L-site ring with N particles (C178's kring.py, copied unmodified), the ground state, and its valley imbalance
M = sum_{cos k>0} n(k) - sum_{cos k<0} n(k) from the JW-fermion momentum distribution. Periodic bosons (theta = 0).
Usage: python3 ring_mag.py wa L N phi1,phi2,... [V2 V3]"""
import numpy as np, sys, time, json
import kring as Kr
wa = float(sys.argv[1]); L = int(sys.argv[2]); N = int(sys.argv[3]); phis = [float(x) for x in sys.argv[4].split(',')]
V2 = float(sys.argv[5]) if len(sys.argv) > 5 else -2.0
V3 = float(sys.argv[6]) if len(sys.argv) > 6 else -2.0 / 3.0
B = Kr.basis(L, N, True)
for phi in phis:
    t0 = time.time()
    w = wa * np.exp(1j * phi)
    rows = []
    best = None
    for m in range(L):
        K = 2 * np.pi * m / L
        H, reps, P, keep = Kr.hamiltonian2(L, N, K, w, V2=V2, V3=V3, B=B)
        if H.shape[0] == 0:
            continue
        e, v = Kr.lowest_par(H, k=1)
        ks, nk = Kr.momentum_distribution(L, N, K, v[:, 0], reps, P, keep)
        M = Kr.valley_imbalance(ks, nk)
        rows.append((m, float(e[0]), float(M)))
        if best is None or e[0] < best[1]:
            best = (m, float(e[0]), float(M))
    rows.sort(key=lambda r: r[1])
    rec = dict(wa=wa, L=L, N=N, phi=phi, V2=V2, V3=V3, gs=best, low=rows[:8], secs=time.time() - t0)
    print(json.dumps(rec), flush=True)
