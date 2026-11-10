"""Low-band level statistics of the full qubit ring (the blockaded band, |P| = Lucas(L) states),
resolved by K and Z-parity, pooled over sectors. Reuses full_levelstat.py's builders.
If the full chain had an exact (dressed) U(1) N, each (K, parity) low band would be a superposition of
independent N blocks and <r> would drop towards Poisson. Also reports the mean |<N>| spread per eigenstate."""
import sys, numpy as np
sys.argv = [sys.argv[0], sys.argv[1], '0'] + sys.argv[2:]
L = int(sys.argv[1])
exec(open('full_levelstat.py').read().split('# translation convention check')[0])
lucas = sum(1 for x in range(dim) if not (x & rotl(x, 1)))
Nop = np.array([bin(x).count('1') for x in range(dim)], float)
for VD in [float(x) for x in sys.argv[3].split(',')]:
    Delta = 1.0; V = VD * Delta; t2 = Delta ** 2 / (2 * V)
    Jz = V / 4; mu = -2.2 * t2 + 3 * Delta ** 2 / V; hz = (mu - 4 * Jz) / 2
    for label, D in (("candidate family D = t2/2", t2 / 2), ("control D = 0", 0.0)):
        H = build(0.5, -0.5, D, Jz, hz)
        rall = []; gaps = []; nvar = []
        for m in range(L):
            if m in (0, L // 4, L // 2, 3 * L // 4):
                continue
            for p in (0, 1):
                # sector basis states and their N (N is constant on orbits)
                K = 2 * np.pi * m / L
                A = sector(H, m, p)
                e, vec = np.linalg.eigh(A)
                R = [int(r) for r in ureps if par[r] == p and abs(np.exp(-1j * K * per_of[int(r)]) - 1) < 1e-9]
                blk = np.array([not (r & rotl(r, 1)) for r in R])
                nb = int(blk.sum())
                gaps.append(e[nb] - e[nb - 1])
                el = e[:nb]
                Nr = Nop[R]
                w2 = np.abs(vec[:, :nb]) ** 2
                Nm = Nr @ w2; N2 = (Nr ** 2) @ w2
                nvar.append(np.mean(N2 - Nm ** 2))
                sp_ = np.diff(el); r = np.minimum(sp_[:-1], sp_[1:]) / np.maximum(sp_[:-1], sp_[1:])
                rall += list(r[int(0.1 * len(r)):int(0.9 * len(r))])
        print(f"L={L} V/Delta={VD:g} {label:26s}: low band <r> = {np.mean(rall):.4f} +- {np.std(rall)/np.sqrt(len(rall)):.4f} (n={len(rall)}), "
              f"min gap to next band = {min(gaps)/t2:.1f} t2, mean Var(N) per low eigenstate = {np.mean(nvar):.3e}", flush=True)
