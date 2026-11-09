"""C178 ring spectroscopy, version 2: twist-resolved and labelled.

For each flux theta in a grid and each momentum sector m (K = 2 pi m / L), the lowest nlev levels of (3.2) or a
control, each with its valley label C1 = <sum_j (b+_{j+1} b_j + h.c.)/2> = sum_k cos(k) n_k (diagonalised inside
degenerate multiplets). M = N_a - N_b is about C1 / cos(s) for valleys at -s and pi + s.
Symmetry used at phi = pi/2 (or w = 0): spec_m(theta) = spec_{NL/2 - m - N}(2 pi - theta) for 0 < theta < 2 pi, and m -> NL/2 - m at theta = 0
(S followed by the large gauge transformation exp(2 pi i sum_j j n_j / L)); C1 -> -C1. Checked in validate2.out.
Usage: python3 spec2.py L N phi wabs model nlev ntheta thetas_to_do out.jsonl
   thetas_to_do: comma list of indices j (theta = 2 pi j / ntheta), or 'all' (or 'half' = 0..ntheta/2 using S)."""
import sys, json, time, numpy as np
import kring as Kr


def params(model, w):
    if model == 'cusp': return dict(w=w, t2=-1.0, V2=-2.0, V3=-2 / 3, blockade=True)
    if model == 'rods': return dict(w=w, t2=0.0, V2=0.0, V3=0.0, blockade=True)
    if model == 'legs': return dict(w=0.0, t2=-1.0, V2=0.0, V3=0.0, blockade=False)
    raise ValueError(model)


def c1_matrix(L, N, K, V, reps, P, keep, blockade, th=0.0):
    # gauge-invariant (kinetic) version: Re[e^{i theta/L} <sum_j b+_{j+1} b_j>]
    newidx = -np.ones(len(reps), np.int64); newidx[keep] = np.arange(len(keep))
    d = V.shape[1]
    Cm = np.zeros((d, d), complex)
    Vc = [np.ascontiguousarray(V[:, a], dtype=np.complex128) for a in range(d)]
    rs = np.array([1, L - 1])
    for a in range(d):
        for b in range(d):
            g = Kr._expect_hops(reps, P, keep, newidx, L, float(K), blockade, False, Vc[b], Vc[a], rs)
            Cm[a, b] = 0.5 * (np.exp(1j * th / L) * g[1] + np.exp(-1j * th / L) * g[L - 1])
    return 0.5 * (Cm + Cm.conj().T)


def main():
    L, N = int(sys.argv[1]), int(sys.argv[2]); phi, wa = float(sys.argv[3]), float(sys.argv[4])
    model, nlev, nth, todo, out = sys.argv[5], int(sys.argv[6]), int(sys.argv[7]), sys.argv[8], sys.argv[9]
    w = wa * np.exp(1j * phi); pr = params(model, w)
    sym = (abs(phi - np.pi / 2) < 1e-12 or model == 'legs') and L % 2 == 0
    if todo == 'all': js = list(range(nth))
    elif todo == 'half': js = list(range(nth // 2 + 1)) if sym else list(range(nth))
    else: js = [int(x) for x in todo.split(',')]
    t0 = time.time()
    B = Kr.basis(L, N, pr['blockade'])
    levels = []  # (j, m, E, C1)
    for j in js:
        th = 2 * np.pi * j / nth
        done = {}
        for m in range(L):
            if sym and (j == 0 or 2 * j == nth):
                mp = (N * L // 2 - m - (N if 2 * j == nth else 0)) % L
                if mp in done:
                    for (E, c) in done[mp]: levels.append((j, m, E, -c))  # S flips cos k -> -cos k
                    continue
            K = 2 * np.pi * m / L
            H, reps, P, keep = Kr.hamiltonian2(L, N, K, pr['w'], t2=pr['t2'], V2=pr['V2'], V3=pr['V3'],
                                               blockade=pr['blockade'], theta=th, B=B)
            if H.shape[0] == 0: continue
            e, v = Kr.lowest_par(H, k=nlev)
            res = []
            i = 0
            while i < len(e):
                g = [i]
                while g[-1] + 1 < len(e) and abs(e[g[-1] + 1] - e[i]) < 1e-8: g.append(g[-1] + 1)
                Cm = c1_matrix(L, N, K, v[:, g], reps, P, keep, pr['blockade'], th)
                cev = np.linalg.eigvalsh(Cm)
                for a, ii in enumerate(g): res.append((float(e[ii]), float(cev[a])))
                i = g[-1] + 1
            done[m] = res
            for (E, c) in res: levels.append((j, m, E, c))
        print(f"theta index {j} done {time.time() - t0:.1f}s", flush=True)
    # fill the S images theta_j -> theta_{nth-j}, m -> NL/2 - m - N (mod L), C1 -> -C1
    if sym and todo == 'half':
        have = set(x[0] for x in levels)
        extra = []
        for (j, m, E, c) in levels:
            jj = (nth - j) % nth
            if jj in have: continue
            extra.append((jj, (N * L // 2 - m - N) % L, E, -c))
        levels += extra
    rec = dict(L=L, N=N, phi=phi, wabs=wa, model=model, nth=nth, nreps=len(B[0]), levels=levels, time=time.time() - t0)
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
    E0 = min(x[2] for x in levels)
    print(f"L={L} N={N} phi={phi:.4f} w={wa} {model}: E0={E0:.8f} ({rec['time']:.0f}s)")


if __name__ == '__main__':
    main()
