"""Ring spectroscopy driver for C178. Momentum-resolved low levels of (3.2) or a control at fixed (L, N),
twists 0 and pi, and for the lowest states at twist 0 the valley imbalance M (eigenvalues of M in each
degenerate multiplet) and, in self-conjugate sectors at phi = pi/2, the S parity.
Usage: python3 spec.py L N phi wabs model nlev nobs out.jsonl     model in {cusp, rods, legs}"""
import sys, json, time, numpy as np
import kring as Kr
import os
FULLMAX = int(os.environ.get('FULLMAX', '42'))

def params(model, w):
    if model == 'cusp': return dict(w=w, t2=-1.0, V2=-2.0, V3=-2/3, blockade=True)
    if model == 'rods': return dict(w=w, t2=0.0, V2=0.0, V3=0.0, blockade=True)
    if model == 'legs': return dict(w=0.0, t2=-1.0, V2=0.0, V3=0.0, blockade=False)
    raise ValueError(model)

def main():
    L, N = int(sys.argv[1]), int(sys.argv[2]); phi, wa = float(sys.argv[3]), float(sys.argv[4])
    model, nlev, nobs, out = sys.argv[5], int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
    w = wa * np.exp(1j * phi); pr = params(model, w)
    t0 = time.time()
    B = Kr.basis(L, N, pr['blockade'])
    rec = dict(L=L, N=N, phi=phi, wabs=wa, model=model, nreps=len(B[0]), levels={})
    sym = (abs(phi - np.pi / 2) < 1e-12 or model == 'legs') and L % 2 == 0
    for th in (0.0, np.pi):
        lev = {}
        nl = nlev if th == 0.0 else max(3, nlev // 2)
        for m in range(L):
            if sym:  # S partner: theta = 0: m' = NL/2 - m ; theta = pi: m' = NL/2 - m - N (verified numerically)
                mp = (N * L // 2 - m - (N if th != 0.0 else 0)) % L
                if mp in lev:
                    lev[m] = lev[mp]; continue
            K = 2 * np.pi * m / L
            H, reps, P, keep = Kr.hamiltonian2(L, N, K, pr['w'], t2=pr['t2'], V2=pr['V2'], V3=pr['V3'], blockade=pr['blockade'], theta=th, B=B)
            if H.shape[0] == 0: continue
            e, _ = Kr.lowest_par(H, k=nl)
            lev[m] = [float(x) for x in e]
        rec['levels'][f"{th:.6f}"] = lev
        print(f"theta={th:.3f} done {time.time()-t0:.1f}s", flush=True)
    # observables at twist 0
    lev = rec['levels'][f"{0.0:.6f}"]
    allst = sorted([(E, m, i) for m, es in lev.items() for i, E in enumerate(es)])
    E0 = allst[0][0]
    want = {}
    for (E, m, i) in allst[:nobs]: want.setdefault(m, []).append(i)
    states = []
    for m, idxs in want.items():
        K = 2 * np.pi * m / L
        H, reps, P, keep = Kr.hamiltonian2(L, N, K, pr['w'], t2=pr['t2'], V2=pr['V2'], V3=pr['V3'], blockade=pr['blockade'], theta=0.0, B=B)
        e, v = Kr.lowest_par(H, k=max(nlev, max(idxs) + 3))
        selfconj = abs(np.exp(1j * (2 * K - N * np.pi)) - 1) < 1e-9 and model == 'cusp' and abs(phi - np.pi / 2) < 1e-12 and L % 2 == 0
        S = Kr.S_matrix_fast(L, N, K, reps, P, keep) if selfconj else None
        done = set()
        for i in sorted(idxs):
            if i in done: continue
            grp = [j for j in range(len(e)) if abs(e[j] - e[i]) < 1e-7]
            done |= set(grp)
            V = v[:, grp]
            Cm = Kr.imbalance_matrix(L, N, K, V, reps, P, keep, blockade=pr['blockade'], full=False)
            Cev = np.linalg.eigvalsh((Cm + Cm.conj().T) / 2)
            if L <= FULLMAX:
                Mm = Kr.imbalance_matrix(L, N, K, V, reps, P, keep, blockade=pr['blockade'])
                Mev = np.linalg.eigvalsh((Mm + Mm.conj().T) / 2)
            else:
                Mev = []
            Sev = None
            if S is not None:
                Ss = V.conj().T @ (S @ V)
                Sev = sorted(np.round(np.linalg.eigvals(Ss).real, 4).tolist())
            states.append(dict(E=float(e[i]), dE=float(e[i] - E0), m=int(m), deg=len(grp), M=[float(x) for x in Mev], C1=[float(x) for x in Cev], S=Sev))
    states.sort(key=lambda d: d['E'])
    rec['states'] = states; rec['time'] = time.time() - t0
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(f"L={L} N={N} phi={phi:.4f} {model}: E0={E0:.8f}  ({rec['time']:.0f}s)")
    for d in states:
        print(f"  dE={d['dE']:.5f} m={d['m']:3d} (K/pi={2*d['m']/L:.3f}) deg={d['deg']} M={np.round(d['M'],2).tolist()} C1={np.round(d['C1'],2).tolist()} S={d['S']}")

if __name__ == '__main__':
    main()
