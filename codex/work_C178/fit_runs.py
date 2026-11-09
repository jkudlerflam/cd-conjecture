"""Fit the twist-relaxed pseudospin levels (analyze2.summarize) of each cusp run to an N-site XXZ ring (xxz_fit).
Prints the best fits with Phi free, and the best fit at the predicted residual twist Phi_s = (k_b - k_a) L (mod 2 pi,
both gauge images Phi_s and Phi_s + pi N). Usage: python3 fit_runs.py files.jsonl"""
import sys, json, numpy as np
import analyze2 as A, xxz_fit as X

def run(r):
    o = A.summarize(r, verbose=False)
    L, N = o['L'], o['N']
    if r['model'] != 'cusp' or N % 2 or N > 10: return
    data = {'D_ps': o['D_ps'], 'D_4': o.get('D_4')}
    for k, v in o['M0_by_kappa'].items():
        if k > 0: data[f'k{k}'] = v[0]
    data = {k: v for k, v in data.items() if v is not None}
    ks, es, cs = A.valley_cos(r['phi'], r['wabs'])
    dk = abs(ks[0] - ks[1])
    Phip = (dk * L) % (2 * np.pi)
    res = X.fit(N, data)
    best = res[0]
    print(f"{r['model']} L={L} N={N} n={N/L:.3f} phi={r['phi']:.4f} w={r['wabs']}: data " + " ".join(f"{k}={v:.4f}" for k, v in data.items()))
    print(f"   free Phi: rms={best[0]:.3f} Delta={best[1]:.3f} Phi/pi={best[2]/np.pi:.3f} J={best[3]:.4f}")
    # range of Delta with rms within 1.5x best
    ok = [x for x in res if x[0] < 1.5 * best[0] + 0.02]
    print(f"   Delta range (rms < 1.5 best + 0.02): [{min(x[1] for x in ok):.2f}, {max(x[1] for x in ok):.2f}]")
    for P in sorted({Phip, (Phip + np.pi * N) % (2 * np.pi), (-Phip) % (2*np.pi)}):
        rr = X.fit(N, data, Phis=np.array([P]))
        if rr:
            print(f"   Phi fixed {P/np.pi:.3f} pi: rms={rr[0][0]:.3f} Delta={rr[0][1]:.3f} J={rr[0][3]:.4f}")

if __name__ == '__main__':
    for fn in sys.argv[1:]:
        for line in open(fn):
            run(json.loads(line))
