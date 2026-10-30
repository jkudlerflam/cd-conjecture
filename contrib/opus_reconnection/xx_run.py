"""Large-L XX checks of the exact decomposition G = R - 2 eps_L + eps_D - H2,
of the clustering limit G - R -> -N log 2, of the candidate G <= R_L, and of
the CFT-regime inequality R_L <= (pi/4) A beta/L with A = c v = 1.
Output: xx_results.json."""
import json, time
from pathlib import Path
import numpy as np
from xx_c123 import xx_parts, xx_R, xx_histories

ROOT = Path(__file__).resolve().parent

def case(n, a, N):
    beta = 2 * a * n * N          # tau M / L = a with M = 1/2
    q = xx_histories(n, beta, N, [0])
    qd = xx_histories(2 * n, beta, N, [0, n])
    x = np.pi / (2 * n)
    defect = 2 * np.sin(x / 2) ** 2 / (np.sin(x) * np.cos(x))
    lzL, lzD = xx_parts(n, beta), xx_parts(2 * n, beta)
    R = xx_R(n, beta)
    G = beta * defect + np.log(q @ q) - np.log(qd.sum())
    epsL = lzL - np.log(q.sum()); epsD = lzD - np.log(qd.sum())
    pi = q / q.sum(); H2 = -np.log(pi @ pi)
    resid = G - (R - 2 * epsL + epsD - H2)
    return dict(L=n, a=a, N=N, beta=beta, R=R, G=float(G), epsL=float(epsL),
                epsD=float(epsD), H2=float(H2), GmR=float(G - R),
                GmR_plus_Nlog2=float(G - R + N * np.log(2)), identity_residual=float(resid),
                G_le_R=bool(G <= R + 1e-10))

def main():
    t = time.time(); rows = []
    for n in [8, 32, 128, 512, 2048, 8192]:
        for a in [0.25, 1.0, 2.0]:
            for N in [1, 2, 4]:
                r = case(n, a, N); rows.append(r)
                assert abs(r['identity_residual']) < 1e-8
        print('L', n, 'done', time.time() - t, flush=True)
    # CFT-regime inequality R_L <= (pi/4) beta/L for XX (A = 1)
    scan = []
    for n in [16, 64, 256, 1024, 4096]:
        for y in np.concatenate([np.linspace(0.05, 2, 40), np.linspace(2.5, 40, 16)]):
            beta = y * n; R = xx_R(n, beta)
            scan.append(dict(L=n, y=float(y), R=R, bound=float(np.pi * y / 4),
                             slack=float(np.pi * y / 4 - R)))
    out = dict(cases=rows, R_scan=scan, seconds=time.time() - t)
    (ROOT / 'xx_results.json').write_text(json.dumps(out, indent=1) + '\n')
    for r in rows:
        if r['L'] in (32, 8192):
            print(f"L={r['L']} a={r['a']} N={r['N']} G-R+Nlog2={r['GmR_plus_Nlog2']:.3e} "
                  f"epsL/N={r['epsL']/r['N']:.5f} epsD/N={r['epsD']/r['N']:.5f} H2/N={r['H2']/r['N']:.5f} G<=R {r['G_le_R']}")
    print('all G<=R:', all(r['G_le_R'] for r in rows))
    print('min slack of R <= pi y/4:', min(s['slack'] for s in scan),
          'at', min(scan, key=lambda s: s['slack']))

if __name__ == '__main__':
    main()
