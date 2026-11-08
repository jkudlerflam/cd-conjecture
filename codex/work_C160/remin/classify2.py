"""Final classification of the re-minimized C160 texture bonds (coordinator, 2026-11-08).
uniform: bulk variance < 1e-6 and bulk energy within 1e-7 of the best ring cell.
single-q: every Blackman-window spectral peak above 1e-3 of the maximum is a harmonic m q of the strongest peak
          (within 2 resolution cells), checked on the saved N = 2000 open-chain bulk.
"""
import json, glob, numpy as np
rows = []
for f in sorted(glob.glob('out/r_*.json')):
    if f.endswith('_local.json') and glob.glob(f.replace('_local', '')):
        continue
    r = json.load(open(f)); r['file'] = f; rows.append(r)
cnt = {}
for r in rows:
    n = np.load(r['file'].replace('.json', '.npy')).astype(float); N = len(n)
    x = n[N // 4: 3 * N // 4]; M = len(x)
    dE = r['e_bulk_best'] - r['e_ring']
    if r['var_bulk'] < 1e-6 and abs(dE) < 1e-7:
        cl, info = 'uniform', ''
    else:
        F = np.fft.rfft((x - x.mean(0)) * np.blackman(M)[:, None], n=8 * M, axis=0)
        P = (abs(F) ** 2).sum(1); k = np.fft.rfftfreq(8 * M) * 2 * np.pi
        loc = [i for i in range(1, len(P) - 1) if P[i] > P[i - 1] and P[i] >= P[i + 1] and P[i] > 1e-3 * P.max()]
        q = k[np.argmax(P)]; res = 2 * np.pi / M
        bad = [round(k[i], 4) for i in loc if min(abs(k[i] - abs(((m * q + np.pi) % (2 * np.pi)) - np.pi)) for m in range(1, 31)) > 2 * res]
        cl = 'single-q' if not bad else 'CANDIDATE'
        info = f'q={q:.4f} npeaks={len(loc)} nonharmonic={bad}'
    cnt[cl] = cnt.get(cl, 0) + 1
    print(f"{r['seed']:3d} {r['t']}  {cl:9s} var={r['var_bulk']:.1e} dE(bulk-ring)={dE:+.1e} {info}")
print(len(rows), 'bonds', cnt)
