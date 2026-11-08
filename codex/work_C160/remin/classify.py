"""Classify the re-minimized C160 texture bonds (coordinator, 2026-11-08)."""
import json, glob
rows = [json.load(open(f)) for f in sorted(glob.glob('out/r_*.json'))]
print(f'{len(rows)} bonds')
print('seed t  sc    class        var_bulk   resid     q_fit   e_bulk-e_uni  e_bulk-e_ring  e_bulk-e_match  old_unexpl  best')
cnt = {}
for r in rows:
    eb = r['e_bulk_best']
    if r['var_bulk'] < 1e-8:
        cl = 'uniform'
    elif r['resid'] < 1e-3:
        cl = 'single-q'
    else:
        cl = 'CHECK'
    cnt[cl] = cnt.get(cl, 0) + 1
    em = r['e_ring_matched']
    print(f"{r['seed']:3d} {r['t']}  {r['sc']:.2f}  {cl:10s}  {r['var_bulk']:.2e}  {r['resid']:.2e}  {r['q_fit']:.4f}  {eb-r['e_uniform']:+.2e}   {eb-r['e_ring']:+.2e}   "
          f"{'   n/a   ' if em is None else f'{eb-em:+.2e}'}    {r['old_unexplained']:.2f}   {r['best_label']}")
print(cnt)
