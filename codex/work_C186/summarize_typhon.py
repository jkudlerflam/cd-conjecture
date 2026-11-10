"""C186 main-job readout (coordinator, 2026-11-10). (A) grid binding B of the ring pair at N_max = 4 and 5 (positive = bound
below the grid threshold); (B) two-size c = 6[Sbar(192) - Sbar(96)]/ln 2 where both sizes exist.
Usage: python3 summarize_typhon.py   (reads typhon_out/)"""
import json, glob, os, math
def rows(f):
    t = open(f).read().strip(); return [json.loads(l) for l in t.split('\n')] if t else []
print('(A) ring binding B at L = 48 (N_max 4 / 5; sign must agree)')
for u in ('0.522', '1.035', '1.531', '2.0'):
    line = f'  u = {u}:'
    for fr in ('0.25', '0.5', '1.0', '1.5', '2.0'):
        b = []
        for N in (4, 5):
            p = f'typhon_out/out/pb_u{u}_f{fr}_N{N}.jsonl'
            b.append(rows(p)[-1]['B'] if os.path.exists(p) and rows(p) else None)
        line += f'  v/vc={fr}: ' + '/'.join('—' if x is None else f'{x:+.4f}' for x in b)
    print(line)
k = 6 / math.log(2)
def sb(S): m = len(S) // 2; return sum(S[m - 1:m + 2]) / 3
print('(B) block DMRG, u = 1.035')
for f in sorted(glob.glob('typhon_out/bm_*.jsonl')):
    R = {r['Lb']: r for r in rows(f)}
    s = '  ' + os.path.basename(f)[:-6] + ': ' + ', '.join(f"L={3*Lb} Sbar={sb(r['S']):.4f} n={r['n_mid']:.3f} trunc={r['hist'][-1][2]:.1e}" for Lb, r in sorted(R.items()))
    if 32 in R and 64 in R: s += f"  ->  c(96->192) = {k * (sb(R[64]['S']) - sb(R[32]['S'])):.3f}"
    else: s += '  (L = 192 not finished)'
    print(s)
