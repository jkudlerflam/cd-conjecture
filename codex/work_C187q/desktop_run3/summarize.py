"""L = 64, 128 (C181 backfill, ref/) and L = 192 (this run, out/): two-size c over several windows, parity gap, controls."""
import json, glob, math, os, re
os.chdir(os.path.dirname(os.path.abspath(__file__)))
recs = {}
for f in glob.glob('ref/bf_*.jsonl') + glob.glob('out/*.jsonl'):
    lines = open(f).read().strip().splitlines()
    if not lines: continue
    r = json.loads(lines[-1]); c = r['cfg']
    done = (r['sweep'] + 1 == len(c['chis']))
    key = re.sub(r'^bf_', '', r['tag'])
    recs[key] = (r, done)
def fam(k): return re.sub(r'_L\d+_N\d+$', '', k)
def sbar(r, h):
    S = r['S']; m = len(S) // 2; w = S[max(0, m - 1 - h): m + h]; return sum(w) / len(w)
print('run                              L    N  done  chi  trunc     S_mid')
for k in sorted(recs):
    r, d = recs[k]; c = r['cfg']
    print(f"{k:32s} {c['L']:4d} {c['N']:4d}  {'yes' if d else 'NO ':3s} {r['chi']:4d}  {r['trunc']:.1e}  {r['S'][c['L']//2-1]:.4f}")
print('\nTwo-size c = 6 dSbar / ln(L2/L1), Sbar over 3, 7, 13 middle cuts, at N = 3L/16 (converged runs only)')
for F in sorted(set(fam(k) for k in recs)):
    pts = {}
    for k, (r, d) in recs.items():
        c = r['cfg']
        if fam(k) == F and d and c['N'] * 16 == 3 * c['L']: pts[c['L']] = r
    Ls = sorted(pts)
    for i in range(len(Ls) - 1):
        a, b = pts[Ls[i]], pts[Ls[i + 1]]
        cs = [6 * (sbar(b, h) - sbar(a, h)) / math.log(Ls[i + 1] / Ls[i]) for h in (1, 3, 6)]
        print(f'  {F:22s} {Ls[i]:3d}->{Ls[i+1]:3d}: ' + ', '.join(f'{x:.2f}' for x in cs))
print('\nParity gap E(N+1)+E(N-1)-2E(N) at N = 3L/16, with L*gap (c = 2: L*gap roughly constant)')
for F in sorted(set(fam(k) for k in recs)):
    if not F.startswith('w1n316'): continue
    for L in (64, 128, 192):
        N = 3 * L // 16; k = [f'{F}_L{L}_N{N+dn}' for dn in (0, -1, 1)]
        if all(x in recs and recs[x][1] for x in k):
            g = recs[k[1]][0]['E'] + recs[k[2]][0]['E'] - 2 * recs[k[0]][0]['E']
            print(f'  {F:22s} L={L}: gap = {g:.6f}, L*gap = {L*g:.3f}')
