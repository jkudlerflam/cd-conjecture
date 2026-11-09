"""Summarize the quick test: two-size c from mean mid-chain entropy (several windows), and the parity gap."""
import json, glob, math, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
last = {}
for f in glob.glob('out/*.jsonl'):
    lines = open(f).read().strip().splitlines()
    if lines:
        r = json.loads(lines[-1]); last[r['tag']] = r
def sbar(r, half):
    S = r['S']; m = len(S) // 2
    w = S[max(0, m - 1 - half): m + half]
    return sum(w) / len(w)
print('tag            L    N  sweeps  chi   trunc     E              S_mid')
for t in sorted(last, key=lambda t: (t.split('_L')[0], int(t.split('_L')[1]))):
    r = last[t]; c = r['cfg']
    print(f"{t:13s} {c['L']:4d} {c['N']:4d}  {r['sweep']+1:3d}/{len(c['chis']):<3d} {r['chi']:4d}  {r['trunc']:.1e}  {r['E']:.10f}  {r['S'][c['L']//2-1]:.4f}")
print('\nTwo-size c = 6 [Sbar(L2) - Sbar(L1)] / ln(L2/L1), Sbar = mean S over the middle 2*half+1 cuts')
for fam in ('fit8', 'sw', 'c1', 'c2'):
    Ls = sorted(int(t.split('_L')[1]) for t in last if t.split('_L')[0] == fam)
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            a, b = last[f'{fam}_L{Ls[i]}'], last[f'{fam}_L{Ls[j]}']
            cs = [6 * (sbar(b, h) - sbar(a, h)) / math.log(Ls[j] / Ls[i]) for h in (1, 3, 6)]
            print(f"  {fam:5s} {Ls[i]:3d}->{Ls[j]:3d}: " + ', '.join(f'{x:.2f}' for x in cs) + '   (windows 3, 7, 13 cuts)')
print('\nParity gap E(N+1) + E(N-1) - 2E(N) at N = L/12. Unbound two-component: L*gap roughly constant. Bound pairs: gap levels off near E_b ~ 0.02, so L*gap grows ~ linearly.')
for fam in ('fit8', 'sw'):
    for L in (48, 96, 144, 192):
        k = [f'{fam}_L{L}', f'{fam}m1_L{L}', f'{fam}p1_L{L}']
        if all(x in last for x in k):
            g = last[k[2]]['E'] + last[k[1]]['E'] - 2 * last[k[0]]['E']
            print(f'  {fam:5s} L={L}: gap = {g:.6f}, L*gap = {L*g:.4f}')
