"""C198 readout: c_two(L -> 2L) = 6[Sbar_W(2L) - Sbar_W(L)]/ln 2 at W = 1, 6, 12 (C176's Sbar),
each size at its last chi and at matched chi, plus the Sbar_6 step of each chi stage.
Run from codex/work_C198/: python3 -E -s c_two198.py typhon_out/out > c_two198.txt"""
import json, sys, numpy as np
D = sys.argv[1] if len(sys.argv) > 1 else 'typhon_out/out'

def Sbar(S, L, W):
    S = np.asarray(S); l = np.arange(1, L); sel = (l >= L // 2 - W) & (l <= L // 2 + W)
    return float(S[sel].mean())

def bychi(tag):
    d = {}
    try:
        for line in open(f'{D}/{tag}.jsonl'):
            r = json.loads(line); d[r['chi']] = r
    except FileNotFoundError:
        pass
    return d

def ctwo(a, La, ca, b, Lb, cb):
    return [6 * (Sbar(b[cb]['S'], Lb, W) - Sbar(a[ca]['S'], La, W)) / np.log(2) for W in (1, 6, 12)]

for fam, sizes in (('sym316', (128, 256, 512, 1024)), ('c2', (512, 1024)), ('c1', (512,))):
    runs = {L: bychi(f'{fam}_L{L}') for L in sizes}
    for L, r in runs.items():
        if not r: continue
        ks = sorted(r)
        steps = ' '.join(f'{ks[i-1]}->{ks[i]}:{Sbar(r[ks[i]]["S"], L, 6) - Sbar(r[ks[i-1]]["S"], L, 6):+.4f}' for i in range(max(1, len(ks) - 4), len(ks)))
        print(f'{fam} L={L}: last chi {ks[-1]}, trunc {r[ks[-1]]["trunc"]:.1e}, Sbar6 steps {steps}')
    for La, Lb in zip(sizes, sizes[1:]):
        a, b = runs[La], runs[Lb]
        if not a or not b: continue
        # the smaller size at its last chi, the larger at each of its last three chi
        for cb in sorted(b)[-3:]:
            v = ctwo(a, La, max(a), b, Lb, cb)
            print(f'  c_two({La}->{Lb}) {La}@{max(a)} {Lb}@{cb}: W=1,6,12 ' + ' '.join(f'{x:.3f}' for x in v))
