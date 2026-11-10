"""R8 readout: the chi = 1024 stage at L = 256 against the chi = 768 runs, C4 at L = 256, and c_two(128 -> 256)
on the ground branch with L = 128 at its registered final chi = 512. Uses sbar and c_prof of ../analyze190.py.
Usage: python3 -E -s analyze_x1024.py"""
import json, glob, math, os, sys
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, '..'))
from analyze190 import sbar, c_prof
G = os.path.join(H, 'typhon_out', 'outG'); X = os.path.join(H, 'typhon_out', 'outX1024')
def load(f): return json.load(open(f))
def Ephys(d): return d['hist'][-1]['E'] * 1 + 0  # same convention within a family of files; compared only like with like
def lastS(d, L): return {W: sbar(d['S'], L, W) for W in (1, 6, 12)}
k = 6 / math.log(2)
def ground(files):
    runs = [load(f) for f in files]
    return min(runs, key=lambda d: d['hist'][-1]['E'])
print('C4 at L = 256 (|Sbar6(1024) - Sbar6(768)| <= 0.01) and the stop rule of R8')
for f in sorted(glob.glob(os.path.join(X, '*.json'))):
    name = os.path.basename(f); dx = load(f); dg = load(os.path.join(G, name)); h = dx['hist']
    s6x, s6g = sbar(dx['S'], 256, 6), sbar(dg['S'], 256, 6)
    print(f"  {name:32s} Sbar6 768 -> 1024: {s6g:.4f} -> {s6x:.4f}  step {s6x - s6g:+.4f}  "
          f"last dE {h[-1]['E'] - h[-2]['E']:+.1e}  last dSbar6 {h[-1]['Sbar6'] - h[-2]['Sbar6']:+.1e}  "
          f"trunc {h[-1]['trunc']:.1e}  res {h[-1]['res_max']:.1e}  sweeps {len(h)}")
print()
for key in ('p153c', 'sym40_n316', 'xx2'):
    g128 = ground(glob.glob(os.path.join(G, f'{key}_L128_P*.json')))
    g256 = ground(glob.glob(os.path.join(X, f'{key}_L256_P*.json')))
    g768 = ground(glob.glob(os.path.join(G, f'{key}_L256_P*.json')))
    a, b, b7 = lastS(g128, 128), lastS(g256, 256), lastS(g768, 256)
    print(f"{key}: ground parity L128 {g128['parity']}, L256 {g256['parity']}; c_prof(256, chi 1024) = {c_prof(g256['S'], 256):.3f}")
    for W in (1, 6, 12):
        print(f"   W = {W:2d}: c_two(128 -> 256) = {k * (b[W] - a[W]):.3f} at chi 1024   ({k * (b7[W] - a[W]):.3f} at chi 768)")
