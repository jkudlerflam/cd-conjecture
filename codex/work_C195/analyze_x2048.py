"""C195x readout: the chi = 2048 stage at L = 512 against L = 256 at chi = 1536 (codex/work_C190/gpu/typhon_out/outX1536).
Usage: python3 -E -s analyze_x2048.py"""
import json, math, os, sys
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, '..', 'work_C190'))
from analyze190 import sbar, c_prof
k = 6 / math.log(2)
G = os.path.join(H, '..', 'work_C190', 'gpu', 'typhon_out')
for tag256, tag512, d1536 in (('p153c_L256_P1_rand1', 'p153c_L512_P0_rand1', 'typhon_out/out512'),
                              ('sym40_n316_L256_P1_rand1', 'sym40_n316_L512_P1_rand1', 'typhon_out_b/out512')):
    a = json.load(open(os.path.join(G, 'outX1536', tag256 + '.json')))
    x = json.load(open(os.path.join(H, 'typhon_out_x2048', tag512 + '.json')))
    p = json.load(open(os.path.join(H, d1536, tag512 + '.json')))
    h = [r for r in x['hist'] if 'Sbar6' in r]
    p15 = [r for r in p['hist'] if 'Sbar6' in r and r['chi'] == 1536][-1]
    print(f"== {tag512}: status {x['status']}, {len(h)} sweeps at chi 2048")
    print(f"   T1: res {h[-1]['res_max']:.1e} (<=1e-8), last dE {h[-1]['E'] - h[-2]['E']:+.1e} (<=2e-8), "
          f"last dS6 {h[-1]['Sbar6'] - h[-2]['Sbar6']:+.1e} (<=1e-3), trunc {h[-1]['trunc']:.1e} (<=1e-7)")
    s15, s20 = sbar(p['S'], 512, 6), sbar(x['S'], 512, 6)
    print(f"   C4 step S6 1536 -> 2048: {s20 - s15:+.4f} (bound 0.01; previous step 1024 -> 1536 reported in REPORT 9.199)")
    e15, e20 = p15['trunc'], h[-1]['trunc']
    for W in (1, 6, 12):
        sa, sb15, sb20 = sbar(a['S'], 256, W), sbar(p['S'], 512, W), sbar(x['S'], 512, W)
        ext = sb20 + (sb20 - sb15) * e20 / (e15 - e20) if e15 > e20 else float('nan')
        print(f"   W={W:2d}: c_two(256@1536 -> 512@2048) = {k * (sb20 - sa):.3f}  (at 1536: {k * (sb15 - sa):.3f}; "
              f"R8 extrapolation through 1536, 2048: {k * (ext - sa):.3f})")
    print(f"   c_prof(512) = {c_prof(x['S'], 512):.3f}; discarded weight 1536 -> 2048: {e15:.1e} -> {e20:.1e}")
