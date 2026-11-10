"""C197 readout: c_eff(48 -> 96) = 6 [S(96) - S(48)] / ln 2 for each (D, hz), with the central-three and
fractional-window entropies of the collaborator's wrapper, on the even sector and on the ground branch (lower
energy of the two parity sectors at each size). Reference: their chi = 192/256 run at D = 1/10, hz = -7/8 gives
1.4185 / 1.4363. Frustration-free point of this family: D = 0, hz = -3/4 (degenerate weighted-Dicke manifold)."""
import json, glob, math, os
rec = {}
for f in glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', '*.json')):
    d = json.load(open(f))
    if d.get('status') != 'done' or d['L'] not in (48, 96): continue
    c = d['couplings']; rec[(c['Dz'], c['hz'], d['L'], d['parity'])] = d
def ce(a, b, key): return 6 * (b[key] - a[key]) / math.log(2)
print('   D      hz    even: c3   cfrac  | ground (P48,P96): c3   cfrac  | E_odd - E_even at 48, 96 | n_bulk(96)')
for D, h in sorted({(k[0], k[1]) for k in rec}, key=lambda t: (t[0], -t[1])):
    e48, e96 = rec.get((D, h, 48, 0)), rec.get((D, h, 96, 0))
    o48, o96 = rec.get((D, h, 48, 1)), rec.get((D, h, 96, 1))
    s = f'{D:5.2f} {h:8.4f}  '
    s += f'{ce(e48, e96, "S_central3"):7.3f} {ce(e48, e96, "S_fractional"):7.3f}  |' if e48 and e96 else '   incomplete     |'
    if e48 and e96 and o48 and o96:
        g48 = e48 if e48['E_physical'] <= o48['E_physical'] else o48
        g96 = e96 if e96['E_physical'] <= o96['E_physical'] else o96
        s += f'  ({g48["parity"]},{g96["parity"]}) {ce(g48, g96, "S_central3"):7.3f} {ce(g48, g96, "S_fractional"):7.3f}  |'
        s += f'  {o48["E_physical"] - e48["E_physical"]:+.1e} {o96["E_physical"] - e96["E_physical"]:+.1e}  | {g96["density_bulk"]:.4f}'
    else: s += '  odd sector incomplete'
    print(s)
