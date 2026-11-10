"""C182 stage 2: the registered reading (RESULT.md Q7), from the main-N runs.
  python3 stage2_read.py [typhon_out]
The job's summary2.txt was made with C176's analyze.py, which files the N + 1 runs under the same series and
size as the main runs, so its L = 144 and 192 rows are the N + 1 states. This script reads the main runs by tag.
"""
import numpy as np, json, glob, os, sys
from analyze import Sbar, c_prof

D = sys.argv[1] if len(sys.argv) > 1 else 'typhon_out'

def last(folder, tag):
    f = f'{D}/{folder}/{tag}.jsonl'
    if not os.path.exists(f): return None
    recs = [json.loads(l) for l in open(f) if l.strip()]
    if not recs: return None
    r = recs[-1]
    r['done'] = r['sweep'] == len(r['cfg']['chis']) - 1
    r['Sprev'] = [x for x in recs if x['chi'] < r['chi']][-1]['S'] if any(x['chi'] < r['chi'] for x in recs) else None
    return r

def ctwo(ra, La, rb, Lb, W):
    return 6 * (Sbar(rb['S'], Lb, W) - Sbar(ra['S'], La, W)) / np.log(Lb / La)

sers = sorted(set(os.path.basename(f).rsplit('_L', 1)[0] for f in glob.glob(f'{D}/out2/*.jsonl')))
ctl = {}
rows = []
for s in sers:
    r = {L: last('out2', f'{s}_L{L}') for L in (96, 144, 192)}
    if not all(r.values()): print(s, 'incomplete'); continue
    a48, a96 = last('out', f'{s}_L48'), last('out', f'{s}_L96')
    c_scr = ctwo(a48, 48, a96, 96, 6) if a48 and a96 else float('nan')
    c = {(a, b): [ctwo(r[a], a, r[b], b, W) for W in (1, 6, 12)] for a, b in ((96, 144), (144, 192), (96, 192))}
    gap = {}
    for L in (96, 144, 192):
        fo = 'out' if L == 96 else 'out2'
        p, q, m = last(fo, f'{s}_L{L}_Np1'), last(fo, f'{s}_L{L}_Nm1'), (a96 if L == 96 else r[L])
        gap[L] = L * (p['E'] + q['E'] - 2 * m['E']) if p and q and m and p['done'] and q['done'] else float('nan')
    m = r[192]['meas']; nm = float(np.mean(m['dens'][48:144])); W = sorted([m['Wa_f'], m['Wb_f']])
    dchi = abs(Sbar(r[192]['S'], 192, 6) - Sbar(r[192]['Sprev'], 192, 6))
    print(f"== {s}  N = {[r[L]['cfg']['N'] for L in (96, 144, 192)]}  done = {[r[L]['done'] for L in (96, 144, 192)]}")
    print("   Sbar6:", ' '.join(f"{Sbar(r[L]['S'], L, 6):.4f}" for L in (96, 144, 192)), " trunc:", ' '.join(f"{r[L]['trunc']:.1e}" for L in (96, 144, 192)),
          f" |Sbar6(512) - Sbar6(384)| at 192 = {dchi:.4f}")
    for k in ((96, 144), (144, 192), (96, 192)):
        print(f"   c_two({k[0]} -> {k[1]}) W = 1, 6, 12: " + ' '.join(f"{x:.3f}" for x in c[k]))
    print(f"   c_two(48 -> 96, W = 6, screen) = {c_scr:.3f}   c_prof(192) = {c_prof(r[192]['S'], 192):.3f}   n_mid(192) = {nm:.4f} (N/L = {r[192]['cfg']['N']/192:.4f})   Wmin/Wmax = {W[0]/max(W[1],1e-12):.3f}")
    print("   L*Delta_P at 96 (screen), 144, 192: " + ' '.join(f"{gap[L]:.3f}" for L in (96, 144, 192)))
    if s.startswith('c1') or s.startswith('c2'): ctl[s[:2]] = c[(96, 192)][1]; continue
    rows.append((s, c_scr, c[(96, 192)][1], c[(144, 192)][1], r[192]['trunc'], gap, nm, r[192]['cfg']['N'] / 192, dchi))
print()
ok_ctl = abs(ctl.get('c1', 9) - 1) <= 0.15 and abs(ctl.get('c2', 9) - 2) <= 0.25
print(f"controls c_two(96 -> 192, W = 6): c1 = {ctl.get('c1', float('nan')):.3f}, c2 = {ctl.get('c2', float('nan')):.3f}  -> {'calibrated' if ok_ctl else 'NOT calibrated'}")
print("registered reading (Q7): candidate iff c_two(96->192, W=6) > 1.2, not below c_two(48->96) by more than 0.3, trunc(192) < 1e-7, controls calibrated, L*Delta_P growth 96 -> 192 at most 50 percent")
for s, cs, c, c2, tr, gap, nm, n0, dchi in rows:
    t1 = c > 1.2; t2 = c >= cs - 0.3; t3 = tr < 1e-7; g96, g192 = gap[96], gap[192]
    t5 = (abs(g192) <= 1.5 * abs(g96)) if np.isfinite(g96) and np.isfinite(g192) else False
    print(f"{s:22s} c(96->192) = {c:.3f} [{'ok' if t1 else 'fail'}]  vs screen {cs:.3f} [{'ok' if t2 else 'fail'}]  trunc {tr:.1e} [{'ok' if t3 else 'fail'}]  "
          f"gap {g96:.3f} -> {g192:.3f} [{'ok' if t5 else 'fail'}]  => {'CANDIDATE' if (t1 and t2 and t3 and t5 and ok_ctl) else 'not a candidate'}"
          f"   (c(144->192) = {c2:.3f}, n_mid - N/L = {nm - n0:+.4f}, chi step {dchi:.4f})")
