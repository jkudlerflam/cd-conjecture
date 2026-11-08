"""C159 summary of scan/, refine/, fes/ (Claude Opus agent, 2026-11-07). usage: python3 summarize_tk.py"""
import json, glob, os
import numpy as np

recs = []
for fn in glob.glob('scan/*.jsonl'):
    for l in open(fn):
        try:
            r = json.loads(l)
        except Exception:
            continue
        if 'c_S' in r:
            recs.append(r)
print('== Stage A: ring scans (12, 16)')
for mo in sorted({r['mode'] for r in recs}):
    R = [r for r in recs if r['mode'] == mo]
    gl = [r for r in R if r['crit'] < 0.15 and r['ndeg'][-1] < 3]
    co = np.array([r['c_obj'] for r in gl]) if gl else np.array([np.nan])
    cs = np.array([r['c_S'] for r in gl]) if gl else np.array([np.nan])
    ising_sym = [r for r in gl if 0.3 < r['c_S'] < 0.7 and list(r['ndeg']) == [1, 1]]
    deg = [r for r in R if r['ndeg'][-1] >= 2]
    print(f"{mo}: n={len(R)} gapless(crit<0.15)={len(gl)} ({len(gl)/max(len(R),1):.3f})  ndeg>=2 at L16: {len(deg)}")
    print(f"   c_S over gapless: pct 5/50/95/max = {np.nanpercentile(cs,5):.3f} {np.nanpercentile(cs,50):.3f} {np.nanpercentile(cs,95):.3f} {np.nanmax(cs):.3f}")
    print(f"   c_obj over gapless: pct 50/95/99/max = {np.nanpercentile(co,50):.3f} {np.nanpercentile(co,95):.3f} {np.nanpercentile(co,99):.3f} {np.nanmax(co):.3f}")
    print(f"   n(c_obj>1.1)={int(np.sum(co>1.1))} n(c_obj>1.2)={int(np.sum(co>1.2))} n(c_obj>1.3)={int(np.sum(co>1.3))}")
    print(f"   symmetric Ising-like (0.3<c_S<0.7, unique GS at 12 and 16): {len(ising_sym)}")
    # c_S histogram over gapless split by parent Delta
    for lo, hi in ((-9, -0.6), (-0.6, 0.6), (0.6, 1.0), (1.0, 99)):
        sub = [r for r in gl if lo <= r['Delta'] < hi]
        if sub:
            c = np.array([r['c_obj'] for r in sub])
            print(f"   parent Delta in [{lo},{hi}): n_gapless={len(sub)} median c_obj={np.median(c):.3f} max={c.max():.3f}")
    best = sorted(gl, key=lambda r: -r['c_obj'])[:5]
    for r in best:
        print('   top', np.round(r['vec'], 3).tolist(), f"c_S={r['c_S']:.3f} c_cas={r['c_cas']:.3f} crit={r['crit']:.3f} Delta={r['Delta']:.2f} th={r['theta']:.2f} harm={r['harm']:.2f}")

print('\n== Stage B: refinement')
for fn in sorted(glob.glob('refine/*.json')):
    d = json.load(open(fn))
    r16 = d.get('ring16_20', {})
    s = f"{os.path.basename(fn)}: scan c_S={d['scan']['c_S']:.3f} c_obj={d['scan']['c_obj']:.3f}"
    if isinstance(r16, dict):
        s += f" | ring16/20 c_S={r16.get('c_S', np.nan):.3f} c_cas={r16.get('c_cas', np.nan):.3f} crit={r16.get('crit', np.nan):.3f} ndeg={r16.get('ndeg')}"
    dm = d.get('dmrg')
    if isinstance(dm, dict):
        s += ' | DMRG ' + ' '.join(f"L{L}: prof={x['c_prof']:.3f} sm={x['c_sm']:.3f} tr={x['trunc']:.0e}" for L, x in dm['L'].items())
        s += f" c_two={dm.get('c_two', np.nan):.3f}"
    print(s)

print('\n== Stage V: VUMPS (transfer-matrix momenta, FES)')
try:
    import fesfit_c158 as F
except Exception:
    F = None
for fn in sorted(glob.glob('fes/*.json')):
    d = json.load(open(fn))
    R = d['results']
    if not R:
        continue
    print(os.path.basename(fn), 'n =', d['args']['n'])
    for r in R:
        lam = np.array([complex(*z) for z in r['lam']])
        mods = np.abs(lam) / abs(lam[0]); args = np.angle(lam) / np.pi
        S = np.mean(r['S'])
        br = ' '.join(f"{m:.4f}@{a:+.4f}" for m, a in zip(mods[1:7], args[1:7]))
        print(f"  chi={r['chi']} S={S:.4f} xi={r['xi']:.2f} |lam|@arg/pi: {br}")
    if F is not None and len(R) >= 3:
        try:
            F.analyse(fn)
        except Exception as ex:
            print('  fesfit failed', ex)
