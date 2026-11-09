#!/usr/bin/env python3
"""C166 summary (Claude Opus 5.5 agent, 2026-11-08). usage: python3 summarize.py [ROOT]
ROOT contains oc/ (open_chain.py outputs) and vc/ (vumps_cells.py outputs); default '.'.
Reference open-chain S_mid at L = 32..256 from C159/C163/C155 is in inputs/ref_open.json.
Prints: (1) VUMPS table and the lowest-energy cell per bond with its FES slope; (2) open chains without field:
c_prof, S_mid, truncation, two-size slopes, S(x) at fixed x, emid against e0, order-parameter correlation;
(3) pinning; (4) the pre-fixed criterion C1-C3 of RESULT.md for each XZ bond, and the FES second stage for the
C163 path passes."""
import os, sys, json, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fesfit_c158 as F
import warnings
warnings.filterwarnings('ignore')

root = sys.argv[1] if len(sys.argv) > 1 else '.'
ref = {k: {int(L): r for L, r in v.items()} for k, v in json.load(open(os.path.join(HERE, 'inputs', 'ref_open.json'))).items()}

# ---------------------------------------------------------------- VUMPS
print('== (1) VUMPS cells (vc/): last chi reached per file')
V = {}
for fn in sorted(glob.glob(os.path.join(root, 'vc', '*.json'))):
    try:
        d = json.load(open(fn))
    except Exception as ex:
        print('  unreadable', fn, ex); continue
    R = d['results']
    if not R:
        continue
    nm = d['args']['name']; n = d['args']['n']; sd = d['args']['seed']
    r = R[-1]
    m = np.array(r['mag'])
    fit = {}
    if len(R) >= 4:
        try:
            fit = F.analyse(fn, 3, quiet=True)
        except Exception as ex:
            fit = {'err': str(ex)}
    half = [x for x in R if x['chi'] <= r['chi'] // 2]
    de = abs(r['e'] - half[-1]['e']) if half else float('nan')
    rec = dict(file=os.path.basename(fn), n=n, seed=sd, chi=r['chi'], e=r['e'], de_half=de, S=float(np.mean(r['S'])),
               Sspread=float(np.ptp(r['S'])), xi=r['xi'], m=m, mstag=float(np.max(np.abs(m))), tr=r.get('tr_radius'),
               bond_e=r.get('bond_e'), c_ls3=fit.get('c_ls3', float('nan')), c_pair=fit.get('c_pair', [])[-1:] or [float('nan')],
               xis=[x['xi'] for x in R][-4:], chis=[x['chi'] for x in R][-4:], argl=float(np.angle(complex(*r['lam'][1])) / np.pi))
    # a non-injective uMPS (odd cell for a Neel state: degenerate transfer matrix, xi -> infinity) makes the
    # energy_AL estimate unreliable; such runs are listed but excluded from the lowest-energy selection
    rec['valid'] = bool(np.isfinite(r['xi']) and r['xi'] < 1e6 and -10 < r['e'] < 10 and r['err'] < 1e-4)
    V.setdefault(nm, []).append(rec)
    tr = rec['tr'] if isinstance(rec['tr'], (int, float)) else float('nan')
    print(f"  {nm:10s} n={n:2d} s={sd} chi={r['chi']:4d} e={r['e']:.12f} (|e-e(chi/2)|={de:.1e}) S={rec['S']:.5f}"
          f" (spread {rec['Sspread']:.1e}) xi={r['xi']:8.2f} |m|max={rec['mstag']:.4f} TR={tr:.5f}"
          f" bond_e spread={np.ptp(r['bond_e']) if r.get('bond_e') else float('nan'):.1e} c_ls3={rec['c_ls3']:.3f}"
          f" last c_pair={rec['c_pair'][0]:.3f}{'' if rec['valid'] else ' [non-injective or unconverged, excluded]'} arg1/pi={rec['argl']:+.3f} xi(last chis {rec['chis']})={np.round(rec['xis'], 1).tolist()}")
best = {}
for nm, L in V.items():
    Lv = [x for x in L if x['valid']]
    if not Lv:
        continue
    b = min(Lv, key=lambda x: x['e'])
    best[nm] = b
    print(f"  lowest e for {nm}: {b['file']} e0={b['e']:.12f} chi={b['chi']} FES c_ls3={b['c_ls3']:.3f} m_site0={np.round(b['m'][:, 0], 4).tolist()}")

# ---------------------------------------------------------------- open chains
print('\n== (2) open chains without edge field (oc/)')
OC = {}
for fn in sorted(glob.glob(os.path.join(root, 'oc', '*.json'))):
    try:
        d = json.load(open(fn))
    except Exception as ex:
        print('  unreadable', fn, ex); continue
    OC.setdefault((d['name'], d['tag']), {})[d['L']] = d
for (nm, tag), R in sorted(OC.items()):
    if tag:
        continue
    print(f"  {nm}")
    e0 = best[nm]['e'] if nm in best else float('nan')
    for L, d in sorted(R.items()):
        C = np.array(d['corr'])[:, :, -1]
        nv = best[nm]['m'][:, 0] if nm in best else np.array([0, 0, 1.0])
        nh = nv / np.linalg.norm(nv)
        sgn = (-1) ** (d['corr_js'][-1] - d['corr_i'])
        op = float(nh @ C @ nh) * sgn
        m2 = float(nv @ nv) if nm in best else float('nan')
        print(f"    L={L:4d} c_prof={d['c_prof']:.3f} c_sm={d['c_sm']:.3f} Smid={d['Smid']:.4f} trunc={d['trunc']:.1e}"
              f" E={d['E']:.8f} emid={d['emid']:.10f} emid-e0={d['emid'] - e0:+.2e} S(1)={d['S'][0]:.4f}"
              f" order corr (n.C.n)(-1)^r at r=L/2: {op:+.4f} (VUMPS m^2={m2:.4f}) t={d['time_s']:.0f}s")
    S = {L: r['Smid'] for L, r in ref.get(nm, {}).items()}
    for L, d in R.items():
        S[L] = d['Smid']
    Ls = sorted(S)
    print('    two-size slopes: ' + '  '.join(f"{Ls[i]}-{Ls[i+1]}: {6 * (S[Ls[i+1]] - S[Ls[i]]) / np.log(Ls[i+1] / Ls[i]):.3f}" for i in range(len(Ls) - 1)))
    Sx = {}
    for L, d in R.items():
        Sx[L] = d['S']
    for L, r in ref.get(nm, {}).items():
        if 'S' in r and L not in Sx:
            Sx[L] = r['S']
    print('    S(x) at fixed x (smoothed over x-1..x+1):')
    for L in sorted(Sx):
        s = np.array(Sx[L])
        print(f"      L={L:4d} " + ' '.join(f"x={x}:{np.mean(s[x - 2:x + 1]):.4f}" for x in (8, 16, 32, 64, 128) if x + 1 < L // 2))

# ---------------------------------------------------------------- pinning
print('\n== (3) pinning (edge fields along the VUMPS order parameter: Lst left only, LRst both ends compatible;'
      ' Lx perpendicular, left only)')
for (nm, tag), R in sorted(OC.items()):
    if not tag:
        continue
    base = OC.get((nm, ''), {})
    for L, d in sorted(R.items()):
        b = base.get(L) or ref.get(nm, {}).get(L)
        bS = b['Smid'] if b else float('nan')
        C = np.array(d['corr'])[:, :, -1]
        m = np.array(d['mag'])
        stag = np.mean([m[:, j] * (-1) ** j for j in range(L // 2 - 4, L // 2 + 4)], axis=0)
        print(f"  {nm:10s} {tag:5s} L={L:4d} c_prof={d['c_prof']:.3f} Smid={d['Smid']:.4f} (no field {bS:.4f}, drop {bS - d['Smid']:+.4f})"
              f" trunc={d['trunc']:.1e} emid={d['emid']:.10f} central staggered m={np.round(stag, 4).tolist()}")
    Ls = sorted(R)
    if len(Ls) >= 2:
        print(f"  {nm:10s} {tag:5s} Smid({Ls[-1]}) - Smid({Ls[0]}) = {R[Ls[-1]]['Smid'] - R[Ls[0]]['Smid']:+.4f}")

# ---------------------------------------------------------------- criterion
print('\n== (4) pre-fixed criterion (RESULT.md)')
for nm in ('hi_T_40', 'dist_T_50', 'cand_60', 'cand_61'):
    R = dict(ref.get(nm, {}))
    R.update(OC.get((nm, ''), {}))
    ok1 = None
    if 256 in R and 512 in R:
        a, b = R[256], R[512]
        sl = 6 * (b['Smid'] - a['Smid']) / np.log(2)
        tr_ok = max(a['trunc'], b['trunc']) < 1e-8
        ok1 = a['c_prof'] > 1.2 and b['c_prof'] > 1.2 and b['c_prof'] >= a['c_prof'] - 0.03 and sl > 1.2 and tr_ok
        print(f"  {nm:10s} C1: c_prof 256/512 = {a['c_prof']:.3f}/{b['c_prof']:.3f}, slope 256-512 = {sl:.3f}, trunc ok = {tr_ok} -> {ok1}")
    else:
        print(f"  {nm:10s} C1: L = 512 missing")
    ok2 = None
    if nm in best:
        bb = best[nm]
        ok2 = bb['chi'] >= 256 and bb['c_ls3'] > 1.2
        print(f"  {nm:10s} C2: lowest cell {bb['file']} chi={bb['chi']} c_ls3={bb['c_ls3']:.3f} -> {ok2}")
    ok3 = None
    if nm in best and 512 in R and 256 in R and 'emid' in R[512] and 'emid' in R[256]:
        e0 = best[nm]['e']
        allow = 2 * abs(R[512]['emid'] - R[256]['emid']) + 1e-6
        ok3 = abs(R[512]['emid'] - e0) <= allow
        print(f"  {nm:10s} C3: emid(512)={R[512]['emid']:.10f} e0={e0:.10f} diff={R[512]['emid'] - e0:+.2e} allowance={allow:.1e} -> {ok3}")
    print(f"  {nm:10s} c > 1 candidate: {bool(ok1) and bool(ok2) and bool(ok3)}")
for nm in ('KI22', 'KI26', 'N34'):
    if nm in best:
        bb = best[nm]
        print(f"  {nm:10s} FES second stage: lowest cell {bb['file']} chi={bb['chi']} S={bb['S']:.4f} xi={bb['xi']:.1f}"
              f" c_ls3={bb['c_ls3']:.3f} -> passes (> 1.2 at chi >= 256): {bb['chi'] >= 256 and bb['c_ls3'] > 1.2}")
