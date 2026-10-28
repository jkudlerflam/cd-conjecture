"""Print a summary table from results/n*.json.  Usage: python3 -B summarize.py [run_label]"""
import json, os, sys, glob, math
HERE = os.path.dirname(os.path.abspath(__file__))
lab = sys.argv[1] if len(sys.argv) > 1 else 'C33_kappa1e-3'
rows = []
for f in sorted(glob.glob(os.path.join(HERE, 'results', 'n*.json')), key=lambda x: int(x.split('/n')[-1][:-5])):
    R = json.load(open(f))
    r = R['runs'].get(lab)
    if r is None:
        continue
    g = lambda k: float(r[k]) if r.get(k) is not None else float('nan')
    cs = [float(b['centre_k']) for b in r['bands']]
    kfs = [float(b['kF']) for b in r['bands']]
    print('%2d p=%s f=%3d lmF=%7.2f | d2 %.2e d4 %.2e dPT %.2e dU1 %.2e dU2 %.2e dF %.2e | w4 %s | lgv %.2f..%.2f lgv/2b2 %.2f vw %.3f | loopW %.1e fieldW %.2f | uni %s steps %s | spread %.1e | kF %.3f..%.3f'
          % (R['n'], r.get('p'), r['f'], r['log10_lambda_min_F'], g('delta2'), g('delta4'), g('dPT'), g('dU1'), g('dU2'), g('dF'),
             r['witness4'], r['log10_vmin'], r['log10_vmax'], r['log10_vmin_over_2beta2'], r['vmin_whitened'],
             r['loop_whitened'], r['field_whitened'], r['unimodal'], r['inertia_steps_ok'], max(cs) - min(cs), min(kfs), max(kfs)))
