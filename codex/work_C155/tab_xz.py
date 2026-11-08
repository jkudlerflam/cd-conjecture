"""Tabulate xz scan outputs. Usage: python3 tab_xz.py DIR"""
import json, glob, sys, numpy as np
for fn in sorted(glob.glob(sys.argv[1] + '/pt_*.json'), key=lambda s: int(s.split('_')[-1].split('.')[0])):
    r = json.load(open(fn)); e1 = r['ed_12_16']; e2 = r['ed_16_20']; d = r['dmrg']
    dm = ' '.join(f"L{L}:{x['c_prof']:.2f}/{x['c_sm']:.2f}(tr{x['trunc']:.0e},dE{x['dE_last']:.0e})" for L, x in d.get('L', {}).items())
    print(f"{r['dir']:8s} d={r['delta']:.2f} ED12/16 S {e1['c_S']:.2f} cas {e1['c_cas']:.2f} crit {e1['crit']:.2f} ndeg {e1['ndeg']} | "
          f"ED16/20 S {e2['c_S']:.2f} cas {e2['c_cas']:.2f} crit {e2['crit']:.2f} | {dm} two {d.get('c_two', np.nan):.2f} "
          f"two128/256 {d.get('c_two_128_256', np.nan):.2f} {r['labels']}")
