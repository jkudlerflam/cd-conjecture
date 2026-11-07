"""Summarize completed search files, retaining failure and accuracy diagnostics."""
import argparse,json
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument('folder',nargs='?',default='results'); args=ap.parse_args()
records=[]
for p in Path(args.folder).glob('*.json'):
 r=json.loads(p.read_text())
 if isinstance(r,dict) and 'Q' in r and 'seed' in r: records.append(r)
for L in (2,4,6,8,3,5,7):
 rr=[r for r in records if r['L']==L]
 if not rr:
  print(f'L={L}: no optimization results'); continue
 r=max(rr,key=lambda x:x['Q'])
 print(json.dumps(dict(L=L,parity=r['parity'],completed=len(rr),converged=sum(x['success'] for x in rr),best_Q=r['Q'],best_seed=r['seed'],best_coefficients_M1=r['coefficients'],best_gauge_bracket=r['M'],max_residual=max(max(x['residual_L'],x['residual_2L']) for x in rr),max_M_bracket_gap=max(x['M']['gap'] for x in rr)),indent=2))
