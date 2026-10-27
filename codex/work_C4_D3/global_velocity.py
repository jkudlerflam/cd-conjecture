#!/usr/bin/env python3
"""Audit the adjacent-momentum minimum over every magnetization sector."""
import argparse,json,numpy as np
from pathlib import Path
from d3 import Sector,measure,append
p=argparse.ArgumentParser();p.add_argument('files',nargs='+');p.add_argument('--output',required=True);p.add_argument('--largest-degenerate-entropy',action='store_true');a=p.parse_args()
seen=set()
if Path(a.output).exists():
 for line in Path(a.output).read_text().splitlines():
  old=json.loads(line);seen.add((old['label'],old['L']))
for fn in a.files:
 for line in open(fn):
  r=json.loads(line);L=r['L'];c=r['coeff'];k0=r['ground_k']
  if (r['label'],L) in seen:continue
  def fold(k):k%=L;return min(k,L-k)
  ks=sorted(set(fold(k0+s) for s in [-1,1]));mins={str(k):[] for k in ks}
  for m in range(L+1):
   sec=Sector('spin1',L,[m],cache_blocks=False)
   for k in ks:
    ev,_,res=sec.solve(c,k,1);sec._blocks.clear()
    if len(ev):mins[str(k)].append({'M':m,'E':float(ev[0]),'residual':res})
  lowest={k:min(es,key=lambda x:x['E']) for k,es in mins.items()}
  gap=min(x['E'] for x in lowest.values())-r['E0']
  row={'label':r['label'],'L':L,'source':fn,'ground_k':k0,'E0':r['E0'],'adjacent_minima':lowest,'global_adjacent_gap':gap,'global_v_current_assumed':L*gap/(2*np.pi),'sector_v_current_assumed':r['v_current_assumed']}
  audit=r.get('charge_audit',r.get('other_charge_minima',[]))
  deg=sorted(set([int(r['charge'][0])]+[x['M'] for x in audit if abs(x['E0']-r['E0'])<1e-8]))
  row['ground_charge_set_nonnegative']=deg
  if a.largest_degenerate_entropy and max(deg)>int(r['charge'][0]):
   chosen=measure(Sector('spin1',L,[max(deg)],cache_blocks=False),c,r['label']+'_largest_ground_M',6)
   row['largest_ground_charge_result']=chosen
  with open(a.output,'a') as f:f.write(json.dumps(row)+'\n')
  print(r['label'],L,'ground M',deg,'v sector/global',r['v_current_assumed'],row['global_v_current_assumed'],flush=True)
