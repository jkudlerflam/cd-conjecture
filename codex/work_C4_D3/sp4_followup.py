#!/usr/bin/env python3
import sys,json,numpy as np
from d3 import Sector,measure,append
mode=sys.argv[1]
if mode=='fine':
 sec=Sector('sp4',8,cartan=True)
 for i in range(768):
  t=2*np.pi*i/768
  append('sp4_fine.jsonl',measure(sec,[np.cos(t),np.sin(t)],f'sp4_fine_{i:04d}',2,False))
else:
 from scipy.optimize import brentq
 for L,lo,hi in [(6,.15,.27),(8,.16,.20),(10,.13,.16),(12,.10,.131)]:
  sec=Sector('sp4',L,cartan=True)
  def diff(t):
   c=[np.cos(t),np.sin(t)];e0=sec.solve(c,0,2)[0][0];epi=sec.solve(c,L//2,2)[0][0]
   return float(epi-e0)
  f_lo=diff(lo);f_hi=diff(hi)
  root=brentq(diff,lo,hi,xtol=1e-9)
  row={'L':L,'theta_crossing':root,'theta_bracket':[lo,hi],'endpoint_differences':[f_lo,f_hi],'difference_at_root':diff(root),'definition':'finite-L E(k=pi)=E(k=0), not a thermodynamic critical point'}
  with open('sp4_crossings.jsonl','a') as f:f.write(json.dumps(row)+'\n')
  print(row,flush=True)
