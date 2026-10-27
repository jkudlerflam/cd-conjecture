import json,numpy as np
from d3 import Sector,measure,append
for L in [8,10,12]:
 sectors=[Sector('spin1',L,[m],cache_blocks=L<=10) for m in range(L+1)]
 for delta in [0,.5,1,1.5]:
  for D in [3,3.5,4,4.5]:
   c=[1,1,delta-1,D];label=f'sat_z{delta}_D{D}';mins=[];mom={str(k):[] for k in range(L//2+1)}
   for m,sec in enumerate(sectors):
    here=[]
    for k in range(L//2+1):
     ev,_,_=sec.solve(c,k,1)
     if not sec.cache_blocks:sec._blocks.clear()
     if len(ev):
      here.append((float(ev[0]),k));mom[str(k)].append({'M':m,'E':float(ev[0])})
    best=min(here);mins.append({'M':m,'E0':best[0],'k':best[1]})
   e0=min(x['E0'] for x in mins);m=min(x['M'] for x in mins if abs(x['E0']-e0)<1e-8)
   result=measure(sectors[m],c,label,6);result['charge_audit']=mins
   result['global_momentum_minima']={k:min(es,key=lambda x:x['E']) for k,es in mom.items()}
   result['ground_charge_set_nonnegative']=[x['M'] for x in mins if abs(x['E0']-e0)<1e-8]
   append('saturation.jsonl',result)
