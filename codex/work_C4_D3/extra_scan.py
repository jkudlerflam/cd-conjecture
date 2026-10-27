#!/usr/bin/env python3
import argparse, json, time
from d3 import Sector,measure,append
p=argparse.ArgumentParser();p.add_argument('mode',choices=['orbital','u1','audit']);p.add_argument('--sizes',nargs='+',type=int,default=[8,10]);p.add_argument('--output',required=True);a=p.parse_args()
for L in a.sizes:
 if a.mode=='orbital':
  sec=Sector('orbital',L)
  vals=[-2,-1,-.75,-.5,-.25,0,.125,.25,.375,.5,.75,1,2]
  for sign in [1,-1]:
   for i,x in enumerate(vals):
    for y in vals[i:]:
     coeff=[sign,sign*y,sign*x];label=f'orb_s{sign}_x{x:g}_y{y:g}'
     append(a.output,measure(sec,coeff,label,4))
 elif a.mode=='u1':
  sectors=[Sector('spin1',L,[m]) for m in range(L+1)]
  for delta in [-1.5,-.5,0,.5,1,1.5]:
   for beta in [-1,0,.5,1,2]:
    for D in [-2,-1,0,1,2]:
     coeff=[1,beta,delta-1,D];label=f'u1_z{delta:g}_b{beta:g}_D{D:g}'
     out=measure(sectors[0],coeff,label,4)
     audit=[]
     for m in range(1,L+1):
      # Spin inversion covers negative charges; time reversal covers -k.
      best=min((sectors[m].solve(coeff,k,1)[0][0],k) for k in range(L//2+1) if len(sectors[m].blocks(k)[0]))
      audit.append({'M':m,'E0':float(best[0]),'k':best[1]})
     out['other_charge_minima']=audit
     winner=min(audit,key=lambda r:r['E0'])
     out['global_ground_M']=winner['M'] if winner['E0']<out['E0']-1e-8 else 0
     if out['global_ground_M']:
      out['M0_E0']=out['E0'];out['M0_entropy_fits']=out['entropy_fits']
      chosen=measure(sectors[winner['M']],coeff,label,4)
      out.update(chosen)
     append(a.output,out)
 else:
  # Audit the 24 U(1)-only deformations in the original spin-1 scan.
  from d3 import cases
  sectors=[Sector('spin1',L,[m]) for m in range(L+1)]
  for label,coeff in cases('spin1',96)[96:]:
   mins=[]
   for m,sec in enumerate(sectors):
    best=min((sec.solve(coeff,k,1)[0][0],k) for k in range(L//2+1) if len(sec.blocks(k)[0]))
    mins.append({'M':m,'E0':float(best[0]),'k':best[1]})
   row={'label':label,'L':L,'charge_minima':mins}
   with open(a.output,'a') as f:f.write(json.dumps(row)+'\n')
   print(label,L,min(mins,key=lambda r:r['E0']),flush=True)
