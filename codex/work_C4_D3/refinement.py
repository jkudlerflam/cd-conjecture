#!/usr/bin/env python3
import argparse,json,numpy as np
from d3 import Sector,measure,append,cases
from analyze import cf
p=argparse.ArgumentParser();p.add_argument('mode',choices=['su2','spin1','orbital','u1','su2near']);a=p.parse_args()
def read(fn):return [json.loads(x) for x in open(fn)]
if a.mode=='su2':
 rows=[r for r in read('scan_su2.jsonl') if r['L']==8]
 chosen=sorted(rows,key=cf,reverse=True)[:28]
 sec=Sector('su2',11,cache_blocks=False)
 for r in chosen:append('refine_su2.jsonl',measure(sec,r['coeff'],r['label'],8))
elif a.mode=='su2near':
 u=np.array([-1,1,-1,1.]);u/=np.linalg.norm(u)
 b1=np.array([1,0,-1,0.])/np.sqrt(2);b2=np.array([0,1,0,-1.])/np.sqrt(2)
 cs=[('su2_uls',u),('su2_tb',np.array([0,2,3,11/3.])),('su2_heis',np.array([-3.75,-2.75,-.75,2.25]))]
 for radius in [.03,.1,.3]:
  for i in range(16):
   t=2*np.pi*i/16;c=u+radius*(np.cos(t)*b1+np.sin(t)*b2);c/=np.linalg.norm(c)
   cs.append((f'su2near_r{radius}_i{i}',c))
 for L in [6,8,10]:
  sec=Sector('su2',L)
  for label,c in cs:append('scan_su2near.jsonl',measure(sec,c,label,6))
elif a.mode=='spin1':
 cs=cases('spin1',96);indices=[8,9,10,11,12,14,16,20,23,24,60,61,62,84]
 for L in [8,10,14]:
  sec=Sector('spin1',L,cache_blocks=L<14)
  for i in indices:
   label,c=cs[i];append('refine_spin1.jsonl',measure(sec,c,label,8))
elif a.mode=='orbital':
 # L=12 is commensurate with the SU4 singlet, unlike L=6,10.
 choices=[(1,0,.25),(1,.125,.5),(1,.125,.375),(1,.125,.25),(1,.25,.25),(1,.25,.5),(1,.375,.375),(1,.5,.5),(1,1,1),(-1,.125,.375),(-1,.25,.25)]
 sec=Sector('orbital',12,cache_blocks=False)
 for sign,x,y in choices:
  append('refine_orbital.jsonl',measure(sec,[sign,sign*y,sign*x],f'orb_s{sign}_x{x:g}_y{y:g}',8))
 for c,label in [([0,1,1],'orb_decoupled'),([0,1,2],'orb_decoupled_unequal')]:
  append('refine_orbital.jsonl',measure(sec,c,label,8))
else:
 rows=[r for r in read('scan_u1.jsonl') if r['L']==10]
 labels=['u1_z0_b1_D2','u1_z0.5_b1_D2','u1_z1_b1_D2','u1_z1.5_b1_D2','u1_z1.5_b2_D1','u1_z0_b2_D0','u1_z0.5_b2_D1','u1_z1_b2_D2','u1_z1.5_b2_D2','u1_z1.5_b2_D0']
 cs=[r for r in rows if r['label'] in labels]
 for L in [12,14]:
  for r in cs:
   c=r['coeff'];label=r['label'];sec=Sector('spin1',L,[0],cache_blocks=False)
   result=measure(sec,c,label,8);mins=[{'M':0,'E0':result['E0'],'k':result['ground_k']}]
   del sec
   for m in range(1,L+1):
    sec=Sector('spin1',L,[m],cache_blocks=False);energies=[]
    for k in range(L//2+1):
     ev,_,res=sec.solve(c,k,1);sec._blocks.clear()
     if len(ev):energies.append((float(ev[0]),k))
    best=min(energies);mins.append({'M':m,'E0':best[0],'k':best[1]});del sec
   winner=min(mins,key=lambda x:x['E0'])
   result['global_ground_M']=winner['M'] if winner['E0']<result['E0']-1e-8 else 0
   if result['global_ground_M']:
    result=measure(Sector('spin1',L,[winner['M']],cache_blocks=False),c,label,8)
    result['global_ground_M']=winner['M']
   result['charge_audit']=mins
   append('refine_u1.jsonl',result)
