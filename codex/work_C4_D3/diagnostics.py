#!/usr/bin/env python3
"""Selected-state total-spin diagnostic. S2 is an expectation, not an irrep label."""
import json,numpy as np
from d3 import Sector,measure,append

def total_spin_sq(sec,psi):
 d=sec.d;s=(d-1)/2;codes=[];amps=[]
 first=int(sec.basis[0]);M=0
 for site in range(sec.L):M+=s-(first//d**site)%d
 for site in range(sec.L):
  digit=(sec.basis//d**site)%d;sel=digit>0;m=s-digit[sel]
  codes.append(sec.basis[sel]-d**site)
  amps.append(psi[sel]*np.sqrt(s*(s+1)-m*(m+1)))
 codes=np.concatenate(codes);amps=np.concatenate(amps)
 if not len(codes):return M*(M+1)
 ix=np.argsort(codes);codes=codes[ix];amps=amps[ix]
 starts=np.r_[0,np.flatnonzero(np.diff(codes))+1];vals=np.add.reduceat(amps,starts)
 return float(np.vdot(vals,vals).real+M*(M+1))

def diagnose(family,L,c,label,cartan=False):
 sec=Sector(family,L,cartan=cartan,cache_blocks=False)
 result=measure(sec,c,label,8)
 ev,vec,_=sec.solve(c,result['ground_k'],8);psi=sec.expand(result['ground_k'],vec[:,0]);sec._blocks.clear()
 result['S2_expectation']=total_spin_sq(sec,psi)
 return result
if __name__=='__main__':
 rows=[json.loads(x) for x in open('scan_su2.jsonl')]
 labels=['su2_s070','su2_s160','su2_s178','su2_s131','su2_s097','su2_s110','su2_s123','su2_s144','su2_s105','su2_s157','su2_s173','su2_s139','su2_s186']
 for r in rows:
  if r['L']==10 and r['label'] in labels:
   append('spin_diagnostics.jsonl',diagnose('su2',10,r['coeff'],r['label']))
 for t in [.06544984695,.1308996939,.16,.1963495408,5.9,6.1]:
  append('spin_diagnostics.jsonl',diagnose('sp4',10,[np.cos(t),np.sin(t)],f'sp4_diag{t}',True))
