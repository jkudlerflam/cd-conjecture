"""Finite checks of selector bounds and exact shortest-path coefficients."""
import math,itertools,json
from fractions import Fraction
from pathlib import Path
import numpy as np
import wrapper_tools as wt

def coefficient(k,p):
 N=k+p;start=tuple(range(p));end=tuple(range(k,N));dist=lambda s:sum(s)-sum(start)
 cfg=list(itertools.combinations(range(N),p));weights={start:Fraction(1)}
 for s in sorted(cfg,key=dist):
  if s==start:continue
  v=Fraction(0)
  for j in range(p):
   if s[j]>0 and s[j]-1 not in s:
    t=s[:j]+(s[j]-1,)+s[j+1:];v+=weights.get(t,0)
  interface=sum((j in s)!=(j+1 in s) for j in range(N-1))
  weights[s]=v if s==end else v/Fraction(interface-1)
 return str(2*weights[end])

def run():
 out=[]
 for n in (1,2,3,4):
  k=n+2;b=k+1;U=1;wv=.1
  a=lambda l:2*(1-math.cos(math.pi/l))
  db=(b+1)*a(b)-b*a(b+1);ec=a(b+1)-a(b);mc=2-db
  et=min((a(b-1)-2*a(b)+a(b+1))/2,(a(b)-2*a(b+1)+a(b+2))/3)
  ts=2*U*(1/et+b/mc+1/abs(ec))
  ch={p:wt.charge(k,p,wv) for p in (1,2)}
  ee={p:float(c['chi']@c['H']@c['chi']) for p,c in ch.items()}
  mu=ts*mc+4*U-(b+1)*ee[1]+b*ee[2];e=ts*ec+ee[2]-ee[1]
  vals=[(l,(a(l)-a(b)-(l-b)*ec)/(1+abs(l-b))) for l in range(2,501) if l not in (b,b+1)]
  assert min(x for l,x in vals)>=et*(1-1e-10) and mu>0 and e<0
  out.append(dict(n=n,k=k,b=b,t_selector=ts,mu_selector=mu,energy_density=e,
    eta_selector=et,wrong_length_coefficient=ts*et-U,illegal_selector_margin=U,
    wrong_occupation_margin=U-wv*(k+1),K_selector=4*ts+9*U,
    shortest_path_coefficients={str(p):coefficient(k,p) for p in (1,2)},
    length_scan_max=500,min_secant_ratio=min(x for l,x in vals)))
 Path('bounds_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
