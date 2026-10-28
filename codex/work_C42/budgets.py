"""C42 exact finite-D/g budgets for the certified n=1,2,3 references.
Uses only standard-library rational arithmetic. Certificates are generated
by the attributed C39 certify.py snapshot. No numerical sign decisions.
"""
from fractions import Fraction as F
from math import log10
import json,sys
from pathlib import Path
from certify import gram,det

def ceil(x):return -((-x.numerator)//x.denominator)
def lg(x):return log10(x.numerator)-log10(x.denominator)
def run(n):
 c=json.loads(Path(f'certificate_n{n}.json').read_text())
 w,v,u,eta,mu=[F(c[k]) for k in ['w','v','u','eta','mu']]
 GA,GB=gram(n,w),gram(n,w*v)
 trA=sum(GA[i][i] for i in range(n));trB=sum(GB[i][i] for i in range(n))
 tmin=det(GA)/trA**(n-1)
 assert mu<tmin
 endpoint=[]
 for bracket in c['root_intervals']:
  signs=[]
  for s in bracket:
   x=F(s);co=(1-x*x)/(1+x*x);si=2*x/(1+x*x)
   H=[[-2*GA[i][j]*co+2*eta*GB[i][j]*si+mu*(w*u)**(4*i)*(i==j) for j in range(n)] for i in range(n)]
   dh=det(H);assert dh
   signs.append(dh>0)
   M=max(F(1),n*max(abs(y) for row in H for y in row))
   endpoint.append(abs(dh)/M**(n-1))
  assert signs[0]!=signs[1]
 tau=min(F(1),tmin/8,min(endpoint)/4)
 # Saved angular lower bounds imply at least half as large a chord bound.
 delta=min(F(c['delta2_lower']),F(c['delta4_lower']))/2
 assert delta>0
 r=F(1,16384);p=2*n+4;B0=128
 D=ceil(max(2/r,8*n*B0*(1+eta)/(tau*r**(p+1))))
 g2max=min(F(1,16),tau/(192*n*(1+eta)*D**p))
 exponent=1
 while F(1,4**exponent)>g2max:exponent+=1
 g=F(1,2**exponent)
 ed=2*n*B0*(1+eta)/(D*r**(p+1))
 eg=48*n*(1+eta)*g*g*D**p
 assert ed<=tau/4 and eg<=tau/4
 B=1+2*trA+2*eta*trB+mu
 speed=(tmin/2)**n*delta**(2*n-1)/(B+1)**(n-1)
 q=(w*u)**2;S=sum(q**(2*j) for j in range(n+2));vis=q**(2*n+2)/S**2
 J0=2*vis*g*g/F(D**p)
 out=dict(n=n,delta_chord=str(delta),tau=str(tau),tmin=str(tmin),D=str(D),g_dyadic_exponent=exponent,
          normalized_speed_lower=str(speed),D_log10=log10(D),g_log10=-exponent*log10(2),
          delta_log10=lg(delta),tau_log10=lg(tau),speed_log10=lg(speed),J0_log10=lg(J0),
          finite_D_error_over_tau=float(ed/tau),finite_g_error_over_tau=float(eg/tau))
 Path(f'budgets_n{n}.json').write_text(json.dumps(out,indent=2))
 print({k:v for k,v in out.items() if k not in ['delta_chord','tau','tmin','D','normalized_speed_lower']},flush=True)
if __name__=='__main__':run(int(sys.argv[1]))
