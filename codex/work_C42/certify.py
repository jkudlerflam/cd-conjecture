"""Exact rational certificates for the limiting C39 Gram pencils, n<=3.
No floating-point signs. Fractions throughout root isolation and angles.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import factorial
import json,time
from functools import lru_cache
def outward(l,r,bits=250):
 d=2**bits
 return F(l.numerator*d//l.denominator,d), F(-((-r.numerator*d)//r.denominator),d)

def inv(A):
 n=len(A);B=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
 for j in range(n):
  piv=next(i for i in range(j,n) if B[i][j]);B[j],B[piv]=B[piv],B[j]
  d=B[j][j];B[j]=[x/d for x in B[j]]
  for i in range(n):
   if i!=j:
    z=B[i][j]
    if z:B[i]=[x-z*y for x,y in zip(B[i],B[j])]
 return [row[n:] for row in B]

def det(A):
 A=[r[:] for r in A];d=F(1);n=len(A)
 for j in range(n):
  p=next((i for i in range(j,n) if A[i][j]),None)
  if p is None:return F(0)
  if p!=j:A[p],A[j]=A[j],A[p];d=-d
  z=A[j][j];d*=z
  for i in range(j+1,n):
   q=A[i][j]/z
   for k in range(j+1,n):A[i][k]-=q*A[j][k]
 return d

def gram(n,w):
 A=[[F(4*(i==j))-w*(abs(i-j)==1) for j in range(n)] for i in range(n)]
 K=[[A[i//n][j//n]*(i%n==j%n)+A[i%n][j%n]*(i//n==j//n) for j in range(n*n)] for i in range(n*n)]
 R=inv(K);z=[r[0] for r in R]
 return [[2*sum(z[i*n+y]*R[i*n+y][x*n+j]*z[x*n+j] for x in range(n) for y in range(n)) for j in range(n)] for i in range(n)]

def trim(p):
 while len(p)>1 and not p[-1]:p.pop()
 return p

def val(p,x):
 out=F(0)
 for a in reversed(p):out=out*x+a
 return out

def rem(p,q):
 p=p[:]
 while len(p)>=len(q) and any(p):
  z=p[-1]/q[-1];k=len(p)-len(q)
  for j,a in enumerate(q):p[k+j]-=z*a
  trim(p)
 return p

def sturm(p):
 ss=[trim(p),trim([i*p[i] for i in range(1,len(p))])]
 while len(ss[-1])>1:
  r=[-x for x in rem(ss[-2],ss[-1])]
  # Positive normalization prevents integer growth while preserving signs.
  scale=abs(r[-1]);ss.append([x/scale for x in r])
 return ss

def variation(ss,x):
 a=[val(p,x) for p in ss];s=[1 if z>0 else -1 for z in a if z]
 return sum(x!=y for x,y in zip(s,s[1:]))

def isolate(p,bits=220):
 ss=sturm(p);out=[];pending=[(F(-2),F(2))];eps=F(1,2**bits)
 assert variation(ss,F(-2))-variation(ss,F(2))==len(p)-1
 while pending:
  l,r=pending.pop();count=variation(ss,l)-variation(ss,r)
  if not count:continue
  if count==1 and r-l<eps:out.append((l,r));continue
  m=(l+r)/2
  assert val(p,m)!=0,'dyadic root needs explicit handling'
  pending.extend([(l,m),(m,r)])
 return sorted(out)

def atan_small(x,K=100):
 # Alternating-series enclosure, valid |x|<1.
 assert abs(x)<1
 s=sum((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(K))
 e=abs(x)**(2*K+1)/F(2*K+1)
 return s-e,s+e

@lru_cache(None)
def pi_interval():
 a,b=atan_small(F(1,5));c,d=atan_small(F(1,239))
 return outward(16*a-4*d,16*b-4*c)

def angles(roots):
 pl,pr=pi_interval();out=[]
 for l,r in roots:
  sign=1 if l>0 else -1
  def point(x):
   y=(x-sign)/(1+sign*x);a,b=atan_small(y)
   pp=(pl/4,pr/4) if sign==1 else (-pr/4,-pl/4)
   return 2*(pp[0]+a),2*(pp[1]+b)
  out.append(outward(point(l)[0],point(r)[1]))
 return out

def circular_lower(l,r):
 pl,pr=pi_interval();bounds=[]
 for k in range(-3,4):
  a,b=(2*k*pl,2*k*pr) if k>=0 else (2*k*pr,2*k*pl)
  ll=l-b;rr=r-a
  bounds.append(F(0) if ll<=0<=rr else min(abs(ll),abs(rr)))
 return min(bounds)

def certificate(n):
 start=time.monotonic();w=F(1,4);v=F(1,1024);u=F(1,64);eta=F(1,100);mu=F(1,10**(8*n+10))
 GA=gram(n,w);GB=gram(n,w*v)
 A=[[GA[i][j]/w**(2*i+2*j) for j in range(n)] for i in range(n)]
 B=[[GB[i][j]/w**(2*i+2*j) for j in range(n)] for i in range(n)]
 E=[[u**(4*i)*F(i==j) for j in range(n)] for i in range(n)]
 # det[2 A(x^2-1)+4 eta B x+mu E(1+x^2)], rational interpolation.
 xs=list(map(F,range(2*n+1)))
 ys=[det([[2*A[i][j]*(x*x-1)+4*eta*B[i][j]*x+mu*E[i][j]*(1+x*x) for j in range(n)] for i in range(n)]) for x in xs]
 Vinv=inv([[x**j for j in range(2*n+1)] for x in xs]);p=[sum(a*b for a,b in zip(row,ys)) for row in Vinv]
 roots=isolate(p);angs=angles(roots)
 d2=min(circular_lower(angs[i][0]-angs[j][1],angs[i][1]-angs[j][0]) for i in range(2*n) for j in range(i))
 pairs=list(combinations_with_replacement(range(2*n),2));d4=None
 for z,(i,j) in enumerate(pairs):
  for k,l in pairs[:z]:
   d=circular_lower(angs[i][0]+angs[j][0]-angs[k][1]-angs[l][1],angs[i][1]+angs[j][1]-angs[k][0]-angs[l][0])
   d4=d if d4 is None else min(d4,d)
 assert d2>0 and d4>0
 # Two arcs give sorted matching flavor order at these tiny fields.
 centers=[((angs[j][0]+angs[n+j][0])/2,(angs[j][1]+angs[n+j][1])/2) for j in range(n)]
 centersep=min((circular_lower(centers[i][0]-centers[j][1],centers[i][1]-centers[j][0]) for i in range(n) for j in range(i)),default=None)
 if centersep is not None:assert centersep>0
 out=dict(n=n,w=str(w),v=str(v),u=str(u),eta=str(eta),mu=str(mu),polynomial=list(map(str,p)),root_intervals=[[str(l),str(r)] for l,r in roots],delta2_lower=str(d2),delta4_lower=str(d4),delta2_display=float(d2),delta4_display=float(d4),center_separation_display=None if centersep is None else float(centersep),centers_display=[float((l+r)/2) for l,r in centers],seconds=time.monotonic()-start)
 with open(f'certificate_n{n}.json','w') as f:json.dump(out,f,indent=2)
 print({k:z for k,z in out.items() if k not in ['polynomial','root_intervals','delta2_lower','delta4_lower']},flush=True)
if __name__=='__main__':
 import sys
 certificate(int(sys.argv[1]))
