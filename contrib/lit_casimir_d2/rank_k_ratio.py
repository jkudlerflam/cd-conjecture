from math import comb
from fractions import Fraction as F
def H(k): return sum(F(1,l) for l in range(1,k+1))
def ratio(N,k): return F(k*(N*N-1), (N+k)*N)/H(k)
print("ratio R(N,k)=sum c v/(2 pi ||h||) = k(N^2-1)/((N+k) N H_k)")
print("N\\k "+" ".join(f"{k:>7d}" for k in range(1,9)))
for N in range(2,13):
    print(f"{N:>3d} "+" ".join(f"{float(ratio(N,k)):7.4f}" for k in range(1,9)))
viol=[(comb(N+k-1,k),N,k,float(ratio(N,k))) for N in range(2,40) for k in range(1,40) if ratio(N,k)>1]
viol.sort()
print("smallest d with R>1:",viol[:6])
print("N->inf limit k/H_k:",[ (k,float(k/H(k))) for k in (1,2,3,4,5,10)])
for N in (2,3,4):
    print(N,"sup_k R:",max((float(ratio(N,k)),k) for k in range(1,400)))
for N in range(5,8):
    print(N,"first k with R>1:",min([k for k in range(1,400) if ratio(N,k)>1] or [None]))
