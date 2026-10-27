from c11_calculation import *
from itertools import combinations,product
# Independent local continuum expansion of the leading q=3 tensor.
r=4; pairs=list(combinations(range(r),2)); K=np.diag([29.5,54.25,54.25,54.25,54.25,29.5]); K[0,5]=K[5,0]=-19.5; K[1,4]=K[4,1]=-5.25; K[2,3]=K[3,2]=5.25
out={}
for i,A in enumerate(pairs):
 for j,B in enumerate(pairs):
  for lr in product(range(2),repeat=4):
   if sum(lr)%2:continue
   c=4*K[i,j]*(-1)**(lr[2]+lr[3]); ix=[fl+r*ch for fl,ch in zip(A+B,lr)]
   if len(set(ix))<4:continue
   c*=(-1)**sum(ix[a]>ix[b] for a in range(4) for b in range(a+1,4)); key=tuple(sorted(ix));out[key]=out.get(key,0)+c
M=np.array([[out.get(A+tuple(x+r for x in B),0) for B in pairs]for A in pairs])
print('RRLL matrix\n',M);print('eigenvalues',np.linalg.eigvalsh(M));print('pure chiral',out[(0,1,2,3)],out[(4,5,6,7)])
# Hodge-star basis on two-forms: 12<->34,13<->-24,14<->23.
U=np.zeros((6,6))
for p,(a,b,sgn) in enumerate([(0,5,1),(1,4,-1),(2,3,1)]):
 U[a,p]=U[a,p+3]=1/np.sqrt(2);U[b,p]=sgn/np.sqrt(2);U[b,p+3]=-sgn/np.sqrt(2)
print('self dual / anti self dual\n',U.T@M@U)
