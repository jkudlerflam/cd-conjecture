from c11_calculation import *
from itertools import combinations,product

def cbasis(gam):
 h=len(gam[0]); out=[((),np.eye(h))]
 for k in range(1,len(gam)+1):
  for ix in combinations(range(len(gam)),k):
   M=reduce(np.matmul,[gam[a] for a in ix]); M=(1j)**(k*(k-1)//2)*M
   out.append((tuple(a+1 for a in ix),M))
 return out

for qp,qm in [(1,1),(1,2)]:
 H,G,F,gam,eta,P=single(qp,qm,.1); H2=bond_sw(sectors(H,G,F),sectors(H,G,F)); h=len(P)
 print('single full coefficients',qp,qm)
 for (i,A),(j,B) in product(cbasis(gam),repeat=2):
  c=np.trace(np.kron(A,B)@H2).real/h**2
  if abs(c)>1e-9: print(i,j,c)
