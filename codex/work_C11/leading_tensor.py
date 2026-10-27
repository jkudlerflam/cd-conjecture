from c11_calculation import *
from itertools import combinations
for qp,qm in [(1,1),(1,2)]:
 q,r,h,gam,eta,P,F,om,ops=ingredients(qp,qm); pairs=list(combinations(range(r),2)); J=np.zeros((r,r,len(pairs)),complex); L=J.copy()
 for a in range(r):
  for b in range(r):
   for z,(c,d) in enumerate(pairs):
    J[a,b,z]=(a==c and b==d)-(a==d and b==c)
    L[a,b,z]=np.trace(F[c]@F[a]@F[b]@F[d]-F[d]@F[a]@F[b]@F[c])/q
 J=J.reshape(r*r,-1); L=L.reshape(r*r,-1)
 K=119/4*J.T@J-7/4*(J.T@L+L.T@J)-.5*L.T@L
 charge=-.5*np.einsum('aaz->z',L.reshape(r,r,-1))/1j
 print('q,r',q,r,'pairs',pairs,'K / (lambda^2 g^4/Delta^5)\n',np.real_if_close(K),'onsite i gamma gamma /(lambda^2 g^2/Delta^3)',np.real_if_close(charge))
