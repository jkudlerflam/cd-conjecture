import numpy as np
from scipy.linalg import eigvalsh
q=2; d=q+1
T=np.zeros((d*d,d*d)); S=[np.zeros((d,d),complex) for _ in range(3)]
S[0][1,2]=S[0][2,1]=.5
S[1][1,2]=-.5j; S[1][2,1]=.5j
S[2][1,1]=.5; S[2][2,2]=-.5
for a in range(1,d):T[a,d*a]=T[d*a,a]=1
B=1.5*T-(8/3)*sum(np.kron(s,s) for s in S)
for L in (2,3,4,5,6):
    H=sum(np.kron(np.kron(np.eye(d**j),B),np.eye(d**(L-j-2))) for j in range(L-1))
    print(L, eigvalsh(H,subset_by_index=(d**L-1,d**L-1))[0]/(L-1))
