import numpy as np
from scipy.linalg import expm
for k in range(1,9):
    S=k/2;d=k+1;m=np.arange(S,-S-1,-1)
    sp=np.zeros((d,d))
    for j in range(1,d):sp[j-1,j]=np.sqrt(S*(S+1)-m[j]*(m[j]+1))
    sx=(sp+sp.T)/2;sy=(sp-sp.T)/(2j);sz=np.diag(m)
    X=sum(np.kron(a,a) for a in [sx,sy,sz]);ev,U=np.linalg.eigh(X)
    J=np.rint((np.sqrt(1+8*(ev+S*(S+1)))-1)/2).astype(int)
    H=np.array([sum(1/i for i in range(1,int(j)+1)) for j in J]);h=(U*H)@U.conj().T
    y=expm(1j*np.pi*sy);v=np.kron(np.eye(d),y);hp=v@h@v.conj().T
    K=(-hp).reshape(d,d,d,d).transpose(0,2,1,3).reshape(d*d,d*d)
    eig=np.linalg.eigvalsh(K)
    print(k,[(round(float(x),10),int(sum(abs(eig-x)<1e-8))) for x in np.unique(np.round(eig,10))])
