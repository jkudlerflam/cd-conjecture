from ks_test import *
for N in [12,16,20,22]:
 o=Orbits(N,'even');d=density(N,'ising',.35); h=total(d,N)
 H=o.matrix(h,0); e,g=lowest(H,2); T=o.matrix(h,2); et,v=lowest(T,10)
 hv=o.matrix(total(d,N,2),0,2)@g[:,0]; w=abs(v.conj().T@hv)**2
 print('N',N,'ground',e[0])
 for i in range(6): print(i,et[i]-e[0],w[i],8*w[i]/(et[i]-e[0])**2)
