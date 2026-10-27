import numpy as np

def calc(N,literal=False):
 k=(2*np.arange(N)+1)*np.pi/N; q=4*np.pi/N; kp=(k-q)%(2*np.pi)
 def eig(k):
  M=np.zeros((len(k),2,2),complex); M[:,0,1]=2j*(1-np.exp(-1j*k));M[:,1,0]=M[:,0,1].conj()
  return np.linalg.eigh(M)
 e,u=eig(k); ep,up=eig(kp)
 M=np.zeros((N,2,2),complex)
 c1=1 if literal else .5*(1+np.exp(-1j*q))
 c2=np.exp(1j*q) if literal else .5*(1+np.exp(1j*q))
 M[:,0,1]=2j*(1-np.exp(-1j*k)*c2); M[:,1,0]=-2j*(1-np.exp(1j*k)*c1)
 F=np.einsum('ni,nij,nj->n',up[:,:,1].conj(),M,u[:,:,0])
 w=.5*abs(F)**2; gap=ep[:,1]-e[:,0]
 m0=w.sum();m1=(w*gap).sum()
 return dict(N=N,c2=8*m0**3/m1**2,m0N2=m0*N*N,m1N3=m1*N**3,cMV=8*sum(w[k<q])/gap[k<q][0]**2)
if __name__=='__main__':
 for n in (8,20,40,100,1000,10000,100000): print(calc(n))
 print('limit',9/98,'m0N2',4*np.pi**2,'m1N3',224*np.pi**3/3)
