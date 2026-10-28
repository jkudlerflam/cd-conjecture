import itertools,json,numpy as np
from momentum_ed import OrbitHamiltonian
from qutrit_observables import entropy_profile,cfit

out={}
for L in [9,12,15,18,21,24,36,48]:
 N=L//3;mom=2*np.pi*np.arange(N)/L
 full=np.exp(1j*np.outer(np.arange(L),mom))/np.sqrt(L)
 prof=[]
 for l in range(1,L//2+1):
  occ=np.linalg.eigvalsh(full[:l]@full[:l].conj().T);occ=np.clip(occ,1e-15,1-1e-15)
  prof.append(dict(cut=l,S=float(-sum(occ*np.log(occ)+(1-occ)*np.log(1-occ))),S2=float(-sum(np.log(occ**2+(1-occ)**2)))))
 out['xx_'+str(L)]=dict(entropy=prof,fit=cfit(prof,L))
 print('XX',L,out['xx_'+str(L)]['fit'],flush=True)
swap=np.zeros((9,9))
for a in range(3):
 for b in range(3):swap[3*a+b,3*b+a]=1
for L in [9,12,15]:
 n=L//3;codes=[];sites=set(range(L))
 for twos in itertools.combinations(range(L),n):
  rem=sorted(sites-set(twos));base=sum(2<<(2*j) for j in twos)
  for ones in itertools.combinations(rem,n):codes.append(base+sum(1<<(2*j) for j in ones))
 states=np.sort(np.asarray(codes,dtype=np.uint64))
 ob=OrbitHamiltonian(L,L,custom_basis=states,custom_h=swap)
 es,vs,row=ob.solve(0,nev=2,tol=1e-11)
 states,v=ob.reconstruct(0,vs[:,0]);prof=entropy_profile(states,v,L)
 out['su3_'+str(L)]=dict(E=float(es[0]),entropy=prof,fit=cfit(prof,L),residual=row['residual'])
 print('SU3',L,out['su3_'+str(L)]['fit'],es.tolist(),flush=True)
 with open('periodic_calibration_results.json','w') as f:json.dump(out,f,indent=2)
