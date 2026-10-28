import numpy as np,json
from simple_dmrg import DMRG
from scipy.linalg import eigh
results={}
# XX hopping with empty and occupied local states, half filling.
splus=np.array([[0,0],[1,0]],float)
h=-(np.kron(splus,splus.T)+np.kron(splus.T,splus))
for L in [24,48]:
 obj=DMRG(L,h,np.zeros((2,2)))
 history=obj.run([8,16,32,48,64,64],output='calibration_xx_L'+str(L))
 Eexact=-sum(2*np.cos(np.pi*np.arange(1,L//2+1)/(L+1)))
 print('XX energy',L,history[-1]['E'],'exact',Eexact,'difference',history[-1]['E']-Eexact,flush=True)
 results['xx'+str(L)]=dict(Eexact=Eexact,E=history[-1]['E'])
# SU(3) antiferromagnetic permutation chain.
d=3;swap=np.zeros((d*d,d*d))
for a in range(d):
 for b in range(d):swap[d*a+b,d*b+a]=1
for L in [24,48]:
 obj=DMRG(L,swap,np.zeros((3,3)))
 history=obj.run([8,16,32,48,64,64,64,64],output='calibration_su3_L'+str(L))
 results['su3'+str(L)]=dict(E=history[-1]['E'])
with open('calibration_summary.json','w') as f:json.dump(results,f,indent=2)
