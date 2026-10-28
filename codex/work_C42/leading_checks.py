"""C42 full-code g^2 quadratic and onsite cancellation checks."""
import json
import numpy as np
import scipy.linalg as la
import model
from checks import normal_quadratic

def run():
 out=[]
 for n in [1,2,3]:
  first_hopping=None
  for phi in [np.pi/2,np.pi]:
   H=model.k2(n,.25,.2,0,phi,True)
   quad,M=normal_quadratic(H,n)
   row=dict(n=n,phi_over_pi=phi/np.pi,dimension=H.shape[0],quadratic_residual=float(la.norm(H-quad)),onsite_residual=float(la.norm(M[:n,:n])+la.norm(M[n:,n:])),hopping_norm=float(la.norm(M[:n,n:])))
   if first_hopping is None:first_hopping=M[:n,n:].copy()
   else:
    row['relative_bus_phase_error']=float(la.norm(M[:n,n:]-1j*first_hopping))
    assert row['relative_bus_phase_error']<1e-16
   assert row['quadratic_residual']<1e-16
   assert row['onsite_residual']<1e-16
   out.append(row)
 with open('leading_checks.json','w') as f:json.dump(out,f,indent=2)
 print(json.dumps(out,indent=2))
if __name__=='__main__':run()
