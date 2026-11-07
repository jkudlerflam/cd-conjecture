"""Prepare the L=9 SU(3) trial costs; the E18 diagonalization stays on SLURM."""
import json
import numpy as np
from checks import ROOT,ground,bond,trial_parameters,optimize_angles,coherent

if __name__=='__main__':
    L=9;d=3
    e,c,v,res,dim=ground(d,L,0.)
    full=np.zeros(d**L);full[c]=v
    a,b,_,_,error=trial_parameters(full,d,L,bond(d,0.))
    fractional,uniform,angles,grad=optimize_angles(L,a,b)
    coh,purities=coherent(full,d,L,a)
    out=dict(status='numerical trial costs, exact doubled energy pending',model='SU3',L=L,E_L=e,
             dimension=dim,residual=res,bare=2*a,uniform=uniform,fractional=fractional,
             coherent=coh,a=a,b=b,angles=angles,gradient=grad,formula_error=error,
             purities=purities,missing='E_18 from the unsubmitted SLURM run')
    (ROOT/'su3_L9_trial.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['angles','purities']},indent=2))
