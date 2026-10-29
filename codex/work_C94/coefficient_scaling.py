"""Second-order half-density checks without forming a joint interacting space."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from inherited_ed import selector
from inherited_second_order import source_vectors,susceptibility

def main():
    rows=[]
    for n in [3,6,9]:
        chi,res=susceptibility(n,n);dist=np.minimum(np.arange(n),n-np.arange(n))
        c=-n*chi[0,n//2];regular=chi[0]+c/n
        good=[3]*(n//3);_,y=source_vectors(good,[False]*len(good),selector(3));bg=float(np.sum(chi*(y.T@y)))
        cases=[('one_pin',good,[True]+[False]*(len(good)-1),1)]
        if n>=6:cases += [('length_pair',[2,4]+good[2:],[False]*len(good),6),('long_cell',[6]+good[2:],[False]*(len(good)-1),6)]
        for label,lengths,pins,w in cases:
            _,y=source_vectors(lengths,pins,selector(3));bs=float(np.sum(chi*(y.T@y)))
            rows.append(dict(N=n,V=n,label=label,W=w,delta_B_over_W=(bs-bg)/w,
                linear_residual=res,constant_choice=c,regular_first_moment=float(np.dot(1+dist,abs(regular))),
                note='The constant is a finite-size diagnostic choice, not a proved asymptotic kernel coefficient.'))
        print(n,res,flush=True)
    Path('coefficient_scaling.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':main()
