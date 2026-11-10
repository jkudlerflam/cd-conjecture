#!/usr/bin/env python3
"""Small-ring algebraic product-eigenstate check, not a thermodynamic claim."""
import json
import math
import analyze_saved as A
import numpy as np
from scipy.linalg import eigvalsh
from check_mappings import product,embed
Z=A.Z

def main():
    rows=[]
    for name,jz,d,current_h in [('R28',1.,1/16,-7/4),('R8.75',5/8,1/10,-7/8)]:
        z=math.sqrt((jz-.5)/(jz+.5));hf=-2*math.sqrt(jz*jz-.25)
        spin=np.array([math.sqrt((1+z)/2),1j*math.sqrt((1-z)/2)])
        for L in [4,6,8]:
            xx=[embed(Z.PX,j,L) for j in range(L)]
            yy=[embed(Z.PY,j,L) for j in range(L)]
            zz=[embed(Z.ZZ_,j,L) for j in range(L)]
            H=sum(.5*xx[j]@xx[(j+1)%L]-.5*yy[j]@yy[(j+1)%L]+jz*zz[j]@zz[(j+1)%L]
                +d*(xx[j]@yy[(j+1)%L]-yy[j]@xx[(j+1)%L]) for j in range(L))+hf*sum(zz)
            v=product([spin.reshape(1,2)]*L).reshape(-1)
            energy=float(np.vdot(v,H@v).real);res=float(np.linalg.norm(H@v-energy*v))
            low=eigvalsh(H,subset_by_index=[0,3]).tolist()
            assert res<1e-12
            rows.append(dict(name=name,L=L,Jz=jz,D=d,current_hz=current_h,factor_hz=hf,
                product_Z=z,product_Y=math.sqrt(1-z*z),E_product=energy,residual=res,
                lowest_energies=low,product_excess=energy-low[0]))
            print(rows[-1],flush=True)
    (A.HERE/'factorization_control.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':main()
