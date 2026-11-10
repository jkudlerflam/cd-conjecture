#!/usr/bin/env python3
"""Finite-ring degenerate perturbation audit, with no large-L inference."""
import json
import math
import analyze_saved as A
import numpy as np
from scipy.linalg import eigh,eigvalsh
from check_mappings import embed
Z=A.Z

def main():
    rows=[]
    for name,jz,d in [('R28',1.,1/16),('R8.75',5/8,1/10)]:
        hf=-2*math.sqrt(jz*jz-.25)
        for L in [4,6,8]:
            xx=[embed(Z.PX,j,L) for j in range(L)];yy=[embed(Z.PY,j,L) for j in range(L)]
            zz=[embed(Z.ZZ_,j,L) for j in range(L)];Qz=sum(zz)
            V=sum(xx[j]@yy[(j+1)%L]-yy[j]@xx[(j+1)%L] for j in range(L))
            H0=sum(.5*xx[j]@xx[(j+1)%L]-.5*yy[j]@yy[(j+1)%L]+jz*zz[j]@zz[(j+1)%L]
                for j in range(L))+hf*Qz
            e,U=eigh(H0);ng=int(np.sum(abs(e-e[0])<1e-10));P=U[:,:ng];Q=U[:,ng:]
            W1=P.conj().T@V@P;Wz=P.conj().T@Qz@P
            B=Q.conj().T@V@P;W2=-B.conj().T@(B/(e[ng:]-e[0])[:,None])
            first=eigvalsh(W1);second=eigvalsh(W2);field=eigvalsh(Wz)
            row=dict(name=name,L=L,Jz=jz,D=d,factor_hz=hf,E0=float(e[0]),ground_dimension=ng,
                gap_above_ground=float(e[ng]-e[0]),DM_first_order_eigenvalues=first.tolist(),
                total_Z_eigenvalues=field.tolist(),DM_second_order_eigenvalues=second.tolist(),
                D_squared_second_order_shift=float(d*d*second[0]),
                exact_shift=float(eigvalsh(H0+d*V,subset_by_index=[0,0])[0]-e[0]))
            rows.append(row);print(name,L,'dim',ng,'gap',row['gap_above_ground'],
                'firstwidth',first[-1]-first[0],'second',row['D_squared_second_order_shift'],
                'exact',row['exact_shift'],flush=True)
    (A.HERE/'factorization_projection.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':main()
