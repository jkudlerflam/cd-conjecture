#!/usr/bin/env python3
"""Tiny dense checks of exact algebraic mappings; no phase assumptions."""
import json
import numpy as np
import analyze_saved as A
Z=A.Z

def product(ops):
    x=np.array([[1.]])
    for op in ops:x=np.kron(x,op)
    return x

def embed(op,j,L):
    return product([op if i==j else Z.I2 for i in range(L)])

def main():
    L=4;jz=1.;d=1/16;h=-7/4
    xx=[embed(Z.PX,j,L) for j in range(L)]
    yy=[embed(Z.PY,j,L) for j in range(L)]
    zz=[embed(Z.ZZ_,j,L) for j in range(L)]
    H=sum(.5*xx[j]@xx[j+1]-.5*yy[j]@yy[j+1]+jz*zz[j]@zz[j+1]
          +d*(xx[j]@yy[j+1]-yy[j]@xx[j+1]) for j in range(L-1))+h*sum(zz)
    U=product([Z.PX if j%2 else Z.I2 for j in range(L)])
    mapped=sum(.5*(xx[j]@xx[j+1]+yy[j]@yy[j+1])-jz*zz[j]@zz[j+1]
          -(-1)**j*d*(xx[j]@yy[j+1]+yy[j]@xx[j+1]) for j in range(L-1))
    mapped+=h*sum((-1)**j*zz[j] for j in range(L))
    R=np.zeros_like(H)
    for s in range(2**L):R[int(f'{s:0{L}b}'[::-1],2),s]=1
    V=product([np.diag(np.exp(-.5j*((-1)**j*np.pi/2)*np.array([1,-1]))) for j in range(L)])
    Qstag=sum((-1)**j*zz[j] for j in range(L))
    Hzero=H-sum(d*(xx[j]@yy[j+1]-yy[j]@xx[j+1]) for j in range(L-1))
    errors=dict(odd_X_rotation=float(np.max(abs(U@H@U.conj().T-mapped))),
        reflection_conjugation=float(np.max(abs(R@H.conj()@R-H))),
        staggered_rotation_conjugation=float(np.max(abs(V@H.conj()@V.conj().T-H))),
        D_zero_staggered_Z_conservation=float(np.max(abs(Hzero@Qstag-Qstag@Hzero))))
    assert max(errors.values())<1e-12
    (A.HERE/'mapping_checks.json').write_text(json.dumps(errors,indent=2)+'\n')
    print(errors)

if __name__=='__main__':main()
