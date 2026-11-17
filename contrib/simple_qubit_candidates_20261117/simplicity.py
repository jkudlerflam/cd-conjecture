#!/usr/bin/env python3
"""Counts and coefficient spreads in two fixed common onsite bases."""
import json
import numpy as np
from common import HERE,Z,coefficients


def main():
    u=np.diag([np.exp(-1j*np.pi/8),np.exp(1j*np.pi/8)]);U=np.kron(u,u);rows=[]
    for V,a in [(40,40*(-19.976047250613814+20)),(16,.96),(12,.96),(8,.96),(6,.96),(4,.96),(4,1.)]:
        b=coefficients(V,a=a);h=Z.pauli_bond(**{k:b[k] for k in ['Jx','Jy','Jz','D']})
        transformed=(.5+b['D'])*np.kron(Z.PX,Z.PY)+(.5-b['D'])*np.kron(Z.PY,Z.PX)+b['Jz']*np.kron(Z.ZZ_,Z.ZZ_)
        err=float(np.max(abs(U@h@U.conj().T-transformed)));assert err<1e-13
        original=[.5,-.5,b['D'],-b['D'],b['Jz'],b['hz']];rot=[.5+b['D'],.5-b['D'],b['Jz'],b['hz']]
        rows.append(dict(V=V,a=a,D=b['D'],hz=b['hz'],t2=b['t2'],original_pauli_monomials=6,rotated_pauli_monomials=4,
            original_spread=max(abs(x) for x in original)/min(abs(x) for x in original),
            rotated_spread=max(abs(x) for x in rot)/min(abs(x) for x in rot),
            small_exchange_difference=2*b['D'],rotation_matrix_error=err))
    out=dict(independent_couplings=4,after_overall_scale=3,
        rotation='common onsite U=exp(-i pi Z/8); all candidates',rows=rows,
        caution='Smaller rotated spread encodes the hopping scale in the difference of two exchanges; no dynamical hierarchy is removed.')
    (HERE/'simplicity.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in rows:print(row)


if __name__=='__main__':main()
