"""A uniform local massive gap does not gap the global bare-Q compression.

Exact independent two-level ancillas, reduced to their permutation-symmetric
sector. This tests a failed elimination shortcut, not a failure of EFT.
"""
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal


def case(L,delta=1.,g=.5):
    omega=np.sqrt(delta*delta+4*g*g)
    e=(delta-omega)/2
    p=(1+delta/omega)/2
    k=np.arange(1,L+1,dtype=float)
    diag=delta*k-L*e
    off=-g*np.sqrt(np.arange(2,L+1)*np.arange(L-1,0,-1))
    low=float(eigh_tridiagonal(diag,off,select='i',select_range=(0,0))[0][0])
    bound=-L*e*p**L/(1-p**L)
    assert low>0 and low<=bound+2e-12
    return dict(L=L,full_massive_gap=float(omega),bare_Q_gap=low,
                variational_upper_bound=float(bound),bare_P_ground_weight=float(p**L))


if __name__=='__main__':
    out=[case(L) for L in (4,8,16,32,64,96,128)]
    with open('massive_projection_checks.json','w') as f:
        json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps(out,indent=2))
