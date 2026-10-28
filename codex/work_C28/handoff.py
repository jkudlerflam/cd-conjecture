"""Save the n=2 spin-chain bond with one common Jordan-Wigner grading."""
import json
import numpy as np
from bus_base import h2
from ff_bus import PHASES
from verify import parameters
from explore import normal_order,roots,scattering

def main():
    n=2;values,K,h=parameters(n);F=np.diag([0,values[1],values[0],sum(values)])
    split=(np.kron(F,np.eye(4))+np.kron(np.eye(4),F))/2
    traceless=lambda A:A-np.trace(A)*np.eye(len(A))/len(A)
    leading=traceless(K[2]+split);quartic=traceless(K[4])
    exact,_=h2(n,g=.02,phases=PHASES)
    full=traceless(exact)/.02**2+traceless(split)
    rr=roots(h,values);R=[r for r in rr if r[1]>0];L=sorted([r for r in rr if r[1]<0],key=lambda r:-r[0])
    pairs=[(i,j) for i in range(4) for j in range(i+1,4)]
    totals=np.array([rr[i][0]+rr[j][0] for i,j in pairs])
    mismatches=[abs(float(np.angle(np.exp(1j*(totals[i]-totals[j]))))) for i in range(6) for j in range(i+1,6)]
    nonzero=[q for q in mismatches if q>1e-8]
    p=[]
    for g in (.2,.1):
        exactg,_=h2(n,g=g,phases=PHASES)
        _,_,U,_=normal_order(exactg,4)
        val=scattering(U,R[1],L[1],R[0],L[0])/g**4
        p.append(dict(g=g,pair_over_g4=[float(val.real),float(val.imag)]))
    np.savez('logical_chain_n2.npz',leading_bond=leading,g4_bond=quartic,
             full_g002_bond=full,mu=values,phases=PHASES,global_scale='lambda^2*g^2/9',basis=np.array(['00','01','10','11']))
    row=dict(quartic_min_nonzero_mismatch=min(nonzero),resonant_pair_differences=len(mismatches)-len(nonzero),
             full_in_g_pair=p,full_g002_minus_through_g4_norm=float(np.linalg.norm(full-leading-.02**2*quartic,2)))
    with open('handoff_results.json','w') as f:json.dump(row,f,indent=2)
    print(json.dumps(row,indent=2))

if __name__=='__main__':main()
