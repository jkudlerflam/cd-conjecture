"""Full physical spectra of the factorized first-order length/XX-payload model.

Markers carry an actual binary payload. A marker crossing the coordinate
origin rotates the canonical payload ordering. This bookkeeping is retained
in both the physical graph and the independent free-spectrum formula.
"""
import itertools
import json
import numpy as np
from scipy import sparse as sp
from length_fluid_checks import marker_states


def physical_joint(b,M,r,n,J,g):
    L=b*M+r
    patterns=sorted(s for s,_ in marker_states(b,M,r))
    bits=[tuple(int(j in occ) for j in range(M))
          for occ in itertools.combinations(range(M),n)]
    states=[(p,z) for p in patterns for z in bits]
    index={s:i for i,s in enumerate(states)}
    rows=[];cols=[];vals=[]
    for i,(p,z) in enumerate(states):
        for j,x in enumerate(p):
            if (x-p[j-1])%L==b and (p[(j+1)%M]-x)%L==b+1:
                changed=list(p);changed[j]=(x+1)%L;changed=tuple(sorted(changed))
                zz=(z[-1],)+z[:-1] if x==L-1 else z
                out=index[(changed,zz)]
                rows.extend([out,i]);cols.extend([i,out]);vals.extend([-J,-J])
        for j in range(M):
            jj=(j+1)%M
            if z[j]!=z[jj]:
                zz=list(z);zz[j],zz[jj]=zz[jj],zz[j]
                rows.append(index[(p,tuple(zz))]);cols.append(i);vals.append(-g)
    H=sp.coo_matrix((vals,(rows,cols)),shape=(len(states),)*2).tocsr()
    assert sp.linalg.norm(H-H.getH())<1e-12
    return H


def fermion_sets(M,n):
    out=[]
    for occ in itertools.combinations(range(M),n):
        q=2*np.pi*(np.asarray(occ)+(n-1)/2)/M
        m=(sum(occ)+n*(n-1)//2)%M
        out.append((q,m))
    return out


def joint_levels(b,M,r,n,J,g):
    L=b*M+r;out=[]
    for q,m in fermion_sets(M,r):
        for p,mp in fermion_sets(M,n):
            for ell in range((-m-mp)%M,L,M):
                out.append(-2*J*np.cos(q+2*np.pi*ell/(L*M)).sum()
                           -2*g*np.cos(p).sum())
    return np.sort(out)


def main():
    rows=[]
    for b,M,r,n in [(4,4,2,2),(3,5,2,2),(4,6,3,3)]:
        J,g=1.,.7
        H=physical_joint(b,M,r,n,J,g)
        e=np.linalg.eigvalsh(H.toarray())
        exact=joint_levels(b,M,r,n,J,g)
        assert len(e)==len(exact)
        error=float(np.max(np.abs(e-exact)))
        assert error<1e-10
        ground=-2*(J*np.sin(np.pi*r/M)+g*np.sin(np.pi*n/M))/np.sin(np.pi/M)
        assert abs(e[0]-ground)<1e-10
        rows.append(dict(b=b,M=M,r=r,payload_number=n,L=b*M+r,
                         full_dimension=len(e),complete_spectrum_error=error,
                         ground_energy=float(e[0]),predicted_ground_energy=float(ground)))
    # Equal velocities J=g=1 at half filling let one Casimir coefficient
    # measure total c=2. This is a projected-model check only.
    b=4;M=np.array([40,48,64,80,96,128,160],float);L=(b+.5)*M
    E=4*M/np.pi-4/np.sin(np.pi/M);velocity=2*(b+.5)
    fit=np.linalg.lstsq(np.stack([np.ones_like(L),L**-2,L**-4],axis=1),E/L,rcond=None)[0]
    c=-6*fit[1]/(np.pi*velocity)
    assert abs(c-2)<1e-5
    out=dict(full_spectra=rows,equal_velocity=velocity,joint_casimir_c=float(c),
             limitation='First-order model with one XX payload, not a microscopic counterexample.')
    with open('joint_fluid_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
