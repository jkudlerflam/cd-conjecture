"""Resolved phase-fluid and full-cell selector diagnostics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import bond16 as b
from checks import ring_matrix, maximum


def compositions(total, n):
    if n == 1:
        yield (total,)
    else:
        for k in range(total+1):
            for tail in compositions(total-k,n-1):
                yield (k,)+tail


def phase_checks():
    out=[]
    for N in (4,5,6):
        L=2*N
        states=list(compositions(N,N)); ix={w:i for i,w in enumerate(states)}
        H=np.zeros((len(states),len(states)))
        for col,w in enumerate(states):
            for j in range(N):
                if w[j]==0:
                    continue
                for k in ((j-1)%N,(j+1)%N):
                    z=list(w); z[j]-=1; z[k]+=1
                    H[ix[tuple(z)],col]-=1
        vals,u=la.eigh(H)
        gs=u[:,0]; E=vals[0]
        A=np.array([w[0]==0 for w in states],float)
        alpha=N/L-np.sin(np.pi*N/L)**2/(N*L*np.sin(np.pi/L)**2)
        source=(A-alpha)*gs
        coeff=u[:,1:].T@source
        chi=float(np.sum(coeff**2/(vals[1:]-E)))
        assert abs(E+2*np.sin(np.pi*N/L)/np.sin(np.pi/L)) < 1e-11
        assert abs(gs@(A*gs)-alpha) < 1e-11
        costs=[]
        for J in (.04,.02,.01):
            energy=la.eigvalsh(H+J*np.diag(A),subset_by_index=(0,0))[0]
            costs.append(dict(J=J,cost=float(energy-E),
                              relaxation_over_J2=float((J*alpha-energy+E)/J**2)))
        out.append(dict(L=L,N=N,dimension=len(states),E0=float(E),
                        E0_formula_error=abs(float(E+2*np.sin(np.pi*N/L)/np.sin(np.pi/L))),
                        alpha=alpha,contact_error=abs(float(gs@(A*gs)-alpha)),
                        chi00=chi,defects=costs))
    return out


def selector_checks():
    # Complete cell Hilbert spaces, never project out forbidden pairs.
    def aa(k):return 2*(1-np.cos(np.pi/(k+1)))
    def dd(k):return (k+1)*aa(k)-k*aa(k+1)
    mu=2-(dd(2)+dd(1))/2
    p=b.params(1,K=1000.,musel=mu,field=0.)
    marker=[i+1 for i,x in enumerate(b.LET) if x[0]=='M']
    data=[i+1 for i,x in enumerate(b.LET) if x[0]!='M']
    out=[]
    for lam in (0.,.005,.01):
        h=b.lift(b.occupied_bond(p,lam),t=0,J=1)
        record=dict(lambda_=lam,parameters=p,sectors=[])
        for places in ((1,3),(0,3)):
            states=list(itertools.product(*[marker if j in places else data for j in range(4)]))
            H=ring_matrix(states,h)
            e,u=eigsh(H,k=1,which='SA',tol=1e-11,ncv=60,
                       v0=np.cos(np.arange(len(states))+.2),maxiter=15000)
            record['sectors'].append(dict(marker_positions=list(places),
                 dimension=len(states),energy=float(e[0]),
                 residual=maximum(H@u-u*e)))
        record['wrong_minus_selected']=record['sectors'][1]['energy']-record['sectors'][0]['energy']
        assert record['wrong_minus_selected'] > 0
        assert max(s['residual'] for s in record['sectors']) < 1e-8
        out.append(record)
    return out


def main():
    out=dict(status='numerical_small_resolved_diagnostics',
             phase_fluid=phase_checks(),full_cell_selector=selector_checks())
    Path('auxiliary_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
