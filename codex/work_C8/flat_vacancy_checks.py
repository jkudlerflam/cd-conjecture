"""Uniform-vacancy wrapper: complete sectors, gaps and closed-form weights."""
import json
from itertools import combinations
import numpy as np
from arrow_wrapper_checks import local_bond,legal_words
from selector_wrapper_checks import sparse_open
from length_filter_checks import filter_data,occupied_sector
from periodic_clean_bond import PeriodicBond
from periodic_bond_checks import small_wrapper_comparison,actual_columns,history_compression,charge_form_identity


def exclusion(N,p,complete=False):
    states=list(combinations(range(N),p));index={s:i for i,s in enumerate(states)}
    H=np.zeros((len(states),len(states)))
    edges=list(combinations(range(N),2)) if complete else [(j,j+1) for j in range(N-1)]
    for i,s in enumerate(states):
        occ=set(s)
        for a,b in edges:
            if (a in occ)==(b in occ):continue
            t=tuple(sorted(occ.symmetric_difference((a,b))))
            H[i,i]+=1;H[index[t],i]-=1
    return H


def gap_checks():
    cases=0;minimum_ratio=1e9;comparison_error=0.
    for N in range(2,10):
        for p in range(1,N):
            H=exclusion(N,p);ev=np.linalg.eigvalsh(H)
            assert abs(ev[0])<1e-12 and ev[1]>=1/N**2-1e-12
            minimum_ratio=min(minimum_ratio,float(ev[1]*N*N));cases+=1
            if N<=6:
                C=exclusion(N,p,True)
                comparison_error=max(comparison_error,max(0.,float(-np.linalg.eigvalsh(N**3*H-C)[0])))
    assert comparison_error<1e-11
    return dict(cases=cases,maximum_N=9,minimum_gap_times_N_squared=minimum_ratio,
                comparison_matrix_negativity=comparison_error)


def sectors(k):
    phase,h,a,b,eps=filter_data(k);out=[]
    for N in (k+1,k+2):
        rows=[];bound=2.
        for r in range(2,N+1):
            H,holes,chi,energy=occupied_sector(N,r,1.,1.,1.,h,a,b)
            ev=np.linalg.eigvalsh(H)
            gap=float(ev[1] if r==k else ev[0])
            if r==N:lower=eps
            else:
                gamma=1/(r*(2*r-1)) if r==k else eps
                aa=(r/N)*gamma;delta=1/N**2;M=r+1
                lower=aa*delta/(aa+M+delta)
            assert gap>=lower-1e-12
            if r==k:assert abs(ev[0])<1e-12
            rows.append(dict(occupied=r,dimension=H.shape[0],gap_or_rejection_energy=gap,lower_bound=lower))
            bound=min(bound,lower)
        # This retains every legal arrow/tag/occupation sector, not just r=k.
        local=local_bond(2,h,a,b,1.,1.,2.,6.)[0]
        words=legal_words(2,N);full=sparse_open(words,local,9)
        ev=np.linalg.eigvalsh(full.toarray())
        assert abs(ev[0]+4)<2e-11 and ev[1]-ev[0]>=bound-2e-11
        out.append(dict(length=N,sectors=rows,complete_legal_dimension=len(words),
                        ground_energy=float(ev[0]),complete_gap=float(ev[1]-ev[0]),
                        complete_gap_lower_bound=bound))
    return dict(k=k,lengths=out)


if __name__=='__main__':
    out=dict(exclusion=gap_checks(),wrappers=[sectors(3),sectors(4)])
    with open('flat_vacancy_checks.json','w') as f:
        json.dump(out,f,indent=2);f.write('\n')
    H=PeriodicBond(flat_vacancies=True)
    oracle=dict(parameters=H.report(),wrapper=small_wrapper_comparison(H),
                columns=actual_columns(H),compression=history_compression(H),
                charge_form_identity=charge_form_identity(H))
    with open('flat_periodic_bond_checks.json','w') as f:
        json.dump(oracle,f,indent=2);f.write('\n')
    print(json.dumps(out,indent=2))
    print(json.dumps(oracle,indent=2))
