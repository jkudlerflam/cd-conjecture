#!/usr/bin/env python3
"""Exact cyclic-color/vacancy test of material translation gluing.

No color dynamics or field theory is used: only perfect cyclic ordered
words and vacancy swaps. Full ring spectra are compared to the allowed
spinless-fermion flux sectors.
"""
from itertools import combinations
from pathlib import Path
import json
import numpy as np


def run(L,M,b):
    assert M%b==0 and M<L
    positions=list(combinations(range(L),M))
    words=[]
    for x in positions:
        for phase in range(b):
            w=[0]*L
            for j,i in enumerate(x):w[i]=1+(j+phase)%b
            words.append(tuple(w))
    lookup={w:i for i,w in enumerate(words)}
    H=np.zeros((len(words),len(words)))
    U=np.zeros_like(H)
    for col,w in enumerate(words):
        U[lookup[(w[-1],)+w[:-1]],col]=1
        for i in range(L):
            j=(i+1)%L
            if (w[i]==0)!=(w[j]==0):
                v=list(w);v[i],v[j]=v[j],v[i]
                H[lookup[tuple(v)],col]-=1
            elif w[i] and w[j]:
                # This is the literal adjacent-color penalty's kernel.
                assert (w[j]-w[i])%b==1%b
    expected=[];sector_energies=[]
    for r in range(b):
        theta=2*np.pi*r/b
        k=(2*np.pi*np.arange(L)-theta-np.pi*(M-1))/L
        eps=-2*np.cos(k)
        eigenvalues=sorted(sum(eps[i] for i in c) for c in combinations(range(L),M))
        expected.extend(eigenvalues)
        sector_energies.append({"r":r,"twist":float(theta),"energy":float(eigenvalues[0])})
    ev,vec=np.linalg.eigh(H)
    discrepancy=float(np.max(np.abs(ev-np.sort(expected))))
    assert discrepancy<1e-10
    gs=vec[:,ev<ev[0]+1e-10]
    # All tested fixed-M rings have a unique physical ground state.
    ground_translation=gs.conj().T@U@gs
    assert gs.shape[1]==1
    return {"L":L,"M":M,"b":b,"dimension":len(words),"full_spectrum_error":discrepancy,"ground_degeneracy":gs.shape[1],"ground_translation_phase":[float(ground_translation[0,0].real),float(ground_translation[0,0].imag)],"sector_ground_energies":sector_energies}


if __name__=='__main__':
    result=[run(*c) for c in [(5,2,2),(6,2,2),(7,4,2),(5,3,3),(7,3,3),(8,6,3)]]
    Path(__file__).with_name('color_gluing_results.json').write_text(json.dumps({"status":"passed","cases":result},indent=2)+'\n')
    for c in result:print(c['L'],c['M'],c['b'],c['dimension'],c['full_spectrum_error'])
