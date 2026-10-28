"""C54 small algebra controls. No source imports, network, or large ED."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path
import numpy as np


def lattice(n):
    N=2*n
    vectors={(0,)*N}
    for k in range(1,4):
        multis=list(combinations_with_replacement(range(N),k))
        counts=[]
        for m in multis:
            x=[0]*N
            for i in m: x[i]+=1
            counts.append(x)
        for x in counts:
            for y in counts:
                z=tuple(a-b for a,b in zip(x,y))
                if sum(a*a for a in z)<=6: vectors.add(z)
    hist=Counter()
    for z in vectors:
        for h in range(-1,2):
            d=F(sum(a*a for a in z),2)+h*h
            if not 0<d<=3: continue
            hist[str(d)]+=1
            if d<=2:
                nz=sorted(a for a in z if a)
                if h==0:
                    assert nz in ([-1,1],[-1,-1,1,1])
                else:
                    assert nz in ([],[-1,1])
    expected_bilinears=N*(N-1)
    expected_quartics=6*math.comb(N,4)
    assert hist['1']==expected_bilinears+2
    assert hist['2']==expected_quartics+2*expected_bilinears
    return dict(sorted(hist.items()))

out={'vertex_counts':{str(n):lattice(n) for n in range(1,5)}}
eta=F(1,4)
assert 9*(1-eta)>4*(1+eta)
out['dimension_bound']={'metric_factor':math.sqrt(3/5),'d3_min':3*math.sqrt(3/5),'d4_min':4*math.sqrt(3/5)}
# Exact C47 normalization; K is an explicitly assumed XXZ benchmark.
t=F(9,32); u=F(243,11648)
delta=2*u/t
assert delta==F(27,182)
K=math.pi/(2*(math.pi-math.acos(float(delta))))
out['xxz_assumed_benchmark']={'delta':str(delta),'K':K,'umklapp_dimension':4*K,'relative_metric_eta_before_d4_reaches_2':(4*K*K-1)/(4*K*K+1)}
# Scalar attractive length chain: h - V(n_i+n_j)/2 is positive for V<=-2t.
length=[]
for V in (-3.,-2.,-1.9):
    t0=1.
    n=np.diag([0.,1.]); eye=np.eye(2)
    a=np.array([[0.,1.],[0.,0.]])
    h=-t0*(np.kron(a.T,a)+np.kron(a,a.T))+V*np.kron(n,n)
    shifted=h-V/2*(np.kron(n,eye)+np.kron(eye,n))
    ev=np.linalg.eigvalsh(shifted)
    assert np.allclose(ev,sorted([0,0,-V/2-t0,-V/2+t0]))
    length.append({'V_over_t':V,'shifted_eigenvalues':ev.tolist()})
out['length_attraction']=length
# Imported rounded C44 diagnostics, not certified parameters.
out['diagnostic_attraction_ratios']=[{'n':n,'V_over_2t':v/(2*t0)} for n,v,t0 in [(1,-.00069091,1.97e-8),(2,-.00034838,4.30e-10)]]
# Finite-dimensional example: weak norms do not control a zero-T coefficient.
# A staggered free mass m has continuum d=1 and exact single-particle gap 2m.
k=np.linspace(-math.pi/2,math.pi/2,3001)
m=.01;t0=1.
energies=np.sqrt(4*t0*t0*np.cos(k)**2+m*m)
assert abs(2*energies.min()-2*m)<1e-12
out['weak_bound_countermodel']={'mass':m,'single_particle_band_gap':float(2*energies.min()),'symmetry_caveat':'breaks one-cell translation; not a C39 counterexample'}
# Exponential winding tail beats the oscillator spacing for every fixed coefficient.
out['winding_M2_exp_minus_M']={str(M):M*M*math.exp(-M) for M in (10,20,40,80)}
Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
