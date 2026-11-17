"""Literal four-coupling family; full endpoint fields; source C190 solver."""
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'codex/work_C190'))
import z2dmrg as Z
import models as Mo
import analyze190 as Ana


def coefficients(V,kappa=1.,a=1.,hz=None):
    hz=-V/2+a/V if hz is None else hz
    return dict(Jx=.5,Jy=-.5,Jz=V/4,D=kappa/(4*V),hz=hz,
                V=V,kappa=kappa,a=V*(hz+V/2),Delta=1.,t2=1/(2*V))


def build(cfg):
    L=cfg['L']
    if cfg.get('legs'):
        legs=cfg['legs']; n=L//legs; N=cfg.get('N_leg',round(.19*n))
        mu=-np.cos(np.pi*N/(n+1))-np.cos(np.pi*(N+1)/(n+1))
        mpo,gph=Mo.xx_control(L,float(mu),legs=legs)
        return mpo,gph,0.,dict(legs=legs,mu=float(mu),N_exact=N*legs)
    b=coefficients(cfg['V'],cfg.get('kappa',1.),cfg['a'],cfg.get('hz'))
    h2=Z.pauli_bond(**{k:b[k] for k in ('Jx','Jy','Jz','D')})
    gph=Z.gauge_phases(L,'stag')
    h2s=[Z.rotate_bond(h2,gph[j],gph[j+1]) for j in range(L-1)]
    assert all(np.isrealobj(h) for h in h2s)
    shift=-(b['Jz']+b['hz'])
    h1=b['hz']*Z.ZZ_+shift*Z.I2
    return Z.MPO(L,h2s,[h1]*L),gph,shift,b


def product_expectation(psi,ops):
    E=np.ones((1,1),complex)
    for j in range(max(ops)+1):E=Z._apply_site(E,psi.M[j],ops.get(j))
    return complex(np.trace(E))


def pair_diagnostics(psi,gph):
    L=psi.L; i=L//4
    create=lambda j:np.conj(gph[j])*Z.SM
    destroy=lambda j:gph[j]*Z.SP
    ai=product_expectation(psi,{i:create(i),i+2:create(i+2)})
    rows=[]
    for r in (4,6,8,12,16,24,32):
        j=i+r
        if j+2>=3*L//4:continue
        aj=product_expectation(psi,{j:destroy(j),j+2:destroy(j+2)})
        c=product_expectation(psi,{i:create(i),i+2:create(i+2),j:destroy(j),j+2:destroy(j+2)})-ai*aj
        rows.append(dict(r=r,connected_pair=[c.real,c.imag]))
    n=np.diag([0.,1.]); mid=L//2
    nn=product_expectation(psi,{mid:n,mid+1:n})
    return dict(anchor=i,separation_two_pair_creation=[ai.real,ai.imag],
                connected_pair_correlations=rows,central_adjacent_density=nn.real)
