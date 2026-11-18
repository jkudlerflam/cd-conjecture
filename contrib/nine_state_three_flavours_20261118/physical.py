"""The five controls are literal nine-state nearest-neighbour operators."""
from pathlib import Path
import sys
import numpy as np
from scipy import sparse

PREVIOUS=Path(__file__).resolve().parents[1]/'ten_state_resonance_audit_20261117'
sys.path.insert(0,str(PREVIOUS))
from literal import LETTERS,CHARGE,MARKER,parts as inherited_parts,cell as inherited_cell
from mechanism import fermions

CONTROL_NUMERATORS=np.array([-3385,9387,-7486,-2560,11999])
CONTROL_DENOMINATOR=500000
CONTROL_NAMES=('mu','eta','boundary','kappa_L','kappa_R')
CONTROLS=CONTROL_NUMERATORS/CONTROL_DENOMINATOR


def parts():
    p=inherited_parts();del p['swap'];ix={x:i for i,x in enumerate(LETTERS)}
    for side in ('L','R'):
        H=np.zeros((81,81));a=9*ix[side+'0']+ix[side+'1'];b=9*ix[side+'1']+ix[side+'0']
        H[a,b]=H[b,a]=1;p['swap_'+side]=H
    p['boundary']=np.diag([int(MARKER[a] and LETTERS[b]=='L1') for a in range(9) for b in range(9)])
    return p


def bond(Delta,mu_sel,K,g,lam,t0,controls=CONTROLS):
    mu,eta,b,l,r=controls;p=parts()
    exc=np.array([x in ('DA','DB','MA','MB') for x in LETTERS])
    number=np.array([x in ('L1','R1') for x in LETTERS]);left=np.array([x=='L1' for x in LETTERS])
    o=np.diag(Delta*exc+mu_sel*MARKER+t0*(mu*number+eta*left))
    return (p['walk']+g*p['convert']+lam*p['port']+K*p['penalty']
            +t0*(b*p['boundary']+l*p['swap_L']+r*p['swap_R'])
            +(np.kron(o,np.eye(9))+np.kron(np.eye(9),o))/2)


def profiles(n):
    s=n+2;v=np.sqrt(2/(s+1))*np.sin(np.pi*np.arange(1,s+1)/(s+1))
    p=np.array([sum(v[i+1:]**2) for i in range(n)])
    L=np.diag([sum(v[i+2:]**2) for i in range(n-1)],1);L+=L.T
    R=np.diag([sum(v[:i+1]**2) for i in range(n-1)],1);R+=R.T
    b=np.zeros((n,n));b[0,0]=p[0]
    return (np.eye(n),np.diag(p),b,L,R)


def onsite(controls=CONTROLS,n=3):
    return sum(c*A for c,A in zip(controls,profiles(n)))


def verify_projections(n=3):
    c=inherited_cell(n);words=c['words'];ix={x:i for i,x in enumerate(LETTERS)};index={w:i for i,w in enumerate(words)}
    G=c['G'];f=fermions(n);errors={};p=parts()
    operators={'mu':c['number'],'eta':c['left']}
    for side in ('L','R'):
        local=p['swap_'+side];H=sparse.lil_matrix((len(words),len(words)))
        for k,w in enumerate(words):
            for j in range(len(w)-1):
                source=ix[w[j]]*9+ix[w[j+1]]
                for dest in np.flatnonzero(local[:,source]):
                    z=w[:j]+(LETTERS[dest//9],LETTERS[dest%9])+w[j+2:]
                    if z in index:H[index[z],k]+=local[dest,source]
        operators['kappa_'+side]=H.tocsr()
    # Every legal cell ends in one of the three counted marker states.
    # Project the actual cross-boundary 81x81 operator on the full pair code.
    pairG=np.kron(G,G)
    diag=np.array([p['boundary'][9*ix[a[-1]]+ix[b[0]],9*ix[a[-1]]+ix[b[0]]] for a in words for b in words])
    boundary_pair=pairG.T@(diag[:,None]*pairG)
    target_boundary=sum(profiles(n)[2][i,j]*(f[i].T@f[j]).toarray() for i in range(n) for j in range(n))
    errors['boundary_pair_full_Fock']=float(np.max(abs(boundary_pair-np.kron(np.eye(1<<n),target_boundary))))
    for name,h in zip(CONTROL_NAMES,profiles(n)):
        if name=='boundary':continue
        target=sum(h[i,j]*(f[i].T@f[j]).toarray() for i in range(n) for j in range(n))
        errors[name]=float(np.max(abs(G.T@operators[name]@G-target)))
    assert max(errors.values())<1e-12
    return errors
