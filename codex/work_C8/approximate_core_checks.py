"""Positive non-FF cores: low-band, angle, and endpoint checks."""
import json
import numpy as np
from vacancy_parent import charge_hamiltonian


def embed(op,j,k,width):
    return np.kron(np.kron(np.eye(2**j),op),np.eye(2**(k-j-width)))


def run_case(p,tau):
    k=3;U=1.;v=.6
    swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]],float)
    num=np.diag([0.,1.])
    h=(np.eye(4)-swap)/2+tau*(np.kron(num,np.eye(2))+np.kron(np.eye(2),num))/2
    a=b=tau*num/2
    boundary=embed(a,0,k,1)+embed(b,k-1,k,1)
    core=boundary.copy()
    for j in range(k-1):core+=embed(h,j,k,2)
    V=np.zeros((2**k,2));V[0,0]=1
    for j in range(k):V[1<<j,1]=1/np.sqrt(k)
    assert np.linalg.norm(core@V-V@np.diag([0,tau]))<1e-12
    gamma=np.linalg.eigvalsh(core)[2]
    holes,hc,chi,ec=charge_hamiltonian(k,p,U,v)
    H=np.kron(hc-ec[0]*np.eye(len(holes)),np.eye(2**k))
    H+=np.kron(np.eye(len(holes)),boundary)
    probs=[]
    for j in range(k-1):
        adjacent=[]
        for hs in holes:
            occupied=[i for i in range(k+p) if i not in hs]
            adjacent.append(int(occupied[j+1]==occupied[j]+1))
        probs.append(float(np.dot(chi*chi,adjacent)))
        H+=np.kron(np.diag(adjacent),embed(h,j,k,2))
    W0=np.kron(chi[:,None],V)
    E,Ue=np.linalg.eigh(H)
    Up=Ue[:,:2]
    left,singular,right=np.linalg.svd(Up.T@W0)
    W=Up@left@right
    angle=np.linalg.norm(W@W.T-W0@W0.T,2)
    delta=ec[1]-ec[0];M=(k-1)*np.linalg.norm(h,2)+np.linalg.norm(a,2)+np.linalg.norm(b,2)
    g0=delta*min(probs)*gamma/((np.sqrt(delta)+np.sqrt(M))**2+min(probs)*gamma)
    bound=min(1.,np.sqrt(M*tau)/(g0-tau))
    assert tau<g0 and E[0]>-1e-12 and E[1]<=tau+1e-12 and E[2]>=g0-1e-12
    assert angle<=bound+1e-12
    X=np.array([[0.,1.],[1.,0.]])
    first_occ=np.diag([int(0 not in hs) for hs in holes])
    F=np.kron(first_occ,embed(X,0,k,1))
    endpoint=np.linalg.norm(W.T@F@W-W0.T@F@W0,2)
    assert endpoint<=2*np.sqrt(2)*angle+1e-12
    return dict(p=p,epsilon=tau,dimension=len(E),gamma=float(gamma),
                predicted_complement_threshold=float(g0),lowest_levels=E[:3].tolist(),
                code_angle=float(angle),code_angle_bound=float(bound),
                endpoint_error=float(endpoint),endpoint_angle_bound=float(2*np.sqrt(2)*angle))


if __name__=='__main__':
    out=[run_case(p,tau) for p in (1,2) for tau in (1e-6,1e-7,1e-8)]
    with open('approximate_core_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
