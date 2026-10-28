"""Build the actual 17-state repeated bond and check its cell restrictions."""
import json
import numpy as np
from scipy import sparse as sp
from spatial_bus import walker_cell,data_ops
from verify_bus import skron


def bond(n,g=.2,lam=.08,mu_sel=0.):
    labels=[(a,b) for a in ('L','R') for b in (0,1)]
    labels += [('X',s,b) for s in ('v','0','+','-') for b in (0,1)]
    labels += [('Me',)]+[('M',s) for s in ('v','0','+','-')]
    ix={x:i for i,x in enumerate(labels)}
    def arrow(x):
        return {'L':('L','L'),'R':('R','R'),'X':('L','R'),
                'Me':('R','L'),'M':('L','L')}[x[0]]
    q=[];mark=[]
    onsite=np.zeros((17,17),complex)
    bare=4-2*np.cos(np.pi/(n+2))
    for i,x in enumerate(labels):
        species=x[1] if x[0] in ('X','M') else 'v'
        bit=x[-1] if x[0] in ('X','L','R') else 0
        q.append(bit+{'v':0,'0':0,'+':1,'-':-1}[species])
        mark.append(int(x[0] in ('M','Me')))
        onsite[i,i]=bare*int(species!='v')+mu_sel*mark[-1]
    for target,source in [(('X','+',0),('X','0',1)),(('X','-',1),('X','0',0))]:
        onsite[ix[target],ix[source]]=g;onsite[ix[source],ix[target]]=g
    h=(np.kron(onsite,np.eye(17))+np.kron(np.eye(17),onsite))/2
    for i,x in enumerate(labels):
        for j,y in enumerate(labels):
            if arrow(x)[1]!=arrow(y)[0]:h[17*i+j,17*i+j]+=8
    phases={'v':0.,'0':.41,'+':.31,'-':-.07}
    def hop(after,before,value):
        row=17*ix[after[0]]+ix[after[1]]
        col=17*ix[before[0]]+ix[before[1]]
        h[row,col]+=value;h[col,row]+=value.conjugate()
    for s in phases:
        for a in (0,1):
            val=-np.exp(1j*phases[s])*(-1 if s=='0' and a else 1)
            for b in (0,1):
                hop((('L',a),('X',s,b)),(('X',s,a),('R',b)),val)
            hop((('L',a),('M',s)),(('X',s,a),('Me',)),val)
    ports=[]
    for typ in ('X','M'):
        fs=[]
        for s in ('0','+','-'):
            F=np.zeros((17,17),complex)
            for b in ((0,1) if typ=='X' else (None,)):
                label=lambda ss:(typ,ss,b) if typ=='X' else (typ,ss)
                F[ix[label(s)],ix[label('v')]]=1
            fs.append(F)
        ports.append((fs[0]+fs[0].conj().T,fs[1]+fs[2].conj().T,fs[2]+fs[1].conj().T))
    X,M=ports
    v=np.kron(M[0],X[0])+1j*np.kron(M[1],X[2])-1j*np.kron(M[2],X[1])
    h+=lam*v
    charge=np.array(q);q2=charge[:,None]+charge[None,:]
    return h,onsite,ports,ix,float(np.linalg.norm((q2.reshape(-1)[:,None]-q2.reshape(-1)[None,:])*h))


def check(n):
    h,on,ports,ix,qerr=bond(n)
    b=n+1;N=2**n
    H=sum((skron([sp.eye(17**j),sp.csr_matrix(h),sp.eye(17**(b-j-2))]) for j in range(b-1)),
          sp.csr_matrix((17**b,17**b)))
    H+=skron([on/2,sp.eye(17**(b-1))])+skron([sp.eye(17**(b-1)),on/2])
    states=[];gauge=[]
    for j,s in [(j,'v') for j in range(b)]+[(j,s) for j in range(b) for s in ('0','+','-')]:
        for z in range(N):
            bits=[(z>>(n-1-i))&1 for i in range(n)]
            word=[('L',a) if i<j else ('X',s,a) if i==j else ('R',a) for i,a in enumerate(bits)]
            word.append(('M',s) if j==n else ('Me',))
            ind=0
            for x in word:ind=17*ind+ix[x]
            states.append(ind)
            gauge.append((-1)**sum(bits[:j]) if s=='0' else 1)
    D=np.diag(gauge)
    restricted=D@H[states][:,states].toarray()@D+2*np.cos(np.pi/(n+2))*np.eye(len(states))
    expected,G,pp,_=walker_cell(n,g=.2,phases=(.41,.31,-.07))
    err=float(np.linalg.norm(restricted-expected))
    perr=0.
    for j in (0,1):
        for a,F in enumerate(ports[j]):
            full=skron([F,sp.eye(17**n)]) if j==0 else skron([sp.eye(17**n),F])
            actual=D@full[states][:,states].toarray()@D
            perr=max(perr,float(np.linalg.norm(actual-pp[j][a])))
    assert err<1e-12 and perr<1e-12 and qerr<1e-12
    return dict(n=n,physical_cell_dimension=17**b,legal_dimension=len(states),
                h_shape=list(h.shape),cell_error=err,port_error=perr,charge_error=qerr,
                hermiticity=float(np.linalg.norm(h-h.conj().T)))


if __name__=='__main__':
    out=[check(n) for n in (1,2)]
    print(json.dumps(out,indent=2))
    with open('literal_bond_results.json','w') as f:json.dump(out,f,indent=2)
