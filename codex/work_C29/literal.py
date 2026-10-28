"""Literal 26-letter nearest-neighbor bond and legal-sector restrictions."""
import json
import numpy as np
from scipy import sparse as sp
from model import cell,parameters,target_data
from ed import kron


def literal(n,g=.2,lam=.04,eta=10.,mu_sel=0.):
    species=['v','A0','A+','A-','B0','B+','B-']
    labels=[(a,b) for a in ('L','R') for b in (0,1)]
    labels += [('X',s,b) for s in species for b in (0,1)]
    labels += [('Me',)]+[('M',s) for s in species]
    ix={x:i for i,x in enumerate(labels)};d0=len(labels)
    def arrow(x):return {'L':('L','L'),'R':('R','R'),'X':('L','R'),'Me':('R','L'),'M':('L','L')}[x[0]]
    on=np.zeros((d0,d0),complex);qs=[]
    pars=parameters(n);alpha=np.sqrt(2/(n+2))*np.sin(np.pi/(n+2))
    T,mu,W=target_data(n,g)
    for i,x in enumerate(labels):
        s=x[1] if x[0] in ('X','M') else 'v'
        bit=x[-1] if x[0] in ('L','R','X') else 0
        qs.append(bit+(1 if s.endswith('+') else -1 if s.endswith('-') else 0))
        if s!='v':on[i,i]=pars[0 if s[0]=='A' else 1][0]-2*np.cos(np.pi/(n+2))+lam*eta
        if x[0] in ('M','Me'):on[i,i]+=mu_sel
    for r in ('A','B'):
        for target,source in [(('X',r+'+',0),('X',r+'0',1)),(('X',r+'-',1),('X',r+'0',0))]:
            on[ix[target],ix[source]]=g;on[ix[source],ix[target]]=g
    h=(np.kron(on,np.eye(d0))+np.kron(np.eye(d0),on))/2
    # Endpoint fields: at n=1 use the first-site field only.
    for i,x in enumerate(labels):
        for j,y in enumerate(labels):
            if arrow(x)[1]!=arrow(y)[0]:h[d0*i+j,d0*i+j]+=8
            if x[0] in ('M','Me') and y[0] in ('L','R','X'):
                h[d0*i+j,d0*i+j]+=lam**2*alpha**4*g*g*mu[0]*y[-1]
            if n==2 and y[0] in ('M','Me') and x[0] in ('L','R','X'):
                h[d0*i+j,d0*i+j]+=lam**2*alpha**4*g*g*mu[1]*x[-1]
    def hop(after,before,val):
        row=d0*ix[after[0]]+ix[after[1]];col=d0*ix[before[0]]+ix[before[1]]
        h[row,col]+=val;h[col,row]+=val.conjugate()
    for s in species:
        if s=='v':bulk=terminal=0.
        else:
            r=0 if s[0]=='A' else 1;a={'0':0,'+':1,'-':2}[s[1]]
            terminal=pars[r][1][-1,a];bulk=pars[r][1][0,a] if n>1 else 0.
        for a in (0,1):
            sign=-1 if s.endswith('0') and a else 1
            for b in (0,1):hop((('L',a),('X',s,b)),(('X',s,a),('R',b)),-np.exp(1j*bulk)*sign)
            hop((('L',a),('M',s)),(('X',s,a),('Me',)),-np.exp(1j*terminal)*sign)
    ports=[[],[]]
    for p,typ in enumerate(('X','M')):
        for r in ('A','B'):
            fs=[]
            for a in ('0','+','-'):
                F=np.zeros((d0,d0),complex)
                for b in ((0,1) if typ=='X' else (None,)):
                    label=lambda s:(typ,s,b) if typ=='X' else (typ,s)
                    F[ix[label(r+a)],ix[label('v')]]=1
                fs.append(F)
            ports[p].append((fs[0]+fs[0].conj().T,fs[1]+fs[2].conj().T,fs[2]+fs[1].conj().T))
    X,M=ports
    for r,(_,_,phi,w) in enumerate(pars):
        h+=lam*w*(np.kron(M[r][0],X[r][0])+np.exp(1j*phi)*np.kron(M[r][1],X[r][2])+np.exp(-1j*phi)*np.kron(M[r][2],X[r][1]))
    q=np.array(qs);q2=(q[:,None]+q).reshape(-1)
    return h,on,ports,ix,float(np.linalg.norm((q2[:,None]-q2)*h))


def check(n):
    # Boundary bonds are cut. Add endpoint onsite halves and the omitted left field explicitly.
    lam=.04;g=.2;eta=10.
    h,on,ports,ix,qerr=literal(n,g,lam,eta);d0=len(ix);b=n+1;N=2**n
    H=sum((kron([sp.eye(d0**j),h,sp.eye(d0**(b-j-2))]) for j in range(b-1)),sp.csr_matrix((d0**b,d0**b)))
    H+=kron([on/2,sp.eye(d0**(b-1))])+kron([sp.eye(d0**(b-1)),on/2])
    states=[];gauge=[]
    order=[(j,'v') for j in range(b)]+[(j,r+a) for r in ('A','B') for j in range(b) for a in ('0','+','-')]
    for j,s in order:
        for z in range(N):
            bits=[(z>>(n-1-i))&1 for i in range(n)]
            word=[('L',a) if i<j else ('X',s,a) if i==j else ('R',a) for i,a in enumerate(bits)]
            word.append(('M',s) if j==n else ('Me',))
            ind=0
            for x in word:ind=d0*ind+ix[x]
            states.append(ind);gauge.append((-1)**sum(bits[:j]) if s.endswith('0') else 1)
    gauge=np.array(gauge)
    restricted=gauge[:,None]*H[states][:,states].toarray()*gauge[None,:]+2*np.cos(np.pi/(n+2))*np.eye(len(states))
    expected,G,pp,alpha,E,q=cell(n,g)
    T,mu,W=target_data(n,g)
    left=np.diag([mu[0]*((z>>(n-1))&1) for z in range(N)])
    restricted+=lam**2*alpha**4*g*g*np.kron(np.eye(len(states)//N),left)
    expected+=lam*eta*E+lam**2*alpha**4*g*g*np.kron(np.eye(len(states)//N),W)
    err=float(np.linalg.norm(restricted-expected));perr=0.
    for end in (0,1):
        for r in range(2):
            for a,F in enumerate(ports[end][r]):
                full=kron([F,sp.eye(d0**n)]) if end==0 else kron([sp.eye(d0**n),F])
                actual=gauge[:,None]*full[states][:,states].toarray()*gauge[None,:]
                perr=max(perr,float(np.linalg.norm(actual-pp[end][r][a])))
    assert err<1e-12 and perr<1e-12 and qerr<1e-12
    return dict(n=n,d0=d0,bond_dimension=d0*d0,physical_cell_dimension=d0**b,legal_cell_dimension=len(states),cell_error=err,port_error=perr,charge_error=qerr,hermiticity=float(np.linalg.norm(h-h.conj().T)))

if __name__=='__main__':
    out=[check(n) for n in (1,2)]
    print(json.dumps(out,indent=2))
    with open('literal_results.json','w') as f:json.dump(out,f,indent=2)
