"""C28: positive 17-letter initialized bus core and exact legal-cell restriction.

The vacuum walker is a Laplacian, replacing C22's adjacency ground state.
The excited rails and ports are unchanged. All energies are in t=1 units.
"""
import numpy as np
from scipy import sparse as sp
from scipy.linalg import block_diag
from bus_base import one_bus, data_ops

PHASES=(np.pi/2+.005,0.,0.)

def cell(n,g=.2,d=4.,phases=PHASES):
    N,m=2**n,n+1
    A=2*np.eye(m)-np.diag(np.ones(m-1),1)-np.diag(np.ones(m-1),-1)
    A[0,0]=A[-1,-1]=1
    Hv=np.kron(A,np.eye(N))
    He,sources=one_bus(n,delta=d,g=g,phases=phases,gauged=True)
    H=block_diag(Hv,He)
    G=np.zeros((len(H),N),complex);G[:m*N]=np.kron(np.ones((m,1))/np.sqrt(m),np.eye(N))
    ports=[]
    for end,j in enumerate((0,n)):
        V=np.zeros((m*N,N),complex);V[j*N:(j+1)*N]=np.eye(N)
        fs=[]
        for a in range(3):
            F=np.zeros_like(H);F[m*N:,:m*N]=sources[end][a]@V.conj().T;fs.append(F)
        ports.append((fs[0]+fs[0].conj().T,fs[1]+fs[2].conj().T,fs[1].conj().T+fs[2]))
    return H,G,ports,1/np.sqrt(m)

def charges(n):
    q=np.array([i.bit_count() for i in range(2**n)])
    return np.concatenate([np.tile(q,n+1),np.concatenate([q+a for _ in range(n+1) for a in (0,1,-1)])])

def field(n,values):
    """Legal-space density matrix. Independent entries are physically used only
    at n<=2 here; at larger n the repeated-readout restriction must be imposed.
    """
    cs,_,_,_=data_ops(n);F=sum(v*c.conj().T@c for v,c in zip(values,cs))
    return np.kron(np.eye(4*(n+1)),F),F

def literal(g=.2,d=4.,phases=PHASES,K=2.):
    # L0,L1,R0,R1,Xv0,Xv1,X00,X01,X+0,X+1,X-0,X-1,Me,Mv,M0,M+,M-
    D=17;arrows=[(0,0)]*2+[(1,1)]*2+[(0,1)]*8+[(1,0)]+[(0,0)]*4
    markers=set(range(12,17));exc=set(range(6,12))|set(range(14,17))
    Q=np.diag([int(i in exc) for i in range(D)])
    C=np.zeros((D,D),complex)
    C[8,7]=C[7,8]=g;C[11,6]=C[6,11]=g
    o=(d-2)*Q+C
    h=np.kron(o/2+Q,np.eye(D))+np.kron(np.eye(D),o/2+Q)
    b=np.zeros_like(h);lap=np.zeros_like(h)
    for s in range(4):
        for a in range(2):
            # X_sa R_b -> L_a X_sb
            for bit in range(2):
                src=(4+2*s+a)*D+2+bit;dst=a*D+4+2*s+bit
                phase=1 if s==0 else np.exp(1j*phases[s-1])*(-1 if s==1 and a else 1)
                if s==0:
                    lap[src,src]+=1;lap[dst,dst]+=1;lap[src,dst]-=1;lap[dst,src]-=1
                else:b[dst,src]-=phase;b[src,dst]-=phase.conjugate()
            src=(4+2*s+a)*D+12;dst=a*D+13+s
            phase=1 if s==0 else np.exp(1j*phases[s-1])*(-1 if s==1 and a else 1)
            if s==0:
                lap[src,src]+=1;lap[dst,dst]+=1;lap[src,dst]-=1;lap[dst,src]-=1
            else:b[dst,src]-=phase;b[src,dst]-=phase.conjugate()
    bad=np.diag([float(arrows[i][1]!=arrows[j][0] or i in markers) for i in range(D) for j in range(D)])
    h+=b+lap+K*bad
    left=o/2+Q+K*np.diag([float(i not in {0,1,4,5,6,7,8,9,10,11}) for i in range(D)])
    right=o/2+Q+K*np.diag([float(i not in markers) for i in range(D)])
    return h,left,right

def legal_embedding(n):
    N,m=2**n,n+1;rows=[];cols=[];vals=[]
    # Match vacuum-first and then position/species ordering of cell().
    for j in range(m):
        for s in range(4):
            for bits in range(N):
                data=[(bits>>(n-1-i))&1 for i in range(n)]
                word=[]
                for i,a in enumerate(data):
                    word.append(a if i<j else 4+2*s+a if i==j else 2+a)
                word.append(12 if j<n else 13+s)
                row=0
                for a in word:row=row*17+a
                col=j*N+bits if s==0 else m*N+(3*j+s-1)*N+bits
                # Ungauge the neutral rail's Jordan-Wigner transport.
                sign=(-1)**sum(data[:j]) if s==1 else 1
                rows.append(row);cols.append(col);vals.append(sign)
    return sp.coo_matrix((vals,(rows,cols)),shape=(17**(n+1),4*m*N)).tocsr()
