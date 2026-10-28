"""Full-sweep extension of binary_clock_checks. Independent exhaustive rules.
Run OPENBLAS_NUM_THREADS=1 python3 sweeping_clock_checks.py
"""
from itertools import product
import json
import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal

ALPHABET=('L','0','1','C0','C1','R0','R1','I','Z','A0','A1','D0','D1','E0','E1')
COUNTER=('C0','C1','R0','R1')
SWEEP=('D0','D1','E0','E1')
HEADS=set(COUNTER+SWEEP+('I',))
BITS=('0','1')
PAIRS={('L','L'),('I','Z'),('Z','Z')}
PAIRS|={('L',h) for h in COUNTER+('I',)}
PAIRS|={(h,b) for h in COUNTER+SWEEP for b in BITS}
PAIRS|={(a,b) for a in BITS for b in BITS}
PAIRS|={('A'+a,'A'+b) for a in BITS for b in BITS}
PAIRS|={('A'+a,h) for a in BITS for h in SWEEP}
LEFT={'L','A0','A1'}|HEADS
RIGHT=set(BITS)|HEADS|{'Z'}
RULES=[]
for b in BITS:
    RULES.append((('C1',b),('L','C'+b),'carry'))
    RULES.append((('L','R'+b),('R0',b),'return'))
    for a in BITS:
        RULES.append((('D'+a,b),('A'+a,'D'+b),'out'))
        RULES.append((('A'+a,'E'+b),('E'+a,b),'back'))
RULES.append((('I','Z'),('L','I'),'initialize'))


def defect_vector(w):
    return (w[0] not in LEFT,)+tuple((a,b) not in PAIRS for a,b in zip(w,w[1:]))+(w[-1] not in RIGHT,)


def forward(w):
    out=[];n=len(w)
    for j in range(n-1):
        for a,b,label in RULES:
            if w[j:j+2]==a:out.append((w[:j]+b+w[j+2:],label,j))
    for j in range(n):
        if w[j]=='C0':out.append((w[:j]+('R1',)+w[j+1:],'turn',j))
    if w[0] in ('R0','R1'):out.append((('D'+w[0][1],)+w[1:],'begin',0))
    if w[0] in ('E0','E1'):out.append((('C'+w[0][1],)+w[1:],'finish',0))
    if w[-1]=='I':out.append((w[:-1]+('R0',),'activate',n-1))
    if w[-1] in ('D0','D1'):out.append((w[:-1]+('E'+w[-1][1],),'bounce',n-1))
    return out


def path(n):
    w=('I',)+('Z',)*(n-1);words=[w];steps=[];seen={w}
    while True:
        fs=forward(w)
        assert len(fs)<=1
        if not fs:break
        w,label,j=fs[0]
        assert w not in seen and not any(defect_vector(w))
        words.append(w);steps.append((label,j));seen.add(w)
    assert len(words)==(2*n+4)*2**n+n-4
    assert sum(label=='finish' for label,j in steps)==2**n
    assert words[-1]==('L',)*(n-1)+('C1',)
    return words,steps


def full_check(n):
    words=list(product(ALPHABET,repeat=n));idx={w:i for i,w in enumerate(words)}
    adj=[[] for w in words];legal=set()
    for i,w in enumerate(words):
        v=sum(defect_vector(w))
        assert sum(h in HEADS for h in w)<=v+1
        if v==0:legal.add(w)
        for z,lab,j in forward(w):
            assert defect_vector(w)==defect_vector(z)
            k=idx[z];adj[i].append(k);adj[k].append(i)
    for w,nb in zip(words,adj):assert len(nb)<=2*sum(h in HEADS for h in w)
    p,steps=path(n)
    assert set(p)==legal
    seen=set();sizes=[];cycles=0
    for i in range(len(words)):
        if i in seen:continue
        comp=[];todo=[i];seen.add(i)
        while todo:
            j=todo.pop();comp.append(j)
            for k in adj[j]:
                if k not in seen:seen.add(k);todo.append(k)
        sizes.append(len(comp))
        cycles+=sum(len(adj[j]) for j in comp)//2 >= len(comp)
    return dict(n=n,full_words=len(words),legal_vertices=len(legal),components=len(sizes),largest_component=max(sizes),components_with_cycles=int(cycles))


def qca_check(n=3):
    # Work site = two qubits (left/right transport lanes).
    # Outward rule swaps R_j with L_(j+1), so these swaps commute.
    # Return rule applies a random fixed on-site F to the departed site.
    wdim=4;D=wdim**n
    rng=np.random.default_rng(125010)
    X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));X=(X+X.conj().T)/2
    ev,U=eigh(X);F=(U*np.exp(.07j*ev))@U.conj().T
    def site(F,j):return np.kron(np.kron(np.eye(wdim**j),F),np.eye(wdim**(n-j-1)))
    def swap_qubits(i,j):
        P=np.zeros((D,D))
        for k in range(D):
            x=[(k>>(2*n-1-b))&1 for b in range(2*n)]
            x[i],x[j]=x[j],x[i]
            v=sum(bit<<(2*n-1-b) for b,bit in enumerate(x));P[v,k]=1
        return P
    G=[swap_qubits(2*j+1,2*j+2) for j in range(n-1)]
    Uq=np.eye(D,dtype=complex)
    for g in G:Uq=g@Uq
    for j in range(n):Uq=site(F,j)@Uq
    words,steps=path(n);V=np.eye(D,dtype=complex);k=0;error=0.
    for label,j in steps:
        if label=='out':V=G[j]@V
        elif label=='back':V=site(F,j+1)@V
        elif label=='finish':
            V=site(F,0)@V;k+=1
            error=max(error,float(np.linalg.norm(V-np.linalg.matrix_power(Uq,k))))
    assert k==2**n and error<1e-11
    return dict(n=n,work_dimension=D,sweeps=k,qca_cycle_error=error)


def spectra_and_convexity():
    def count(n):return (2*n+4)*2**n+n-4
    def a(n):return 4*np.sin(np.pi/(2*(count(n)+1)))**2
    err=0.
    for n in range(1,7):
        C=count(n)
        e=eigh_tridiagonal(np.zeros(C),-np.ones(C-1),eigvals_only=True)
        err=max(err,float(np.max(abs(e+2*np.cos(np.pi*np.arange(1,C+1)/(C+1))))))
    for n in range(2,35):
        lower=(n+2)*a(n)-(n+1)*a(n+1)
        upper=(n+1)*a(n-1)-n*a(n)
        assert 0<lower<upper
        delta=(lower+upper)/2
        values=[(a(m)-delta)/(m+1) for m in range(1,80)]
        assert np.argmin(values)+1==n
    return dict(max_spectral_n=6,max_vertices=count(6),spectrum_error=err,selected_n_max=34)


if __name__=='__main__':
    out={'all_sectors':[full_check(n) for n in range(1,5)],'qca':qca_check(),'spectra':spectra_and_convexity()}
    with open('sweeping_clock_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
