"""Constant-alphabet initialized ripple clock. Only numpy/scipy, no network.
Run OPENBLAS_NUM_THREADS=1 python3 binary_clock_checks.py
"""
from itertools import product
import json
import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal

ALPHABET = ('L', '0', '1', 'C0', 'C1', 'R0', 'R1', 'I', 'Z')
HEADS = frozenset(('C0', 'C1', 'R0', 'R1', 'I'))
BITS = ('0', '1')
PAIRS = {('L', 'L'), ('Z', 'Z')}
PAIRS |= {('L', h) for h in HEADS}
PAIRS |= {(h, b) for h in HEADS - {'I'} for b in BITS}
PAIRS |= {(a, b) for a in BITS for b in BITS}
PAIRS.add(('I', 'Z'))
LEFT = {'L'} | HEADS
RIGHT = set(BITS) | HEADS | {'Z'}
# An ordered rule is a forward clock move; the Hamiltonian also contains its inverse.
RULES = []
for b in BITS:
    RULES.append((('C1', b), ('L', 'C'+b), 'carry'+b))
    RULES.append((('L', 'R'+b), ('R0', b), 'return'+b))
RULES.append((('I', 'Z'), ('L', 'I'), 'initialize'))


def defects(w):
    return int(w[0] not in LEFT) + int(w[-1] not in RIGHT) + sum((a,b) not in PAIRS for a,b in zip(w,w[1:]))


def defect_vector(w):
    return (w[0] not in LEFT,) + tuple((a,b) not in PAIRS for a,b in zip(w,w[1:])) + (w[-1] not in RIGHT,)


def forward(w):
    n = len(w)
    out = []
    for j in range(n-1):
        for a,b,label in RULES:
            if w[j:j+2] == a:
                out.append((w[:j]+b+w[j+2:], label, j))
    for j in range(n):
        if w[j] == 'C0':
            out.append((w[:j]+('R1',)+w[j+1:], 'turn', j))
    if w[0] in ('R0','R1'):
        out.append((('C'+w[0][1],)+w[1:], 'restart', 0))
    if w[-1] == 'I':
        out.append((w[:-1]+('R0',), 'activate', n-1))
    return out


def legal_words(n):
    words = []
    for j in range(n):
        words.append(('L',)*j+('I',)+('Z',)*(n-j-1))
        for h in ('C0','C1','R0','R1'):
            for tail in product(BITS, repeat=n-j-1):
                words.append(('L',)*j+(h,)+tail)
    return words


def ordered_path(n):
    w = ('I',)+('Z',)*(n-1)
    path = [w]
    steps = []
    seen = {w}
    while True:
        nxt = forward(w)
        assert len(nxt) <= 1
        if not nxt:
            break
        w,label,j = nxt[0]
        assert w not in seen and defects(w) == 0
        path.append(w)
        steps.append((label,j))
        seen.add(w)
    assert seen == set(legal_words(n))
    assert w == ('L',)*(n-1)+('C1',)
    assert len(path) == n + 4*(2**n-1)
    return path,steps


def enumerate_full(n):
    words = list(product(ALPHABET, repeat=n))
    index = {w:i for i,w in enumerate(words)}
    adj = [[] for w in words]
    for i,w in enumerate(words):
        assert sum(s in HEADS for s in w) <= defects(w)+1
        for v,label,j in forward(w):
            assert defect_vector(w) == defect_vector(v)
            k = index[v]
            adj[i].append(k)
            adj[k].append(i)
    for w,ns in zip(words,adj):
        assert len(ns) <= 2*sum(s in HEADS for s in w)
    seen = set()
    comps = []
    for i in range(len(words)):
        if i in seen:
            continue
        todo = [i]
        seen.add(i)
        comp = []
        while todo:
            a = todo.pop()
            comp.append(a)
            for b in adj[a]:
                if b not in seen:
                    seen.add(b)
                    todo.append(b)
        comps.append(comp)
    legal = [c for c in comps if defects(words[c[0]]) == 0]
    assert len(legal) == 1 and len(legal[0]) == n+4*(2**n-1)
    # Small full component spectra also test illegal cycles and multi-head sectors.
    min_illegal_ff = np.inf
    min_illegal_hop = np.inf
    if n <= 4:
        for comp in comps:
            if defects(words[comp[0]]) == 0:
                continue
            idx = {a:k for k,a in enumerate(comp)}
            A = np.zeros((len(comp),len(comp)))
            for a in comp:
                for b in adj[a]:
                    A[idx[a],idx[b]] += 1
            v = defects(words[comp[0]])
            min_illegal_ff = min(min_illegal_ff, eigh(np.diag(A.sum(axis=1))/2-A/2+8*v*np.eye(len(comp)), subset_by_index=(0,0))[0][0])
            min_illegal_hop = min(min_illegal_hop, eigh(-A+8*v*np.eye(len(comp)),subset_by_index=(0,0))[0][0])
        assert min_illegal_ff >= 8-1e-10
        assert min_illegal_hop >= 4-1e-10
    return dict(n=n, full_words=len(words), components=len(comps), largest_component=max(map(len,comps)), legal_vertices=len(legal[0]), illegal_ff_min=None if n>4 else float(min_illegal_ff), illegal_hop_min=None if n>4 else float(min_illegal_hop))


def convex_check():
    # Stable trigonometric form avoids cancellation at large n.
    def a(x):
        count = x+4*(2.**x-1)
        return 4*np.sin(np.pi/(2*(count+1)))**2
    records=[]
    for n in range(2,40):
        lower=(n+2)*a(n)-(n+1)*a(n+1)
        upper=(n+1)*a(n-1)-n*a(n)
        assert 0 < lower < upper
        delta=(lower+upper)/2
        energies=np.array([(a(m)-delta)/(m+1) for m in range(1,100)])
        assert np.argmin(energies)+1 == n
        records.append(dict(n=n,delta=delta,selected_density=float(energies[n-1]),crossing_width=upper-lower))
    return records


def clock_spectra():
    records=[]
    for n in range(1,9):
        path,steps=ordered_path(n)
        C=len(path)
        exact_ff=1-np.cos(np.pi*np.arange(C)/C)
        diag=np.ones(C);diag[[0,-1]]=.5
        ff=eigh_tridiagonal(diag,-.5*np.ones(C-1),eigvals_only=True)
        hop=eigh_tridiagonal(np.zeros(C),-np.ones(C-1),eigvals_only=True)
        exact_hop=-2*np.cos(np.pi*np.arange(1,C+1)/(C+1))
        err=max(np.max(abs(ff-exact_ff)), np.max(abs(hop-exact_hop)))
        assert err<2e-12
        records.append(dict(n=n,vertices=C,spectrum_error=float(err),ff_gap=float(ff[1]),hop_ground=float(hop[0])))
    return records


def initialized_quantum_check():
    # q=4 per site = one uninitialized logical qubit and one initialized ancilla.
    # Noncommuting parity-preserving gates dress all counter edges; I sweeps are idle.
    n=2;q=4;D=q**n
    path,steps=ordered_path(n);C=len(path)
    rng=np.random.default_rng(51027)
    parity_local=np.diag([1,1,-1,-1])
    def unitary(dim,parity):
        A=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim))
        A=(A+A.conj().T)/2
        A=(A+parity@A@parity)/2
        ev,V=eigh(A)
        return (V*np.exp(.19j*ev))@V.conj().T
    gate={name:unitary(D,np.kron(parity_local,parity_local)) for name in ('carry0','carry1','return0','return1')}
    gate['turn']=unitary(q,parity_local)
    gate['restart']=unitary(q,parity_local)
    gate['activate']=np.eye(q)
    p=np.diag([0,1,0,1])
    H=np.zeros((C*D,C*D),complex)
    W=[np.eye(D,dtype=complex)]
    for t,(label,j) in enumerate(steps):
        if label=='initialize': U=np.eye(D)
        elif label in ('turn','restart','activate'):
            U=np.kron(gate[label],np.eye(q)) if j==0 else np.kron(np.eye(q),gate[label])
        else: U=gate[label]
        a=slice(t*D,(t+1)*D);b=slice((t+1)*D,(t+2)*D)
        H[a,a]+=.5*np.eye(D);H[b,b]+=.5*np.eye(D)
        H[b,a]-=.5*U;H[a,b]-=.5*U.conj().T
        W.append(U@W[-1])
    for t,w in enumerate(path):
        for j,letter in enumerate(w):
            if letter=='I':
                pin=np.kron(p,np.eye(q)) if j==0 else np.kron(np.eye(q),p)
                H[t*D:(t+1)*D,t*D:(t+1)*D]+=pin
    B=np.eye(D)[:,[0,2,8,10]]
    V=np.vstack([w@B for w in W])/np.sqrt(C)
    ev=eigh(H,eigvals_only=True)
    bound=min(1/(2*C),1/(4*C*C))
    assert np.count_nonzero(abs(ev)<1e-10)==4
    assert ev[4]>=bound-1e-12
    residual=np.linalg.norm(H@V)
    assert residual<1e-12
    assert np.linalg.norm(V.conj().T@V-np.eye(4))<1e-12
    return dict(n=n,clock_vertices=C,physical_legal_dimension=C*D,ground_multiplicity=4,measured_gap=float(ev[4]),proved_gap_lower=bound,history_residual=float(residual))


if __name__=='__main__':
    out={'all_sectors':[enumerate_full(n) for n in range(1,6)],'spectra':clock_spectra(),'selection':convex_check(),'quantum_initialization':initialized_quantum_check()}
    with open('binary_clock_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps({k:(v if k=='quantum_initialization' else v[-1]) for k,v in out.items()},indent=2))
