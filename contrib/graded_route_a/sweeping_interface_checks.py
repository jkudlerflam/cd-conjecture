"""Independent full-sweep graph and ordered-gate interface checks."""
import itertools
import json
from pathlib import Path
import numpy as np

ALPHABET=('L','0','1','C0','C1','R0','R1','I','Z','A0','A1','D0','D1','E0','E1')
COUNTER=('C0','C1','R0','R1')
SWEEP=('D0','D1','E0','E1')
HEADS=set(COUNTER+SWEEP+('I',))
BITS=('0','1')
PAIRS={('L','L'),('I','Z'),('Z','Z')}
PAIRS|={('L',h) for h in COUNTER+('I',)}
PAIRS|={(h,b) for h in COUNTER+SWEEP for b in BITS}
PAIRS|=set(itertools.product(BITS,repeat=2))
PAIRS|={('A'+a,'A'+b) for a in BITS for b in BITS}
PAIRS|={('A'+a,h) for a in BITS for h in SWEEP}
LEFT={'L','A0','A1'}|HEADS
RIGHT=set(BITS)|HEADS|{'Z'}
RULES=[]
for b in BITS:
    RULES.extend([(('C1',b),('L','C'+b),'carry'),(('L','R'+b),('R0',b),'return')])
    for a in BITS:
        RULES.extend([(('D'+a,b),('A'+a,'D'+b),'out'),
                      (('A'+a,'E'+b),('E'+a,b),'back')])
RULES.append((('I','Z'),('L','I'),'init'))


def defects(w):
    return (w[0] not in LEFT,)+tuple(pair not in PAIRS for pair in zip(w,w[1:]))+(w[-1] not in RIGHT,)


def moves(w):
    out=[]
    for j in range(len(w)-1):
        for before,after,label in RULES:
            if w[j:j+2]==before:
                out.append((w[:j]+after+w[j+2:],label,j))
    for j,x in enumerate(w):
        if x=='C0':out.append((w[:j]+('R1',)+w[j+1:],'turn',j))
    if w[0] in ('R0','R1'):out.append((('D'+w[0][1],)+w[1:],'begin',0))
    if w[0] in ('E0','E1'):out.append((('C'+w[0][1],)+w[1:],'finish',0))
    if w[-1]=='I':out.append((w[:-1]+('R0',),'activate',len(w)-1))
    if w[-1] in ('D0','D1'):out.append((w[:-1]+('E'+w[-1][1],),'bounce',len(w)-1))
    return out


def path(n):
    word=('I',)+('Z',)*(n-1)
    words=[word];steps=[]
    while moves(word):
        choices=moves(word);assert len(choices)==1
        word,label,j=choices[0]
        words.append(word);steps.append((label,j))
    assert len(set(words))==len(words)==(2*n+4)*2**n+n-4
    assert sum(label=='finish' for label,j in steps)==2**n
    return words,steps


def graph_check(n):
    legal=set();transitions=0
    for word in itertools.product(ALPHABET,repeat=n):
        dw=defects(word)
        if not any(dw):legal.add(word)
        for after,label,j in moves(word):
            assert dw==defects(after)
            transitions+=1
    words,steps=path(n)
    assert set(words)==legal
    return dict(n=n,all_words=15**n,legal_words=len(legal),all_directed_transitions=transitions,full_sweeps=2**n)


def circuit_check():
    n=3;rng=np.random.default_rng(71028)
    def unitary(d):
        h=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));h=h+h.conj().T
        ev,u=np.linalg.eigh(h)
        return (u*np.exp(.23j*ev))@u.conj().T
    gate=unitary(4);onsite=unitary(2)
    gs=[np.kron(gate,np.eye(2)),np.kron(np.eye(2),gate)]
    fs=[np.kron(np.eye(2**j),np.kron(onsite,np.eye(2**(n-j-1)))) for j in range(n)]
    update=np.eye(8,dtype=complex)
    for g in gs:update=g@update
    for f in fs:update=f@update
    _,steps=path(n)
    actual=np.eye(8,dtype=complex);cycles=0;error=0.
    for label,j in steps:
        if label=='out':actual=gs[j]@actual
        if label=='back':actual=fs[j+1]@actual
        if label=='finish':
            actual=fs[0]@actual;cycles+=1
            error=max(error,float(np.linalg.norm(actual-np.linalg.matrix_power(update,cycles))))
    assert error<1e-11
    return dict(n=n,work_dimension=8,cycles=cycles,
                overlapping_gate_commutator=float(np.linalg.norm(gs[0]@gs[1]-gs[1]@gs[0])),
                maximum_cycle_error=error)


if __name__=='__main__':
    result={'graph':[graph_check(n) for n in range(1,5)],'ordered_circuit':circuit_check()}
    Path(__file__).with_name('sweeping_interface_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
