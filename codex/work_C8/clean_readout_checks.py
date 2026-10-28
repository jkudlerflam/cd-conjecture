"""Finite checks of the actual autonomous readout machine and all clock slices."""
from itertools import product
from collections import Counter
import json
import numpy as np
from clean_readout_machine import Machine,local_op,record_gate,ROLES,STATIC_ALLOWED


def grade(z,q):
    data,ready,v,p,x=z
    return data+x+int(q//2<=v<q)


def local_checks():
    M=Machine();cases=[]
    for q in (4,8,12):
        states=list(product(range(2),range(2),range(q+1),range(2),range(2)))
        count=0
        for op in ('I','PARITY','TAKE','LOAD','ODD_P','EVEN_P','FLAG'):
            for inverse in (False,True):
                outputs=[local_op(op,z,q,inverse) for z in states]
                assert len(set(outputs))==len(states)
                for z,w in zip(states,outputs):
                    assert grade(z,q)==grade(w,q)
                    assert local_op(op,w,q,not inverse)==z
                    count+=1
        for g in ('I','SWAP','EVEN_ADD','ODD_ADD','SUM','INC'):
            for inverse in (False,True):
                states2=list(product(range(q+1),repeat=2))
                outputs=[record_gate(g,*z,q,inverse) for z in states2]
                assert len(set(outputs))==len(states2)
                for z,w in zip(states2,outputs):
                    assert sum(q//2<=v<q for v in z)==sum(q//2<=v<q for v in w)
                    assert record_gate(g,*w,q,not inverse)==z
                    count+=1
        cases.append(dict(q=q,checked_columns=count))
    # Inverse completion is checked for every actual control-role-cursor block.
    assert len(M.table)==928
    for i,(j,op,inv) in enumerate(M.table):
        a,iop,iinv=M.inverse[j]
        assert a==i and iop==op and iinv!=inv
    return dict(cases=cases,complete_blocks=len(M.table),prescribed_blocks=len(M.rules))


def active_check(c,q,bits):
    """Checks exact clean factorization on a basis input, returns active labels."""
    m=len(bits);s=q//2;n=len(c.roles);out=[]
    cursors=[i for i,t in enumerate(c.tape) if t[0]]
    for side,port,color,far in ((0,0,1,2),(1,n-1,n-2,n-3)):
        if not c.tape[port][2] or c.pos==port:continue
        assert len(cursors)==1
        f=cursors[0]-3;p=sum(bits[:f] if side==0 else bits[f+1:])%2
        a,b=divmod(2*f,s)
        assert c.p==p and c.x==0
        expected=a if bits[f]==0 else s+b+p
        assert c.tape[port][3]==expected,(side,f,c.snapshot(),bits)
        assert c.tape[color][3]==(a+b+p)%s
        assert c.tape[far][3]==b
        for j,bit in enumerate(bits):assert c.tape[j+3][1]==(0 if j==f else bit)
        out.append((side,f,int(c.pos!=color)))
    # Product microscopic parity includes the carried quantum bit.
    parity=(sum(t[1]+int(s<=t[3]<q) for t in c.tape)+c.x)%2
    assert parity==sum(bits)%2
    return out


def cycle(q,bits,full_snapshots=False):
    M=Machine();c=M.initialize(q,bits);initial=c.snapshot();cnt=Counter();events=[]
    limit=100*len(bits)*len(c.roles)
    for tick in range(1,limit+1):
        M.step(c,q)
        for side,f,density in active_check(c,q,bits):
            cnt[side,f,'port']+=1
            if density:cnt[side,f,'density']+=1
            if full_snapshots:
                # Erase only the physical endpoint record. Everything else,
                # including the carried bits, must agree for x and x^f.
                port=0 if side==0 else len(c.roles)-1
                tape=[tuple(t) for t in c.tape]
                t=list(tape[port]);t[3]=-1;tape[port]=tuple(t)
                events.append(((tick,side,f), (tuple(tape),c.pos,c.control,c.p,c.x,c.phase)))
        if c.snapshot()==initial:break
    else:raise AssertionError('No cycle')
    assert len(cnt)==4*len(bits)
    # Exact derived full-cycle duration for target m=q^2/8.
    m=len(bits);assert tick==11*m*m+105*m,(m,tick)
    return tick,cnt,events


def cycle_checks():
    cases=[]
    for q in (4,8):
        m=q*q//8;baseline=None;clean={};num=0
        for bits in product(range(2),repeat=m):
            tick,cnt,events=cycle(q,bits,full_snapshots=True)
            if baseline is None:baseline=cnt
            assert cnt==baseline
            for (t,side,f),env in events:
                key=(t,side,f,bits[:f]+bits[f+1:])
                if key in clean:assert clean[key]==env
                else:clean[key]=env
            num+=1
        cases.append(dict(q=q,m=m,inputs=num,cycle_steps=tick,
                          active_port=[baseline[side,0,'port'] for side in range(2)],
                          clean_pairs=len(clean)))
    rng=np.random.default_rng(20261028)
    q=12;m=18;baseline=None
    for _ in range(24):
        bits=tuple(map(int,rng.integers(0,2,m)));tick,cnt,_=cycle(q,bits)
        if baseline is None:baseline=cnt
        assert cnt==baseline
    cases.append(dict(q=q,m=m,inputs=24,cycle_steps=tick))
    return cases


def clock_stream(n):
    """O(n) memory generation of exactly the existing 15-letter legal path.
    Uses its independent forward rule evaluator, not a different clock.
    """
    from sweeping_clock_checks import forward
    w=('I',)+('Z',)*(n-1)
    yield w,None,None
    while True:
        nexts=forward(w)
        assert len(nexts)<=1
        if not nexts:break
        w,lab,j=nexts[0]
        yield w,lab,j


def clock_check(q,bits,full_input=False):
    M=Machine();c=M.initialize(q,bits);direct=M.initialize(q,bits)
    n=len(c.roles);cnt=Counter();sweeps=0;vertices=0;byf=[]
    for w,lab,j in clock_stream(n):
        if lab=='out':M.bond(c,j,q)
        elif lab=='back':M.onsite(c,j+1,q)
        elif lab=='finish':
            M.onsite(c,0,q);M.step(direct,q);sweeps+=1
            assert c.snapshot()==direct.snapshot()
        for side,f,density in active_check(c,q,bits):
            cnt[side,f,'port']+=1
            if density:cnt[side,f,'density']+=1
        vertices+=1
    assert sweeps==2**n
    assert vertices==(2*n+4)*2**n+n-4
    assert len(cnt)==4*len(bits) and min(cnt.values())>0
    return dict(q=q,m=len(bits),n=n,clock_vertices=vertices,sweeps=sweeps,
                port_counts=[[cnt[s,f,'port'] for f in range(len(bits))] for s in range(2)],
                density_counts=[[cnt[s,f,'density'] for f in range(len(bits))] for s in range(2)])


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--clock-q8',action='store_true');args=p.parse_args()
    out=dict(local=local_checks(),cycles=cycle_checks())
    histories=[clock_check(4,bits) for bits in product(range(2),repeat=2)]
    assert all(z==histories[0] for z in histories)
    out['clock']=[histories[0]]
    if args.clock_q8:out['clock'].append(clock_check(8,(1,0,1,1,0,0,1,0)))
    with open('clean_readout_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
