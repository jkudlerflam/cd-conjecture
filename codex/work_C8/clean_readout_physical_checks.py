"""Direct physical compression, all-work-space primitives, and input-format tests."""
from itertools import product
import json
import numpy as np
from clean_readout_machine import Machine,ROLES,STATIC_ALLOWED
from clean_readout_checks import clock_stream,grade


def compress_history(q=4):
    M=Machine();m=q*q//8;inputs=list(product(range(2),repeat=m));D=2**m
    configs=[M.initialize(q,x) for x in inputs];n=m+6;s=q//2
    accum=np.zeros((2,m,2,D,D),dtype=complex)
    fields=np.zeros((2,m,D,D),dtype=float);count=0
    for w,lab,j in clock_stream(n):
        for c in configs:
            if lab=='out':M.bond(c,j,q)
            elif lab=='back':M.onsite(c,j+1,q)
            elif lab=='finish':M.onsite(c,0,q)
        keys={c.snapshot():i for i,c in enumerate(configs)}
        assert len(keys)==D
        for x,c in enumerate(configs):
            for side,port,color in ((0,0,1),(1,n-1,n-2)):
                if not c.tape[port][2] or c.pos==port:continue
                old=c.tape[port][3]
                for f,p in product(range(m),range(2)):
                    a,b=divmod(2*f+p,s);sg=(-1)**p
                    if old not in (a,s+b):continue
                    new=s+b if old==a else a
                    c.tape[port][3]=new
                    y=keys.get(c.snapshot())
                    c.tape[port][3]=old
                    if y is not None:
                        xv=sg;yv=1j*sg*(1 if old==a else -1)
                        accum[side,f,0,y,x]+=xv if side==0 else yv
                        accum[side,f,1,y,x]+=yv if side==0 else -xv
                    if c.pos!=color and c.tape[color][3]==(a+b)%s:
                        fields[side,f,x,x]+=1 if old==a else -1
        count+=1
    scales=[];errors=[]
    for side,f in product(range(2),range(m)):
        target=[]
        for comp in range(2):
            A=np.zeros((D,D),complex)
            for x,bits in enumerate(inputs):
                y=x^(1<<(m-1-f));p=sum(bits[:f] if side==0 else bits[f+1:])%2
                X=(-1)**p;Y=1j*X*(-1)**bits[f]
                A[y,x]=(X if comp==0 else Y) if side==0 else (Y if comp==0 else -X)
            target.append(A)
        N=round(float(np.vdot(target[0],accum[side,f,0]).real/D))
        for comp in range(2):errors.append(float(np.max(abs(accum[side,f,comp]-N*target[comp]))))
        Z=np.diag([(-1)**bits[f] for bits in inputs]);B=round(float(np.trace(Z@fields[side,f])/D))
        errors.append(float(np.max(abs(fields[side,f]-B*Z))))
        assert N>0 and B>0
        scales.append(dict(side=side,f=f,port_count=N,density_count=B))
    assert max(errors)==0
    return dict(q=q,m=m,clock_vertices=count,logical_dimension=D,integer_matrix_error=max(errors),scales=scales)


def format_checks():
    # Count the entire immutable-role language with the actual two-site
    # allowed-pair projector and endpoint role projectors.
    counts=[]
    for n in range(1,31):
        paths={r:int(r=='L') for r in ROLES}
        for _ in range(n-1):paths={b:sum(paths[a] for a in ROLES if (a,b) in STATIC_ALLOWED) for b in ROLES}
        assert paths['R']==int(n>=7)
        counts.append(paths['R'])
    # I-sweep pins permit exactly the stated one or two onsite work columns.
    # On a fixed role, enumerate cursor,data,READY,record and the full direct
    # sum of ordinary + (control,p,x,phase) head states.
    M=Machine();q=4;heads=[None]+list(product(M.names,range(2),range(2),range(2)))
    allowed={};rejected=0
    for role in ROLES:
        total=0
        for cursor,data,ready,v,head in product(range(2),range(2),range(2),range(q+1),heads):
            good=(cursor==int(role=='S') and ready==0 and
                  (role in ('S','D') or data==0) and
                  v==(0 if role in ('A','B','E','F') else q) and
                  head==(('L0',0,0,0) if role=='L' else None))
            total+=good;rejected+=not good
        assert total==(2 if role in ('S','D') else 1)
        allowed[role]=total
    return dict(role_lengths_tested=30,accepted_columns_by_role=allowed,rejected_onsite_columns=rejected,
                initialized_dimension_q4_n8=4)


def whole_space_checks():
    M=Machine();rng=np.random.default_rng(28);count=0
    from clean_readout_machine import Config
    for q in (4,8):
        for _ in range(6000):
            c=Config(list(rng.choice(ROLES,2)),[[int(rng.integers(2)),int(rng.integers(2)),int(rng.integers(2)),int(rng.integers(q+1))] for _ in range(2)],
                     int(rng.integers(2)),str(rng.choice(M.names)),int(rng.integers(2)),int(rng.integers(2)),int(rng.integers(2)))
            initial=c.snapshot();parity=(sum(t[1]+int(q//2<=t[3]<q) for t in c.tape)+c.x)
            M.bond(c,0,q)
            assert (sum(t[1]+int(q//2<=t[3]<q) for t in c.tape)+c.x)==parity
            M.bond(c,0,q);assert c.snapshot()==initial
            j=c.pos;M.onsite(c,j,q)
            assert (sum(t[1]+int(q//2<=t[3]<q) for t in c.tape)+c.x)==parity
            M.onsite(c,j,q);assert c.snapshot()==initial;count+=1
    return dict(random_whole_space_GF_involution_tests=count)

if __name__=='__main__':
    out=dict(physical_compression=compress_history(),format=format_checks(),whole_space=whole_space_checks())
    with open('clean_readout_physical_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
