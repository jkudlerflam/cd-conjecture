"""Exact finite table and all-clock-vertex checks for endpoint_shuttle.py."""
from itertools import product
import json
import numpy as np
from endpoint_shuttle import Shuttle, PHASES
from sweeping_clock_checks import path


def cycles():
    s=Shuttle();cases=[]
    for n in range(8,33):
        for a,b in ((2,4),(2,n-4),(n-6,n-4)):
            if a+1>=b:continue
            c=s.initial(n,a,b);initial=c.snapshot();ticks=0;active=[0,0]
            while True:
                s.step(c);ticks+=1
                if c.flags[0]:assert c.records[0]==a and c.records[1]==a+1;active[0]+=1
                if c.flags[-1]:assert c.records[-1]==b+1 and c.records[-2]==b;active[1]+=1
                if c.snapshot()==initial:break
                assert ticks<20*n
            assert min(active)>0
            assert ticks==10*n+8*a-4*b-12
            cases.append(dict(n=n,a=a,b=b,steps=ticks,active=active))
    return dict(cases=len(cases),steps=sum(c['steps'] for c in cases),largest=cases[-1])


def check_clock(n,a,b):
    s=Shuttle();c=s.initial(n,a,b);direct=s.initial(n,a,b)
    words,steps=path(n);counts=np.zeros((2,2),dtype=int);sweeps=0
    # At every vertex, including after every partial outward and return move.
    def check_vertex():
        for side,port,color,psrc,csrc in ((0,0,1,a,a+1),(1,n-1,n-2,b+1,b)):
            if c.flags[port] and c.position!=port:
                assert c.records[port]==psrc
                counts[side,0]+=1
                if c.position!=color:
                    assert c.records[color]==csrc
                    counts[side,1]+=1
        assert sorted(c.records)==list(range(n))
        # Record transport cannot change the total parity for any grade assignment.
    check_vertex()
    for lab,j in steps:
        if lab=='out':s.bond(c,j)
        elif lab=='back':s.onsite(c,j+1)
        elif lab=='finish':
            s.onsite(c,0);s.step(direct);sweeps+=1
            assert c.snapshot()==direct.snapshot()
        check_vertex()
    assert sweeps==2**n and np.min(counts)>0
    return dict(n=n,clock_vertices=len(words),sweeps=sweeps,port_active=counts[:,0].tolist(),density_active=counts[:,1].tolist())


def local_checks():
    s=Shuttle()
    assert all(s.symbols[i%s.B][0]==s.symbols[j%s.B][0] for i,j in enumerate(s.table))
    assert all(s.inverse[j]==i for i,j in enumerate(s.table))
    # Exhaust both fresh/updated G branches for every control and every pair
    # of two-state records, every neutral flag and static local mark.
    from endpoint_shuttle import Configuration
    total=0
    for k,phase,pos,records,marks,flags in product(range(s.Q),range(2),range(2),product(range(2),repeat=2),product(('L','B','R'),repeat=2),product(range(2),repeat=2)):
        c=Configuration(list(marks),list(flags),list(records),pos,k,phase)
        initial=c.snapshot();parity=sum(c.records)%2
        s.bond(c,0);assert sum(c.records)%2==parity
        s.bond(c,0);assert c.snapshot()==initial;total+=1
    return dict(collision_entries=len(s.table),prescribed=len(s.partial),bond_involution_cases=total)


def code_columns(q,m):
    assert q%4==0 and m<=q*q//8
    s=q//2;blank=q;columns=[]
    for bits in product(range(2),repeat=m):
        col={}
        for side,f in product(range(2),range(m)):
            parity=sum(bits[:f] if side==0 else bits[f+1:])%2
            edge=2*f+parity;a,b=divmod(edge,s)
            port=a if bits[f]==0 else s+b;color=(a+b)%s
            records=(port,color,blank,blank) if side==0 else (blank,blank,color,port)
            env=(side,f,bits[:f]+bits[f+1:])
            key=records+(env,)
            assert key not in col
            col[key]=1
            assert ((port>=s and port<q)+sum(env[2]))%2==sum(bits)%2
        columns.append(col)
    return columns


def compression_checks(q,m):
    cols=code_columns(q,m);D=2**m;s=q//2;den=2*m
    index={}
    for x,col in enumerate(cols):
        for key in col:assert key not in index;index[key]=x
    err=0.;density_err=0.;nentries=0
    for side,f,component in product(range(2),range(m),range(2)):
        mat=np.zeros((D,D),dtype=complex)
        for x,col in enumerate(cols):
            for key in col:
                pidx=0 if side==0 else 3;port=key[pidx]
                for parity in range(2):
                    a,b=divmod(2*f+parity,s);b+=s;sgn=(-1)**parity
                    if port not in (a,b):continue
                    flip=b if port==a else a
                    # Left X,Y. Right Y,-X.
                    if (side,component)==(0,0):v=sgn
                    elif (side,component)==(1,1):v=-sgn
                    else:v=sgn*(1j if port==a else -1j)
                    kp=list(key);kp[pidx]=flip;y=index.get(tuple(kp))
                    if y is not None:mat[y,x]+=v/den
        target=np.zeros_like(mat)
        for x in range(D):
            bits=tuple((x>>(m-1-j))&1 for j in range(m));y=x^(1<<(m-1-f))
            parity=sum(bits[:f] if side==0 else bits[f+1:])%2
            sg=(-1)**parity
            if (side,component)==(0,0):v=sg
            elif (side,component)==(1,1):v=-sg
            else:v=1j*sg*(-1)**bits[f]
            target[y,x]=v/den
        err=max(err,float(np.max(abs(mat-target))));nentries+=np.count_nonzero(mat)
    for side,f in product(range(2),range(m)):
        values=[]
        for x,col in enumerate(cols):
            total=0
            for key in col:
                pi,ci=(0,1) if side==0 else (3,2)
                port,color=key[pi],key[ci]
                for parity in range(2):
                    a,b=divmod(2*f+parity,s)
                    if color==(a+b)%s:
                        if port==a:total+=1
                        if port==s+b:total-=1
            expected=(-1)**((x>>(m-1-f))&1)/den
            density_err=max(density_err,abs(total/den-expected))
    assert err<1e-14 and density_err<1e-14
    return dict(q=q,m=m,logical_dimension=D,code_nonzeros=len(index),port_error=err,density_error=density_err,nonzero_port_entries=int(nentries))

if __name__=='__main__':
    out=dict(local=local_checks(),cycles=cycles(),clock=[check_clock(n,2,n-4) for n in (8,9,10)],
             code=[compression_checks(q,m) for q,m in ((4,2),(8,4),(8,8))])
    with open('endpoint_shuttle_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
