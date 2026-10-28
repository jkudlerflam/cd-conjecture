"""Compute all clock visibility counts without traversing its exponential path."""
from math import gcd
import json
import numpy as np
from clean_readout_machine import Machine


def counts(q):
    M=Machine();m=q*q//8;n=m+6;period=11*m*m+105*m;T=1<<n
    c=M.initialize(q,(0,)*m);initial=c.snapshot();a=[];b=[None]*period
    def visible():
        z=np.zeros((2,m,2),dtype=np.int64)
        cursor=next((j-3 for j,t in enumerate(c.tape) if t[0]),None)
        for side,port,color in ((0,0,1),(1,n-1,n-2)):
            if c.tape[port][2] and c.pos!=port:
                assert cursor is not None
                z[side,cursor,0]=1
                z[side,cursor,1]=int(c.pos!=color)
        return z.reshape(-1)
    b[0]=visible()
    for k in range(1,period+1):
        val=visible()  # begin clock vertex
        for j in range(n-1):M.bond(c,j,q);val+=visible()
        val+=visible()  # bounce clock vertex
        for j in range(n-1,-1,-1):M.onsite(c,j,q);val+=visible()
        a.append(val);b[k%period]=visible()
    assert c.snapshot()==initial
    # Sweep k uses a[(k-1)%period]. Counter dwell after it uses b[k%period].
    wa=[T//period+int(j<T%period) for j in range(period)]
    R=T-1;wb=[R//period]*period
    for k in range(1,R%period+1):wb[k%period]+=1
    for t in range(1,n):
        stride=pow(2,t,period);p=period//gcd(1<<t,period);steps=(1<<(n-t))-1
        whole,rem=divmod(steps,p)
        for j in range(1,p+1):wb[(stride*j)%period]+=2*(whole+int(j<=rem))
    wb[T%period]+=n-1
    result=(np.array(wa,dtype=object)@np.array(a,dtype=object)+np.array(wb,dtype=object)@np.array(b,dtype=object)).reshape(2,m,2)
    C=(2*n+4)*T+n-4
    # Count all identity clock intervals as an independent arithmetic check.
    assert sum(wa)*(2*n+1)+sum(wb)+2*n==C
    assert min(result.flat)>0
    return dict(q=q,m=m,n=n,period=period,clock_vertices=C,
                port_counts=result[:,:,0].tolist(),density_counts=result[:,:,1].tolist())

if __name__=='__main__':
    out=[counts(q) for q in (4,8,12)]
    with open('clean_readout_checks.json') as f:old=json.load(f)['clock']
    for x,y in zip(out,old):
        for k in ('q','m','n','clock_vertices','port_counts','density_counts'):assert x[k]==y[k]
    with open('clean_readout_counts.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps([dict(q=x['q'],m=x['m'],n=x['n'],period=x['period'],clock_vertices=x['clock_vertices'],port_min=min(sum(x['port_counts'],[])),density_min=min(sum(x['density_counts'],[]))) for x in out],indent=2))
