"""Exact rational sixth-order lifting of the L=12 dark-state cancellation.

Requires sympy in addition to the ordinary numerical dependencies. Uses
unnormalized translation/reflection orbits, so all matrix entries are rational
complex numbers. The final physical normalization is sqrt(30).
"""
from quartet import Sector
import sympy as sy
import json
import time
from pathlib import Path


def exact_basis(sec):
    cols=[]; seen=set(); weights=[]
    for j,s in enumerate(sec.reps):
        if j in seen: continue
        xs=[x for x in range(sec.L) if (s>>x)&1]
        out=sum(1<<((-x)%sec.L) for x in xs)
        i,_=sec.lookup[out]
        sign=(-1)**(sec.N//2+sum(xs))
        seen.update((i,j))
        if i==j:
            if sign==1:
                cols.append({j:1}); weights.append(int(sec.periods[j]))
        else:
            cols.append({j:1,i:sign}); weights.append(2*int(sec.periods[j]))
    W=sy.zeros(len(sec.reps),len(cols))
    left=sy.zeros(len(cols),len(sec.reps))
    for col,entries in enumerate(cols):
        for row,val in entries.items():
            W[row,col]=val; left[col,row]=sy.Rational(val,len(entries))
    return W,left,weights


def exact_matrix(source,target,transitions):
    mat=sy.zeros(len(target.reps),len(source.reps))
    for col,state in enumerate(source.reps):
        for out,amp in transitions(state):
            row,_=target.lookup[out]
            if row>=0:
                mat[row,col]+=amp*sy.Rational(int(source.periods[col]),int(target.periods[row]))
    return mat


def model(L=12):
    sectors=[Sector(L,n,0,False) for n in range(0,L+1,2)]
    bases=[exact_basis(s) for s in sectors]
    sizes=[x[0].cols for x in bases]
    ends=[0]
    for n in sizes: ends.append(ends[-1]+n)
    D=sy.zeros(ends[-1]); hop=D.copy(); A=D.copy()
    numbers=[]; weights=[]
    for block,(s,(W,left,metric)) in enumerate(zip(sectors,bases)):
        numbers.extend([s.N]*W.cols); weights.extend(metric)
        def adj(state):
            yield state, sum(((state>>j)&1)*((state>>((j+1)%L))&1) for j in range(L))
        def move(state):
            for j in range(L):
                if not (state>>j)&1: continue
                for direction,amp in ((1,sy.I),(-1,-sy.I)):
                    k=(j+direction)%L
                    if not (state>>k)&1:
                        yield state^(1<<j)^(1<<k),amp
        cut=slice(ends[block],ends[block+1])
        D[cut,cut]=left*exact_matrix(s,s,adj)*W
        hop[cut,cut]=left*exact_matrix(s,s,move)*W
        if block+1<len(sectors):
            nxt=sectors[block+1]; Wn,Ln,_=bases[block+1]
            def create(state):
                for j in range(L):
                    bits=(1<<j)|(1<<((j+1)%L))
                    if not state&bits: yield state|bits,1
            def destroy(state):
                for j in range(L):
                    bits=(1<<j)|(1<<((j+1)%L))
                    if state&bits==bits: yield state^bits,1
            cutn=slice(ends[block+1],ends[block+2])
            A[cutn,cut]=Ln*exact_matrix(s,nxt,create)*W
            A[cut,cutn]=left*exact_matrix(nxt,s,destroy)*Wn
    G=sy.diag(*weights); N=sy.diag(*numbers)
    for X in (D,hop,A): assert G*X == X.conjugate().T*G
    return D,hop,A,N,G,numbers


def calculate():
    D,hop,A,N,G,nums=model()
    B=hop/2-2*N
    p=[i for i in range(D.rows) if D[i,i]==0]
    q=[i for i in range(D.rows) if D[i,i]!=0]
    pp=lambda X:X.extract(p,p)
    qq=lambda X:X.extract(q,q)
    Apq=A.extract(p,q); Aqp=A.extract(q,p)
    Bpq=B.extract(p,q); Bqp=B.extract(q,p)
    Di=sy.diag(*[1/D[i,i] for i in q]); Aq=qq(A); Bq=qq(B)
    e0=sy.Integer(-12)
    def coefficients(mu2=sy.Integer(0),e2=sy.Integer(0)):
        C=-mu2*N/2
        X={1:Aq,2:Bq-e0*sy.eye(len(q)),4:qq(C)-e2*sy.eye(len(q))}
        R=[-Di]
        for order in range(1,5):
            R.append(-Di*sum((X[j]*R[order-j] for j in X if j<=order),sy.zeros(len(q))))
        F={}
        for order in range(2,7):
            F[order]=Apq*R[order-2]*Aqp
            if order>=3: F[order]+=Apq*R[order-3]*Bqp+Bpq*R[order-3]*Aqp
            if order>=4: F[order]+=Bpq*R[order-4]*Bqp
        F[2]+=pp(B); F[4]+=pp(C)
        return F
    F=coefficients(); Gp=pp(G)
    vacuum=sy.zeros(len(p),1); vacuum[[nums[i] for i in p].index(0)]=1
    n4=[j for j,i in enumerate(p) if nums[i]==4]
    dark4=(F[2].extract(n4,n4)-e0*sy.eye(len(n4))).nullspace()
    assert len(dark4)==1
    dark=sy.zeros(len(p),1)
    for j,coef in zip(n4,dark4[0]): dark[j]=coef
    # Fix the uniform orbit amplitude to one, with no guessed normalization.
    dark/=next(x for x in dark if x!=0)
    T=vacuum.row_join(dark)
    metric=T.conjugate().T*Gp*T
    assert metric==sy.diag(1,30)
    left=metric.inv()*T.conjugate().T*Gp
    Pt=T*left
    R=(e0*sy.eye(len(p))-F[2]+Pt).inv()-Pt
    assert R*T==sy.zeros(len(p),2)
    M2=left*(F[4]+F[3]*R*F[3])*T
    assert M2[1,0]==M2[0,1]==0
    mu2=sy.simplify((M2[1,1]-M2[0,0])/2); e2=M2[0,0]
    F=coefficients(mu2,e2)
    X,Y,Z,W=[F[j] for j in (3,4,5,6)]
    M4=left*(W+Z*R*X+X*R*Z+Y*R*Y
             +Y*R*X*R*X+X*R*Y*R*X+X*R*X*R*Y
             +X*R*X*R*X*R*X-e2*X*R*R*X)*T
    physical=sy.simplify(sy.sqrt(30)*M4[1,0])
    assert sy.simplify(abs(physical)-6*sy.sqrt(30)/5)==0
    return dict(dimension=D.rows, blockaded_dimension=len(p),
                leading_energy=str(e0), energy_epsilon4=str(e2),
                mu_eff_epsilon2=str(mu2),
                leading_off_diagonal=str(M2[1,0]),
                epsilon6_off_diagonal=str(physical),
                gap_over_2_t_eps4=str(2*abs(physical)),
                dark_state_norm_squared=str(metric[1,1]),
                method='Exact rational complex matrices; all symmetry sectors stated in notes.')


if __name__=='__main__':
    start=time.perf_counter(); out=calculate(); out['seconds']=time.perf_counter()-start
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
