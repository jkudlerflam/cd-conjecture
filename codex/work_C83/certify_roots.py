# C83 local copy from C80. Original attribution follows.
# C80 local copy of C75 source, used as a reference or scalar certificate.
# C75 reuse of C33 scalar-certificate implementation. See PROOF.md.
"""Rational endpoint inertia and exact chord margins for all Fermi roots.
Numerical eigenvectors only propose a rational congruence. Every acceptance
is an integer/rational diagonal-dominance or inequality check.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,math,time,json
from fractions import Fraction as Q
from functools import cmp_to_key
from rational_bands import integer_ldl

def dyadic_below(x):
    k=max(1,math.ceil(-math.log2(float(x))))
    while Q(1,1<<k)>x:k+=1
    return k

def main(n=23,width=1e-8):
    start=time.monotonic();a=np.load(f'full_loop_n{n}.npz',allow_pickle=True)
    F,M,EF,EM=[a[k] for k in ['F','M','EF','EM']];bits=int(a['bits']);detC=int(a['detC'])
    b=np.load(f'derivative_n{n}.npz',allow_pickle=True);tb=int(b['bits']);scale=1<<tb
    R,ds=integer_ldl(b['F_numerator']+b['F_numerator'].T)
    powers=[(int(d).bit_length()-1)//2 for d in ds];extra=1<<(2*(bits-tb))
    UF=R@F@R.T;UM=R@M@R.T;UI=R@R.T
    E=abs(R)@(EF+EF.T+EM)@abs(R).T
    rad=max(sum((Q(int(E[i,j]),extra*(1<<(powers[i]+powers[j]+1))) for j in range(n)),Q(0)) for i in range(n))
    GH=UF+UF.T
    EH=abs(R)@(EF+EF.T)@abs(R).T
    herm_margin=min(Q(int(GH[i,i])-int(EH[i,i]),extra*(1<<(2*powers[i]+2)))-sum((Q(abs(int(GH[i,j]))+int(EH[i,j]),extra*(1<<(powers[i]+powers[j]+2))) for j in range(n) if j!=i),Q(0)) for i in range(n))
    assert herm_margin>0
    bands=json.load(open(f'full_loop_bands_n{n}.json'));fp=-bands['field_power'];qbits=46;hbits=80
    def inertia(t):
        p,q=t.numerator,t.denominator;c=q*q-p*p;s=2*p*q;d=q*q+p*p
        Hr=-(UF+UF.T)*c+UI*(d*(1<<(2*bits+1-fp)))
        Hi=-(UF-UF.T)*s+UM*d
        rr=np.empty((n,n),dtype=object);ii=rr.copy()
        for i in range(n):
            for j in range(n):
                den=2*extra*d*(1<<(powers[i]+powers[j]))
                rr[i,j]=(int(Hr[i,j])*(1<<hbits))//den
                ii[i,j]=(int(Hi[i,j])*(1<<hbits))//den
        rr=rr+rr.T;ii=ii-ii.T # denominator 2^(hbits+1)
        H=np.array(rr,dtype=float)/(2.**(hbits+1))+1j*np.array(ii,dtype=float)/(2.**(hbits+1))
        _,U=np.linalg.eigh(H)
        qr=np.array(np.rint(U.real*(2.**qbits)),dtype=object).astype(object)
        qi=np.array(np.rint(U.imag*(2.**qbits)),dtype=object).astype(object)
        qr=np.array([[int(x) for x in row] for row in qr],dtype=object)
        qi=np.array([[int(x) for x in row] for row in qi],dtype=object)
        gr=qr.T@qr+qi.T@qi;gi=qr.T@qi-qi.T@qr
        assert min(int(gr[i,i])-sum(abs(int(gr[i,j]))+abs(int(gi[i,j])) for j in range(n) if i!=j) for i in range(n))>0
        assert max(int(gr[i,i]) for i in range(n))<=2*(1<<(2*qbits))
        xr=rr@qr-ii@qi;xi=rr@qi+ii@qr
        zr=qr.T@xr+qi.T@xi;zi=qr.T@xi-qi.T@xr
        assert np.array_equal(zr,zr.T) and np.array_equal(zi,-zi.T)
        gm=min(abs(int(zr[i,i]))-sum(abs(int(zr[i,j]))+abs(int(zi[i,j])) for j in range(n) if i!=j) for i in range(n))
        gap=Q(gm,1<<(hbits+1+2*qbits))-2*n*(rad+Q(2*n,1<<hbits))
        assert gap>0,('inertia enclosure failed',float(t),float(gap),float(rad))
        return sum(int(zr[i,i])<0 for i in range(n)),gap
    intervals=[];smallest=None;counts=[]
    for idx,k in enumerate(bands['roots']):
        den=1<<48
        lo=Q(math.floor(math.tan((k-width)/2)*den),den)
        hi=Q(math.ceil(math.tan((k+width)/2)*den),den)
        il,gl=inertia(lo);ir,gr=inertia(hi)
        assert abs(il-ir)==1,(idx,il,ir)
        counts.append([il,ir]);intervals.append((lo,hi))
        smallest=min([gl,gr]+([] if smallest is None else [smallest]))
        if idx%10==0:print(json.dumps(dict(stage='inertia',completed=idx+1,seconds=time.monotonic()-start)),flush=True)
        if time.monotonic()-start>170:raise TimeoutError('root certification time cap')
    for i in range(len(intervals)-1):assert intervals[i][1]<intervals[i+1][0]
    def phase(t):
        return (1-t*t)/(1+t*t),2*t/(1+t*t)
    zs=[];errs=[]
    for lo,hi in intervals:
        t=(lo+hi)/2;zs.append(phase(t));errs.append(hi-lo) # 2* half-width, atan Lipschitz
    def cmp(a,b):
        ra,ia=a[0];rb,ib=b[0]
        ha=0 if ia>0 or (ia==0 and ra>=0) else 1
        hb=0 if ib>0 or (ib==0 and rb>=0) else 1
        if ha!=hb:return -1 if ha<hb else 1
        cr=ra*ib-ia*rb
        return -1 if cr>0 else (1 if cr<0 else 0)
    def margin(items):
        items=sorted(items,key=cmp_to_key(cmp));lower=None;witness=None
        for i,u in enumerate(items):
            v=items[(i+1)%len(items)]
            ru,iu=u[0];rv,iv=v[0]
            assert ru*iv-iu*rv>=0 # each adjacent oriented gap is at most pi
            d=max(abs(ru-rv),abs(iu-iv))-u[1]-v[1]
            assert d>0,('pair interval overlap',u[2],v[2],float(d))
            if lower is None or d<lower:lower=d;witness=[u[2],v[2]]
        return lower,witness
    d2,_=margin([(zs[i],errs[i],[i]) for i in range(2*n)])
    items=[]
    for i in range(2*n):
        for j in range(i,2*n):
            ri,ii=zs[i];rj,ij=zs[j]
            items.append(((ri*rj-ii*ij,ri*ij+ii*rj),errs[i]+errs[j],[i,j]))
    d4,witness=margin(items)
    # A further normalized Bloch perturbation with norm <=budget in the
    # dyadic congruence keeps every endpoint inertia certificate strict.
    budget=smallest/(4*n);budget_power=dyadic_below(budget)
    traceW=sum((Q(scale*scale*sum(int(x)**2 for x in R[i,:]),1<<(2*powers[i])) for i in range(n)),Q(0))
    coeff=2*n*detC*detC*traceW
    gpower=max(1,math.ceil((math.log2(float(traceW))+math.log2(2*n)+2*math.log2(detC)+budget_power)/2))
    while coeff/Q(1<<(2*gpower))>Q(1,1<<budget_power):gpower+=1
    result=dict(n=n,status='proved_integer_interval_certificate',root_count=2*n,
        D=int(a['D']),a=4,loop_eta_power=-int(a['eta_power']),field_normalized_power=-fp,
        delta2_lower=float(d2),delta4_lower=float(d4),delta4_witness=witness,
        minimum_endpoint_congruence_margin=float(smallest),whitened_matrix_error_bound=float(rad),
        additional_normalized_perturbation_power=-budget_power,finite_g_power=-gpower,
        hermitian_margin_in_dyadic_metric=float(herm_margin),
        log10_finite_g_hermitian_margin_normalized=math.log10(float((herm_margin-Q(1,1<<budget_power))/traceW)),
        seconds=time.monotonic()-start,inertia_changes=counts)
    print(json.dumps(result),flush=True);json.dump(result,open(f'root_certificate_n{n}.json','w'),indent=2)
    # Exact rational brackets and lower bounds are retained as integer strings.
    exact=dict(intervals=[[[str(x.numerator),str(x.denominator)] for x in pair] for pair in intervals],
        delta2=[str(d2.numerator),str(d2.denominator)],delta4=[str(d4.numerator),str(d4.denominator)],
        perturbation_budget_power=budget_power,g_power=gpower)
    json.dump(exact,open(f'root_certificate_exact_n{n}.json','w'),indent=2)
if __name__=='__main__':main()
