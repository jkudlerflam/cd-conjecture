"""C128 small local checks; larger transfer matrices require --cluster under SLURM."""
from pathlib import Path
import os, json, argparse, time
import numpy as np
from scipy.linalg import expm, eigvalsh
from scipy.special import logsumexp
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent

def bond(delta):
    return np.array([[-delta/4,0,0,-.5],[0,delta/4,0,0],
                     [0,0,delta/4,0],[-.5,0,0,-delta/4]])

def layers(L,delta):
    h=bond(delta); hs=[np.zeros((2**L,2**L)) for _ in range(2)]
    for s in range(2**L):
        for j in range(L):
            k=(j+1)%L; a=(s>>j)&1; b=(s>>k)&1
            for aa in range(2):
                for bb in range(2):
                    t=s+(aa-a)*2**j+(bb-b)*2**k
                    hs[j%2][t,s]+=h[2*aa+bb,2*a+b]
    return hs

def transfer(delta,beta,m):
    g=expm(-beta*bond(delta)/m).reshape(2,2,2,2)
    paths=((np.arange(2**(2*m))[:,None]>>np.arange(2*m))&1)
    ks=[]
    for parity in range(2):
        K=np.ones((len(paths),len(paths)))
        for j in range(parity,2*m,2):
            jj=(j+1)%(2*m)
            K*=g[paths[:,jj,None],paths[None,:,jj],paths[:,j,None],paths[None,:,j]]
        ks.append(K)
    w,v=np.linalg.eigh(ks[0]); root=(v*np.sqrt(np.maximum(w,0)))@v.T
    U=root@ks[1]@root
    ev=eigvalsh(U); scale=max(ev); vals=np.maximum(ev/scale,0)
    return ks,ev,vals

def xx_logz(L,beta):
    out=[]
    for shift,parity in [(.5,0),(0.,1)]:
        en=np.cos(2*np.longdouble(np.pi)*(np.arange(L,dtype=np.longdouble)+shift)/L)
        even,odd=np.longdouble(0),np.longdouble(-np.inf)
        for e in en:even,odd=np.logaddexp(even,odd-beta*e),np.logaddexp(odd,even-beta*e)
        out.append((even,odd)[parity])
    return float(np.logaddexp(*out))

def xx_cross(L,beta,cut):
    n=np.arange(cut)
    o=np.arcsinh(np.pi*(2*n+1)/beta); e=np.arcsinh(2*np.pi*(n+1)/beta)
    plus=4*np.log1p(np.exp(-L*o)).sum();minus=4*np.log1p(-np.exp(-L*o)).sum()
    # Transform y=beta*cos(k), resolving the low-temperature endpoint.
    if beta<80:
        pc=2/np.pi*quad(lambda k:np.logaddexp(beta*np.cos(k)/2,-beta*np.cos(k)/2),0,np.pi/2,epsabs=1e-12)[0]
        gap=2/np.pi*quad(lambda k:-np.log(np.tanh(beta*np.cos(k)/2)),0,np.pi/2,epsabs=1e-12)[0]
    else:
        pc=beta/np.pi+2/(np.pi*beta)*quad(lambda y:np.log1p(np.exp(-y))/np.sqrt(1-(y/beta)**2),0,40,epsabs=1e-12)[0]
        gap=2/(np.pi*beta)*quad(lambda y:-np.log(np.tanh(y/2))/np.sqrt(1-(y/beta)**2),0,40,epsabs=1e-12)[0]
    ztilde=np.logaddexp(np.logaddexp(plus,minus)-np.log(2),np.log(2)-L*gap+4*np.log1p(np.exp(-L*e)).sum())
    # Upper bound on omitted logarithmic products; applies to L>1.
    tail=np.exp(np.log(4)+L*np.log(beta/(2*np.pi))+(1-L)*np.log(cut)-np.log(L-1))
    return dict(logtilde=float(ztilde),logZ=float(L*pc+ztilde),plus=float(plus),tail_bound=float(tail))

def bose_linear(x):
    if x==0:return 0.
    k=np.arange(1,max(2,int(np.ceil(40*x)))+1)
    return float(-2*np.log1p(-np.exp(-k/x)).sum())

def run():
    ap=argparse.ArgumentParser();ap.add_argument('--cluster',action='store_true');args=ap.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Large computation requires SLURM')
    start=time.time();out={};rows=[]
    for delta in [0.,.5,1.,2.]:
        M=(1+max(1.,delta))/4
        hs={L:layers(L,delta) for L in [2,4,6]}
        specs={L:eigvalsh(sum(x)) for L,x in hs.items()}
        for beta in [.2,2.,6.,8.]:
            for m in ([1,2,3,4,5] if args.cluster else [1,2,3,4]):
                ks,ev,vals=transfer(delta,beta,m)
                for L in [2,4,6]:
                    temporal=np.trace(np.linalg.matrix_power(expm(-beta*hs[L][0]/m)@expm(-beta*hs[L][1]/m),m))
                    spatial=np.trace(np.linalg.matrix_power(ks[0]@ks[1],L//2))
                    g=np.log(np.sum(vals**(L/2)));g2=np.log(np.sum(vals**L));R=2*g-g2
                    linear=np.log(3)+2*(-np.log1p(-np.exp(-1)))+bose_linear(2*M*beta/L)
                    rows.append(dict(delta=delta,beta=beta,m=m,L=L,relative_contraction_error=float(abs(spatial-temporal)/temporal),
                        positive_trace_error=float(abs(np.sum(ev**(L//2))-spatial)/spatial),min_K_eigenvalue=float(min(eigvalsh(ks[0]))),
                        logZ_Trotter_error=float(np.log(spatial)-logsumexp(-beta*specs[L])),logtilde=float(g),R=float(R),
                        sandwich_lower_margin=float(R-g),sandwich_upper_margin=float(2*g-R),
                        linear_Fock_margin=float(linear-g),XX_special_margin=float(np.log(3)+beta*(1/np.sin(np.pi/L)-L/np.pi)-g)))
    out['transfer']=rows
    count_rows=[];gate_rows=[]
    for delta in [0.,.5,1.,2.]:
        H=sum(layers(4,delta)); W=np.zeros_like(H)
        for state in range(16):
            for j,k in [(0,1),(3,0)]:
                if ((state>>j)&1)==((state>>k)&1):
                    W[state^(1<<j)^(1<<k),state]+=.5
        for beta in [.2,2.,8.]:
            base=logsumexp(-beta*eigvalsh(H))
            for fugacity in [1.25,2.,4.]:
                logmgf=logsumexp(-beta*eigvalsh(H-(fugacity-1)*W))-base
                count_rows.append(dict(delta=delta,beta=beta,fugacity=fugacity,logMGF=float(logmgf),
                    margin=float(beta*(fugacity-1)-logmgf)))
        for dt in [.001,.1,1.,4.]:
            a=np.exp(dt*delta/4)*np.cosh(dt/2);b=np.exp(-dt*delta/4);c=np.exp(dt*delta/4)*np.sinh(dt/2)
            explicit=np.array([[a,0,0,c],[0,b,0,0],[0,0,b,0],[c,0,0,a]])
            gate_rows.append(dict(delta=delta,dt=dt,error=float(np.max(abs(explicit-expm(-dt*bond(delta))))),
                                  realignment_min=float(min(a-b,c)),free_fermion_residual=float(a*a-b*b-c*c)))
    out['count_MGF']=count_rows;out['gates']=gate_rows
    assert min(r['margin'] for r in count_rows)>=-1e-10
    assert max(r['error'] for r in gate_rows)<1e-10
    exact=[]
    for L in [2,4,6,16,64,256,1024]:
        for aspect in [.01,.1,1.,4.,32.]:
            beta=L*aspect;z=xx_cross(L,beta,20000);zz=xx_cross(L,beta,40000)
            R=2*xx_logz(L,beta)-xx_logz(2*L,beta)
            B=np.log(3)+beta*(1/np.sin(np.pi/L)-L/np.pi)
            exact.append(dict(L=L,beta=beta,R=R,logtilde=z['logtilde'],tail_bound=z['tail_bound'],
                crossed_logZ_error=z['logZ']-xx_logz(L,beta),cutoff_change=zz['logZ']-z['logZ'],
                exact_bound_margin=B-z['logtilde'],product_margin=np.log(3)+z['plus']-z['logtilde'],
                midpoint_margin=beta*(1/np.sin(np.pi/L)-L/np.pi)-z['plus'],ratio_bound_margin=2*B-R,
                C119_margin=3*np.log(2)+(5*np.pi/3)*beta/L-R,
                linear_Fock_margin=np.log(3)+2*(-np.log1p(-np.exp(-1)))+bose_linear(beta/L)-z['logtilde']))
    out['xx_exact']=exact
    quadrature=[]
    for L in [2.,3.,4.,6.,16.,64.]:
        I=quad(lambda v:np.cosh(v/L)*np.log1p(np.exp(-v))/L,0,100,epsabs=1e-12)[0]
        exactI=np.pi/(2*np.sin(np.pi/L))-L/2
        quadrature.append(dict(L=L,I=I,identity_error=I-exactI,elementary_margin=L/(L*L-1)-I))
    out['quadrature']=quadrature
    poisson=[]
    for delta in [0.,.5,1.,2.]:
        for dt in [.1,.01,.001,.0001]:
            g=expm(-dt*bond(delta));pair=g[3,0]/g[0,0]
            poisson.append(dict(delta=delta,dt=dt,normalized_pair_weight=float(pair),independent_weight=(dt/2)**2,
                                domination_ratio=float(pair/(dt/2)**2)))
    out['paired_events']=poisson
    out['bridge']=[dict(remaining=r,rate=float(.5/np.tanh(.5*r))) for r in [1.,.1,.01,.001]]
    out['frequency']=[dict(beta=b,k=k,omega=2*np.pi*k/b,
        gap=float(np.arcsinh(np.pi/b)+np.arcsinh((2*k-1)*np.pi/b)),
        linear_rate=2*np.pi*k/b/(2*np.e),
        nonlinear_rate=2*np.pi*k/b/(2*np.e) if 2*np.pi*k/b<=2*np.e else np.log(np.pi*k/b))
        for b in [1.,10.,100.] for k in [1,2,10,100,1000]]
    # Test the conditional nonlinear bosonic sum estimate, independent of a lattice claim.
    nonlinear=[];V=2*np.e
    for L in [2,4,8,16]:
        for beta in [.1,1.,10.,100.]:
            w=2*np.pi*np.arange(1,20001)/beta
            a=np.where(w<=V,w/V,np.log(w/2))
            summ=float(-np.log1p(-np.exp(-L*a)).sum())
            bound=beta*V/(2*np.pi*L)*(np.pi**2/6+2/np.expm1(2))
            nonlinear.append(dict(L=L,beta=beta,sum=summ,bound=bound,margin=bound-summ))
    out['conditional_nonlinear']=nonlinear
    assert max(r['relative_contraction_error'] for r in rows)<1e-10
    assert max(r['positive_trace_error'] for r in rows)<1e-10
    assert min(r['min_K_eigenvalue'] for r in rows)>-1e-9
    assert min(r['sandwich_lower_margin'] for r in rows)>-1e-10
    for key in ['exact_bound_margin','product_margin','midpoint_margin','ratio_bound_margin','C119_margin','linear_Fock_margin']:
        assert min(r[key] for r in exact)>-1e-9,key
    for r in exact:
        assert abs(r['crossed_logZ_error'])<=r['tail_bound']+1e-7
        assert r['exact_bound_margin']>=r['tail_bound']-1e-9
        assert r['linear_Fock_margin']>=r['tail_bound']-1e-9
    ground=[]
    for L in [2,4,6,16,64,256,1024]:
        D=2/np.sin(np.pi/L)-1/np.sin(np.pi/(2*L))
        C=(8-16/np.pi) if L==2 else (16*np.sqrt(2)-64/np.pi)
        ground.append(dict(L=L,scaled_defect=2*L*D,bound=C,margin=C-2*L*D))
    assert min(r['margin'] for r in ground)>0
    out['xx_ground']=ground
    E4=min(eigvalsh(sum(layers(4,1.))))
    out['Heisenberg_XX_extension_failure']=dict(E4=float(E4),
        true_slope=float(4*(.25-np.log(2))-E4),claimed_slope=float(np.sqrt(2)-4/np.pi))
    assert out['Heisenberg_XX_extension_failure']['true_slope']>out['Heisenberg_XX_extension_failure']['claimed_slope']
    assert max(abs(r['identity_error']) for r in quadrature)<1e-12
    assert min(r['margin'] for r in nonlinear)>0
    assert all(r['gap']>=r['nonlinear_rate']-1e-12 for r in out['frequency'])
    out['runtime_seconds']=time.time()-start
    filename='checks_cluster.json' if args.cluster else 'checks.json'
    (ROOT/filename).write_text(json.dumps(out,indent=2)+'\n')
    summary=dict(runtime=out['runtime_seconds'],transfer_cases=len(rows),xx_cases=len(exact),
        max_contraction_error=max(r['relative_contraction_error'] for r in rows),
        max_positive_trace_error=max(r['positive_trace_error'] for r in rows),
        max_crossed_residual_after_tail=max(max(0,abs(r['crossed_logZ_error'])-r['tail_bound']) for r in exact),
        min_interacting_Fock_margin=min(r['linear_Fock_margin'] for r in rows if r['delta']!=0),
        finite_Trotter_XX_special_failure_count=sum(r['XX_special_margin']<-1e-8 for r in rows),
        pair_witness=poisson[-1],bridge_witness=out['bridge'][-1])
    print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
