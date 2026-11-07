"""C125 bounded local checks. All output remains in this directory.
Use --cluster for larger finite-size checks, only within SLURM.
"""
from pathlib import Path
from collections import defaultdict
import argparse, json, os, time
import numpy as np
from scipy.linalg import expm, eigvalsh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import logsumexp, digamma
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent
sx=np.array([[0.,1.],[1.,0.]])/2
sy=np.array([[0.,-1j],[1j,0.]])/2
sz=np.diag([.5,-.5])

def bond(delta):return (np.kron(sx,sx)+np.kron(sy,sy)+delta*np.kron(sz,sz)).real

def ring_spectrum(d,L,delta=0,ground=False):
    powers=d**np.arange(L,dtype=np.int64);sectors=defaultdict(list)
    for state in range(d**L):
        digits=(state//powers)%d
        sectors[tuple(np.count_nonzero(digits==a) for a in range(d-1))].append(state)
    result=[]
    for basis in sectors.values():
        ind={s:i for i,s in enumerate(basis)};rows=[];cols=[];vals=[]
        for col,state in enumerate(basis):
            digits=(state//powers)%d;diag=0.
            for i in range(L):
                j=(i+1)%L;a,b=int(digits[i]),int(digits[j])
                if d==3:
                    if a==b:diag+=1
                    else:
                        target=state+(b-a)*powers[i]+(a-b)*powers[j]
                        rows.append(ind[target]);cols.append(col);vals.append(1.)
                else:
                    diag+=delta*(.25 if a==b else -.25)
                    if a!=b:
                        target=state+(b-a)*powers[i]+(a-b)*powers[j]
                        rows.append(ind[target]);cols.append(col);vals.append(.5)
            rows.append(col);cols.append(col);vals.append(diag)
        H=coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis))).tocsr()
        if ground and len(basis)>60:
            ev=eigsh(H,k=1,which='SA',return_eigenvectors=False,tol=1e-11)
        else:ev=eigvalsh(H.toarray(),check_finite=False)
        result.extend([float(min(ev))] if ground else ev.tolist())
    return np.sort(result)

def xx_logz(L,beta):
    answers=[]
    for shift,parity in [(.5,0),(0.,1)]:
        en=np.cos(2*np.longdouble(np.pi)*(np.arange(L,dtype=np.longdouble)+shift)/L)
        even,odd=np.longdouble(0),np.longdouble(-np.inf)
        for e in en:even,odd=np.logaddexp(even,odd-beta*e),np.logaddexp(odd,even-beta*e)
        answers.append((even,odd)[parity])
    return np.logaddexp(*answers)

def transfer(d,h,beta,m):
    g=expm(-beta*h/m).reshape(d,d,d,d)
    paths=np.array(np.unravel_index(np.arange(d**(2*m)),(d,)*(2*m))).T
    mats=[]
    for parity in [0,1]:
        K=np.ones((len(paths),len(paths)),complex)
        for j in range(parity,2*m,2):
            jj=(j+1)%(2*m)
            K*=g[paths[:,jj,None],paths[None,:,jj],paths[:,j,None],paths[None,:,j]]
        mats.append(K)
    return mats[0]@mats[1]

def xx_cross_logz(L,beta,nmax=20000):
    # Exact infinite-product spectrum, evaluated with a stated cutoff.
    odd=np.arcsinh(np.pi*(2*np.arange(nmax)+1)/beta)
    even=np.arcsinh(2*np.pi*(np.arange(nmax)+1)/beta)
    plus=4*np.log1p(np.exp(-L*odd)).sum()
    minus=4*np.log1p(-np.exp(-L*odd)).sum()
    # Integrable logarithmic endpoint in the second integral.
    if beta<80:
        pc=2/np.pi*quad(lambda k:np.logaddexp(beta*np.cos(k)/2,-beta*np.cos(k)/2),0,np.pi/2,epsabs=1e-12)[0]
        gap=2/np.pi*quad(lambda k:-np.log(np.tanh(beta*np.cos(k)/2)),0,np.pi/2,epsabs=1e-12)[0]
    else:
        pc=beta/np.pi+2/(np.pi*beta)*quad(lambda y:np.log1p(np.exp(-y))/np.sqrt(1-(y/beta)**2),0,40,epsabs=1e-12)[0]
        gap=2/(np.pi*beta)*quad(lambda y:-np.log(np.tanh(y/2))/np.sqrt(1-(y/beta)**2),0,40,epsabs=1e-12)[0]
    evensector=np.logaddexp(plus,minus)-np.log(2)
    other=np.log(2)-L*gap+4*np.log1p(np.exp(-L*even)).sum()
    return L*pc+np.logaddexp(evensector,other),gap

def curvature(L,E,e0):
    q=L*L;theta=-L*E;slopes=np.diff(theta)/np.diff(q)
    ca=6/np.pi*(e0*q+theta)
    return dict(L=L.tolist(),E=E.tolist(),cA=ca.tolist(),max_cA_increase=float(max(np.diff(ca))),
                min_theta_slope_increase=float(min(np.diff(slopes))))

def run():
    parser=argparse.ArgumentParser();parser.add_argument('--cluster',action='store_true');args=parser.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Large checks require SLURM')
    started=time.time();out={};models=[('XX',2,0.),('XXZ_m09',2,-.9),('XXZ_m05',2,-.5),('XXZ_p05',2,.5),('Heisenberg',2,1.),('XXZ_p2',2,2.),('SU3',3,0.)]
    spectra={};ground_data={};tracechecks=[]
    for name,d,delta in models:
        lengths=range(2,9) if d==2 else range(2,9)
        if args.cluster:lengths=range(2,13) if d==2 else range(2,10)
        spectra[name]={L:ring_spectrum(d,L,delta) for L in lengths}
        maxL=12 if d==2 else 9
        grounds={L:float(ev[0]) for L,ev in spectra[name].items()}
        for L in range(2,maxL+1):
            if L not in grounds and (L%2==0 if d==2 else L%3==0):grounds[L]=float(ring_spectrum(d,L,delta,True)[0])
        ground_data[name]=grounds
        betas=[.05,.2,1.,2.,8.,32.,128.]
        for beta in betas:
            logs={L:float(logsumexp(-beta*ev)) for L,ev in spectra[name].items()}
            rs={L:2*logs[L]-logs[2*L] for L in logs if 2*L in logs}
            for L,r in rs.items():
                M=1. if d==3 else (max(eigvalsh(bond(delta)))-min(eigvalsh(bond(delta))))/2
                tracechecks.append(dict(model=name,L=L,beta=beta,R=r,
                    trivial_margin=float(2*L*np.log(d)+4*beta*L*M-r),
                    gluing_margin=float(4*beta*M-abs(r))))
        if name=='XX':
            out['xx_ED_max_error']=max(abs(float(xx_logz(L,b))-float(logsumexp(-b*ev))) for L,ev in spectra[name].items() if L%2==0 for b in betas)
    out['grounds']=ground_data;out['finite_trace_checks']=tracechecks
    out['doubling_ground']={name:{str(L):grounds[2*L]-2*e for L,e in grounds.items() if 2*L in grounds} for name,grounds in ground_data.items()}
    out['R_even_monotonicity']={}
    for name,_,_ in models:
        rows=[r for r in tracechecks if r['model']==name and r['L']%2==0]
        diffs=[]
        for beta in sorted(set(r['beta'] for r in rows)):
            v=sorted((r['L'],r['R']) for r in rows if r['beta']==beta)
            diffs.extend(v[i+1][1]-v[i][1] for i in range(len(v)-1))
        out['R_even_monotonicity'][name]=max(diffs) if diffs else None
    # Realignment positivity and the complete small-Trotter spectrum.
    qt=[]
    for name,d,delta in models:
        h=bond(delta) if d==2 else np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
        for beta in [.2,2.,8.]:
            for m in ([1,2,3,4] if d==2 else [1,2,3]):
                ev=np.linalg.eigvals(transfer(d,h,beta,m));scale=max(abs(ev));en=ev/scale
                z4=np.sum(ev**2);exact=logsumexp(-beta*spectra[name][4])
                qt.append(dict(model=name,beta=beta,m=m,min_real_relative=float(min(en.real)),
                    max_imag_relative=float(max(abs(en.imag))),negative_count=int(sum(en.real<-1e-10)),
                    complex_count=int(sum(abs(en.imag)>1e-10)),logZ4_error=float(np.log(z4.real)-exact)))
    out['qtm']=qt
    # cA and theta on exact rings. XXZ e0 by the zero-field Bethe integral.
    cdata={}
    for name,d,delta in models:
        if d==3:e0=1+2*(digamma(1/3)+np.euler_gamma)/3;Ls=np.array([3,6,9.])
        else:
            Ls=np.array([4,6,8,10,12.])
            if delta==1:e0=.25-np.log(2)
            elif abs(delta)<1:
                gamma=np.arccos(delta)
                def f(x):
                    if x==0:return (np.pi-gamma)/np.pi
                    # stable form of sinh[(pi-gamma)x]/(sinh(pi x) cosh(gamma x))
                    return 2*np.exp(-2*gamma*x)*(-np.expm1(-2*(np.pi-gamma)*x))/((-np.expm1(-2*np.pi*x))*(1+np.exp(-2*gamma*x)))
                e0=delta/4-np.sin(gamma)*quad(f,0,np.inf,epsabs=1e-12)[0]
            else:
                # Theta convexity is independent of e0; omit cA for massive sample.
                e0=0.
        E=np.array([ground_data[name][int(L)] for L in Ls]);cdata[name]=curvature(Ls,E,e0)
        if delta>1:cdata[name].pop('cA');cdata[name].pop('max_cA_increase')
    out['ring_curvatures']=cdata
    # Requested existing Bethe data, read only its N,k,L,E fields for this test.
    source=ROOT/'../../contrib/casimir_survey/astra_bethe/bethe_exact_energies.json'
    raw=json.loads(source.read_text());bethe=[]
    for fam in raw['families']:
        N,k=fam['N'],fam['k']
        if k!=2 or N not in [2,3,4,5,7,8]:continue
        entries=[(r['L'],r['E']) for r in fam['data']]
        if not any(L==N for L,E in entries):entries=[(N,0.)]+entries
        L=np.array([r[0] for r in entries],float);E=np.array([r[1] for r in entries])
        Hk=sum(1/j for j in range(1,k+1));e0=Hk-sum(digamma(1+(j-1)/N)-digamma(j/N) for j in range(1,k+1))/N
        r=curvature(L,E,e0);r.update(N=N,k=k,M=Hk/2,A_over_M=k*(N*N-1)/(N+k)*np.pi/N/(Hk/2),cA_first_over_M=r['cA'][0]/(Hk/2));bethe.append(r)
    out['bethe_sym2']=bethe
    # Large free-fermion rings and the crossed spectrum.
    big=[]
    for L in [4,6,16,64,256,1024,4096,8192]:
        for aspect in [.01,.1,1.,4.,32.,128.]:
            beta=aspect*L;r=2*xx_logz(L,beta)-xx_logz(2*L,beta)
            big.append(dict(L=L,beta=beta,R=float(r),bound_margin=float(2*np.log(2)+5*np.pi/3*beta/L-r)))
    out['xx_large']=big
    cross=[]
    for L,beta in [(4,.2),(4,2),(6,8),(16,64),(64,256),(256,4096),(1024,32768)]:
        z,g=xx_cross_logz(L,beta);z2,_=xx_cross_logz(L,beta,40000)
        cross.append(dict(L=L,beta=beta,error=float(z-xx_logz(L,beta)),cutoff_change=float(z2-z),gap=g,scaled_gap=beta*g))
    out['xx_crossed_products']=cross
    out['matsubara']=[dict(beta=b,k=k,omega=2*np.pi*k/b,rate=float(np.arcsinh(np.pi/b)+np.arcsinh((2*k-1)*np.pi/b)),linear_claim=2*np.pi*k/b/(2*np.e)) for b in [1.,10.,100.] for k in [1,2,10,100,1000]]
    # A compact-boson level counting test for the finite-rank countermodel.
    # Levels j=2 Delta. One boson: two integer zero modes and two oscillator partitions.
    jmax=500
    p=np.zeros(jmax+1,dtype=object);p[0]=1
    for k in range(1,jmax//2+1):
        for j in range(2*k,jmax+1):p[j]+=p[j-2*k]
    osc=np.convolve(p,p)[:jmax+1];zero=np.zeros(jmax+1,dtype=object)
    for m in range(-int(np.sqrt(jmax)),int(np.sqrt(jmax))+1):
        for w in range(-int(np.sqrt(jmax)),int(np.sqrt(jmax))+1):
            if m*m+w*w<=jmax:zero[m*m+w*w]+=1
    one=np.convolve(osc,zero)[:jmax+1];coef=np.zeros(jmax+1,dtype=object);coef[0]=1;counter=[]
    for n in range(1,9):
        coef=np.convolve(coef,one)[:jmax+1]
        if n not in [1,2,4,8]:continue
        for beta in [2.,8.,16.,32.]:
            # Complete energy shells only, rank no larger than floor exp(beta/2).
            cum=0;last=0
            for j,c in enumerate(coef):
                if np.log(float(cum+int(c)))>beta/2:break
                cum+=int(c);last=j
            logmult=np.array([np.log(float(c)) if c else -np.inf for c in coef[:last+1]])
            energies=np.pi*np.arange(last+1)/beta
            for L in [1.,2.,8.,32.]:
                g1=logsumexp(logmult-L*energies);g2=logsumexp(logmult-2*L*energies);r=2*g1-g2
                counter.append(dict(n=n,beta=beta,L=L,cutoff_level=last,enumeration_exhausted=last==jmax,log_rank=float(np.log(float(cum))),R=float(r),LR_over_beta=float(L*r/beta),rank_margin=float(np.log(float(cum))-r)))
    out['countermodel_finite']=counter
    # Analytic thermodynamic limit, B=pi n v/6 and b=1/2.
    asym=[]
    for n in [1,4,16,64,256,1024]:
        B=np.pi*n/6;b=.5;ell=2*B/b
        def g(L):return b-L*b*b/(4*B) if L<ell else B/L
        for f in [.1,.5,1.,2.,10.]:
            L=f*ell;asym.append(dict(n=n,L=L,R_over_beta=2*g(L)-g(2*L),LR_over_beta=L*(2*g(L)-g(2*L)),A=n))
    out['countermodel_limit']=asym
    out['runtime_seconds']=time.time()-started
    (ROOT/('checks_cluster.json' if args.cluster else 'checks.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(runtime=out['runtime_seconds'],xx_ED_error=out['xx_ED_max_error'],
        xx_crossed_max_error=max(abs(r['error']) for r in cross),
        curvatures={k:{a:v for a,v in x.items() if a.startswith('min_') or a.startswith('max_')} for k,x in cdata.items()},
        bethe=[{k:r[k] for k in ['N','max_cA_increase','min_theta_slope_increase','cA_first_over_M','A_over_M']} for r in bethe]),indent=2))
if __name__=='__main__':run()
