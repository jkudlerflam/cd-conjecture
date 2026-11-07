"""Independent C127 checks. Larger matrices require --cluster inside SLURM."""
from pathlib import Path
import argparse,json,os
import numpy as np
from scipy.linalg import expm,eigh
from scipy.special import logsumexp
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent
SX=np.array([[0,1],[1,0]])/2
SY=np.array([[0,-1j],[1j,0]])/2
SZ=np.diag([.5,-.5])

def bond(delta,reflected=False):
    signs=(-1,1,-1) if reflected else (1,1,1)
    return sum(s*c*np.kron(a,a) for s,c,a in zip(signs,[1,1,delta],[SX,SY,SZ])).real

def boson_log(t,S=1):
    if t<=0:return 0.,0.
    n=max(1,int(np.ceil(40*t)));x=np.arange(1,n+1)/t
    val=-S*np.log(-np.expm1(-x)).sum()
    tail=S*np.exp(-(n+1)/t)/((-np.expm1(-1/t))*(-np.expm1(-(n+1)/t)))
    return float(val),float(tail)

def rotation_bases(paths,m):
    D=len(paths);pows=2**np.arange(2*m-1,-1,-1)
    perm=np.roll(paths,2,axis=1)@pows;unseen=set(range(D));orbits=[]
    while unseen:
        a=min(unseen);orb=[a];b=int(perm[a])
        while b!=a:orb.append(b);b=int(perm[b])
        unseen.difference_update(orb);orbits.append(orb)
    bases=[]
    for k in range(m):
        valid=[o for o in orbits if (k*len(o))%m==0];Q=np.zeros((D,len(valid)),complex)
        for j,o in enumerate(valid):Q[o,j]=np.exp(-2j*np.pi*k*np.arange(len(o))/m)/np.sqrt(len(o))
        bases.append(Q)
    return perm,bases

def trotter(delta,beta,m,lengths):
    g=expm(-beta*bond(delta,True)/m).reshape(2,2,2,2)
    paths=np.array(np.unravel_index(np.arange(2**(2*m)),(2,)*(2*m))).T;mats=[]
    for parity in [0,1]:
        K=np.ones((len(paths),len(paths)))
        for j in range(parity,2*m,2):
            jj=(j+1)%(2*m)
            K*=g[paths[:,jj,None],paths[None,:,jj],paths[:,j,None],paths[None,:,j]]
        mats.append(K)
    Ke,Ko=mats;scale=np.linalg.norm(Ke,2);Ke=Ke/scale;Ko=Ko/scale
    w,W=eigh(Ke);root=(W*np.sqrt(np.maximum(w,0)))@W.T
    U=root@Ko@root;U=(U+U.T)/2;perm,Qs=rotation_bases(paths,m)
    comm=np.linalg.norm(U[np.ix_(perm,perm)]-U)/np.linalg.norm(U)
    blocks=[];allvals=[]
    for k,Q in enumerate(Qs):
        vals,vec=eigh(Q.conj().T@U@Q);blocks.append((k,vals,Q@vec));allvals.extend(vals)
    top=max(allvals);rows=[];sources=[]
    _,_,vects=max(blocks,key=lambda b:max(b[1]));v0=vects[:,-1]
    if beta>=2 and m>=3:
        inv=(W*np.where(w>1e-13,1/np.sqrt(np.maximum(w,1e-300)),0))@W.T
        for k in range(m):
            diag=sum(np.exp(2j*np.pi*k*j/m)*(.5-paths[:,2*j]) for j in range(m))/np.sqrt(m)
            src=inv@(diag*(root@v0));norm=np.vdot(src,src).real
            if norm<1e-16:continue
            groups=[]
            for kk,vals,vec in blocks:
                weights=abs(vec.conj().T@src)**2/norm
                for val,weight in zip(vals,weights):
                    if weight>1e-8 and val>top*1e-13:groups.append((kk,float(-.5*np.log(val/top)),float(weight)))
            distinct=[]
            for kk,gap,weight in groups:
                if not any(kk==q and abs(gap-e)<1e-6 for q,e in distinct):distinct.append((kk,gap))
            sources.append(dict(k=k,rank_input=1,distinct_detected_energies=len(distinct),spectral_weights=groups))
    for k,vals,vec in blocks:
        kk=k if k<=m//2 else k-m
        gaps=[float(-.5*np.log(v/top)) for v in vals if v>top*1e-13]
        rows.append(dict(k=k,signed_k=kk,p=2*np.pi*kk/beta,gaps=gaps,discarded=len(vals)-len(gaps)))
    positive=np.array([v for v in allvals if v>0])/top
    M=np.ptp(eigh(bond(delta),eigvals_only=True))/2;V=4*np.e*M;comparisons=[]
    for L in lengths:
        z=float(logsumexp((L/2)*np.log(positive)));z2=float(logsumexp(L*np.log(positive)));R=2*z-z2
        B,tail=boson_log(V*beta/(2*np.pi*L),6)
        comparisons.append(dict(L=L,logZ=z,R=R,Klog3_boson_margin=np.log(3)+B-z,
            Kstar_boson_margin=np.log(3)+8/(3*np.e)+B-z,K0_boson_margin=B-z,boson_tail=tail,
            ratio_lower_margin=R-z,ratio_upper_margin=2*z-R))
    Hs=[np.zeros((16,16)) for _ in range(2)];h0=bond(delta)
    for state in range(16):
        digs=(state//(2**np.arange(4)))%2
        for site in range(4):
            j=(site+1)%4;a,b=digs[site],digs[j]
            for aa in range(2):
                for bb in range(2):
                    target=state+(aa-a)*2**site+(bb-b)*2**j
                    Hs[site%2][target,state]+=h0[2*aa+bb,2*a+b]
    original=np.trace(np.linalg.matrix_power(expm(-beta*Hs[0]/m)@expm(-beta*Hs[1]/m),m))
    spatial=np.sum(np.array(allvals)**2)*scale**4
    edlog=logsumexp(-beta*eigh(Hs[0]+Hs[1],eigvals_only=True))
    return dict(delta=delta,beta=beta,m=m,M=float(M),dimension=len(paths),rotation_commutator=float(comm),
        minimum_relative_eigenvalue=float(min(allvals)/top),gate_min=float(min(w)),
        trace_contraction_relative_error=float(abs(spatial/original-1)),logZ4_Trotter_error=float(np.log(spatial)-edlog),
        sectors=rows,insertion=sources,comparisons=comparisons)

def xx_logz(L,beta):
    ans=[]
    for shift,parity in [(.5,0),(0.,1)]:
        en=np.cos(2*np.longdouble(np.pi)*(np.arange(L,dtype=np.longdouble)+shift)/L)
        even,odd=np.longdouble(0),np.longdouble(-np.inf)
        for e in en:even,odd=np.logaddexp(even,odd-beta*e),np.logaddexp(odd,even-beta*e)
        ans.append((even,odd)[parity])
    return float(np.logaddexp(*ans))

def xx_cross(L,beta,n=20000):
    odd=np.arcsinh(np.pi*(2*np.arange(n)+1)/beta);even=np.arcsinh(2*np.pi*(np.arange(n)+1)/beta)
    if beta<80:
        pc=2/np.pi*quad(lambda k:np.logaddexp(beta*np.cos(k)/2,-beta*np.cos(k)/2),0,np.pi/2,epsabs=1e-12)[0]
        gap=2/np.pi*quad(lambda k:-np.log(np.tanh(beta*np.cos(k)/2)),0,np.pi/2,epsabs=1e-12)[0]
    else:
        pc=beta/np.pi+2/(np.pi*beta)*quad(lambda y:np.log1p(np.exp(-y))/np.sqrt(1-(y/beta)**2),0,40,epsabs=1e-12)[0]
        gap=2/(np.pi*beta)*quad(lambda y:-np.log(np.tanh(y/2))/np.sqrt(1-(y/beta)**2),0,40,epsabs=1e-12)[0]
    pp=4*np.log1p(np.exp(-L*odd)).sum();pm=4*np.log1p(-np.exp(-L*odd)).sum()
    z=np.logaddexp(np.logaddexp(pp,pm)-np.log(2),np.log(2)-L*gap+4*np.log1p(np.exp(-L*even)).sum())
    logtail=np.log(4/(L-1))+L*np.log(beta/(4*np.pi))+(1-L)*np.log(n-.5)
    return float(z),float(pc),float(gap),float(np.exp(logtail))

def analytic_checks():
    out={};xx=[]
    for L in [4,6,16,64,256,1024]:
        for ratio in [.01,.1,1,4,32]:
            beta=L*ratio;z,pc,g,tail=xx_cross(L,beta);B,_=boson_log(2*np.e*beta/(2*np.pi*L),6)
            xx.append(dict(L=L,beta=beta,logZ=z,tail_bound=tail,ring_log_error=z+L*pc-xx_logz(L,beta),
                integral_bound_margin=np.log(3)+2*beta/np.pi*L/(L*L-1)-z,
                simple_bound_margin=np.log(3)+8*beta/(3*np.pi*L)-z,
                product_bound_margin=np.log(3)+8/(3*np.e)+B-z,K0_margin=B-z,delta=g))
    # L=2 uses the exact physical trace; retain explicit product tail diagnostics.
    for beta in [.2,2,8,32]:
        z,pc,g,tail=xx_cross(2,beta);zexact=xx_logz(2,beta)-2*pc;B,_=boson_log(np.e*beta/(2*np.pi),6)
        xx.append(dict(L=2,beta=beta,logZ=zexact,tail_bound=tail,ring_log_error=z+2*pc-xx_logz(2,beta),
            integral_bound_margin=np.log(3)+4*beta/(3*np.pi)-zexact,
            simple_bound_margin=np.log(3)+4*beta/(3*np.pi)-zexact,
            product_bound_margin=np.log(3)+8/(3*np.e)+B-zexact,K0_margin=B-zexact,delta=g))
    out['XX']=xx
    rank=[]
    for beta in [20,100,500,2000]:
        j=np.arange(10000);nu=(2*j+1)*np.pi/beta;omega=2*np.pi/beta
        gaps=np.arcsinh(nu)+np.arcsinh(abs(omega-nu))
        rank.append(dict(beta=beta,omega=omega,energy_cutoff=.5,distinct_pairs_below_cutoff=int(sum(gaps[1:]<=.5)),
            lowest_gap=float(gaps[0]),low_frequency_bound=omega/(2*np.e)))
    out['XX_fixed_momentum_count']=rank;cone=[]
    for beta in [1,10,100]:
        for k in [1,2,10,100,1000]:
            w=2*np.pi*k/beta;E=np.arcsinh(np.pi/beta)+np.arcsinh((2*k-1)*np.pi/beta)
            alpha=w/(2*np.e) if w<=2*np.e else np.log(w/2)
            cone.append(dict(beta=beta,k=k,omega=w,gap=E,linear_margin=E-w/(2*np.e),corrected_margin=E-alpha))
    out['cone']=cone;product=[]
    for t in np.geomspace(.001,1000,50):
        B,tail=boson_log(t);lower=0 if t<1 else -np.log1p(-np.exp(-1))*t/2
        product.append(dict(t=float(t),logB=B,upper_margin=np.pi**2*t/6-B-tail,lower_margin=B-lower))
    out['product_bounds']=product;theta=[]
    for r in [1,2,10,100]:
        rho=4/r
        for x in [.001,.01,.1,1,10,100]:
            a=rho*x;n=max(1,int(np.ceil(np.sqrt(80*a))))
            lt=r*np.log1p(2*np.exp(-np.arange(1,n+1,dtype=float)**2/(2*a)).sum())
            theta.append(dict(rank=r,x=x,logtheta=float(lt),budget_margin=float(16*x-lt)))
    out['theta']=theta;soft=[]
    for L in [2,4,6,16,64]:
        for beta in [.2,2,8,100]:
            M=.5;V=4*np.e*M;n=20000;w=2*np.pi*np.arange(1,n+1)/beta
            alpha=np.where(w<=V,w/V,np.log(w/(4*M)));B=-np.log(-np.expm1(-L*alpha)).sum()
            c=2*M*beta/np.pi;q=(c/(n+1))**L;tail=c**L*n**(1-L)/(L-1)/(1-q)
            coeff=np.e*np.pi/3+4/(np.pi*np.e*(1-np.exp(-2)))
            soft.append(dict(L=L,beta=beta,logB=float(B),tail=float(tail),margin=float(coeff*M*beta/L-B-tail)))
    out['softened_reference']=soft;decay=[]
    for beta in [2,10,100]:
        for k in [1,2,5]:
            w=2*np.pi*k/beta;V=2*np.e
            if w>=V:continue
            j=np.arange(-30000,30001);nu=(2*j+1)*np.pi/beta
            rates=np.arcsinh(abs(nu))+np.arcsinh(abs(nu-w))
            cs=np.sign(nu*(nu-w))/np.sqrt((1+nu**2)*(1+(nu-w)**2))/beta
            for x in [2,4,8,16]:
                G=abs(np.sum(cs*np.exp(-x*rates)));bd=.5*V/(w*(V-w))*np.exp(-w*x/V)
                decay.append(dict(beta=beta,k=k,x=x,correlator=float(G),bound=float(bd),margin=float(bd-G)))
    out['XX_decay']=decay;potts=[]
    for d in [2,3,5]:
        for beta in [4,8,16,24]:
            L=2*max(1,int(np.exp(beta/2)/2));r=-np.expm1(-beta)/(1+(d-1)*np.exp(-beta))
            z=np.log1p((d-1)*np.exp(L*np.log(r)));B,_=boson_log(2*np.e*beta/(2*np.pi*L),2*(d*d-1))
            potts.append(dict(d=d,beta=beta,L=L,logZ=float(z),logd=float(np.log(d)),K0_margin=float(B-z)))
    out['classical_zero_sector']=potts
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cluster',action='store_true');a=ap.parse_args()
    if a.cluster and not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Large checks require SLURM')
    out=analytic_checks();out['transfer']=[];ms=[2,3,4,5,6] if a.cluster else [1,2,3,4]
    for d in [0,.5,1,2]:
        for b in [.2,2,8]:
            for m in ms:out['transfer'].append(trotter(d,b,m,[4,6,16,64]))
    assert all(abs(r['ring_log_error'])<r['tail_bound']+2e-7 for r in out['XX'])
    assert min(r['integral_bound_margin'] for r in out['XX'])>=-1e-8
    assert min(r['product_bound_margin'] for r in out['XX'])>=-1e-8
    assert min(r['corrected_margin'] for r in out['cone'])>=-1e-12
    for key in ['upper_margin','lower_margin']:assert min(r[key] for r in out['product_bounds'])>=-1e-9
    assert min(r['budget_margin'] for r in out['theta'])>=-1e-12
    for key in ['softened_reference','XX_decay']:assert min(r['margin'] for r in out[key])>=-1e-9
    assert max(r['trace_contraction_relative_error'] for r in out['transfer'])<1e-9
    assert max(r['rotation_commutator'] for r in out['transfer'])<1e-9
    assert min(r['minimum_relative_eigenvalue'] for r in out['transfer'])>-1e-10
    for key in ['Kstar_boson_margin','ratio_lower_margin','ratio_upper_margin']:
        assert min(c[key] for r in out['transfer'] for c in r['comparisons'])>=-1e-9
    (ROOT/('checks_cluster.json' if a.cluster else 'checks.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(XX_error_L4plus=max(abs(r['ring_log_error']) for r in out['XX'] if r['L']>=4),
        transfer_error=max(r['trace_contraction_relative_error'] for r in out['transfer']),
        Klog3_min_margin=min(c['Klog3_boson_margin'] for r in out['transfer'] for c in r['comparisons']),
        Kstar_min_margin=min(c['Kstar_boson_margin'] for r in out['transfer'] for c in r['comparisons']),
        max_insertion_cyclic_dimension=max(s['distinct_detected_energies'] for r in out['transfer'] for s in r['insertion'])),indent=2))
if __name__=='__main__':main()
