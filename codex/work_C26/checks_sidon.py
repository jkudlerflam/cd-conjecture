"""C26 finite-cell checks. Python 3, numpy/scipy; all outputs local."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import itertools, json, time
import numpy as np
import scipy.linalg as la
import sidon_cell as sc

def tl(a):
    return a-np.trace(a)/len(a)*np.eye(len(a))

def norm(a):return float(la.norm(a,2))

def residual(a,n,cells=2):
    gs,_,p=sc.modes(n); N=2**n
    gam=[sc.kron_all([p]*j+[g]+[np.eye(N)]*(cells-j-1)) for j in range(cells) for g in gs]
    rem=tl(a).copy()
    for i in range(len(gam)):
        for j in range(i):
            op=1j*gam[i]@gam[j]
            rem-=np.trace(op@a)/len(a)*op
    return norm(rem)

def margins(n):
    d=sc.sidon_data(n); p=np.array(d['points']); circ=lambda x:abs((x+np.pi)%(2*np.pi)-np.pi)
    d2=min([circ(2*x) for x in p]+[circ(p[i]+s*p[j]) for i,j in itertools.combinations(range(2*n),2) for s in [-1,1]])
    d4=min([circ(sum(s*p[i] for s,i in zip(ss,ids))) for ids in itertools.combinations(range(2*n),4) for ss in itertools.product([-1,1],repeat=4)] or [float('inf')])
    t,th,mu,pts,_=sc.target_parameters(n)
    d.update(delta2=float(d2),delta4=float(d4) if np.isfinite(d4) else None,
             theta=th.tolist(),mu=mu.tolist(),vmin=float(np.min(2*np.sin((pts[:,1]-pts[:,0])/2))),
             root_error=float(np.max(abs(-2*np.cos(pts+th[:,None])+mu[:,None]))))
    return d

def side_embedding(c,side):
    # The unused Bell pair remains Omega. The used pair has four arbitrary states.
    N,k,D,K=c['N'],c['k'],c['D'],c['K']
    E=np.zeros((D,4*(k+1)*N),complex)
    for p in range(2):
      for s in range(2):
       for f in range(k+1):
        for z in range(N):
         col=((2*p+s)*(k+1)+f)*N+z
         for b in range(2):
          pl,sl,pr,sr=(p,s,b,b) if side=='L' else (b,b,p,s)
          row=(2*pl+pr)*K+f*4*N+(2*sl+sr)*N+z
          E[row,col]=1/np.sqrt(2)
    return E

def spectral(c,side=None):
    E=side_embedding(c,side) if side else np.eye(c['D'])
    h=E.conj().T@c['H']@E; g=E.conj().T@c['G']; N=c['N']
    vals,u=la.eigh(h)
    u=np.column_stack([g,u[:,N:]])
    ev=np.r_[np.zeros(N),np.rint(vals[N:])]
    physical=E@u
    ports={s:[physical.conj().T@f@physical for f in c['F'+s]] for s in ['L','R']}
    err=norm(c['H']@physical-physical*ev)
    return ev,ports,err

def pair_engine(c,reduced=True):
    er,fr,err=spectral(c,'R' if reduced else None)
    el,fl,erl=spectral(c,'L' if reduced else None)
    D=len(er); N=c['N']; low=N*N; dim=D*D
    def V(x):
        a=x.reshape(D,D,-1); out=np.zeros_like(a)
        for r,l in zip(fr['R'],fl['L']):
            first=(r@a.reshape(D,-1)).reshape(D,D,-1)
            out+=(l@first.transpose(1,0,2).reshape(D,-1)).reshape(D,D,-1).transpose(1,0,2)
        return out.reshape(dim,-1)
    idx=np.array([i*D+j for i in range(N) for j in range(N)])
    g=np.zeros((dim,low),complex);g[idx,np.arange(low)]=1
    energy=(er[:,None]+el[None,:]).ravel(); inv=np.zeros(dim)
    np.divide(1,energy,out=inv,where=energy>0)
    b=V(g);w=inv[:,None]*b; y=V(w)
    A=b.conj().T@w;T=w.conj().T@w
    k2=-A; k3=w.conj().T@y
    irreduc=-y.conj().T@(inv[:,None]*y)
    folded=(T@A+A@T)/2
    k4=irreduc+folded
    return dict(K2=k2,K3=k3,K4=k4,irreduc=irreduc,folded=folded,V=V,
                g=g,idx=idx,energy=energy,inv=inv,b=b,dimension=dim,spectral_error=max(err,erl))

def exact_graph(e,lam):
    V,g,idx,inv,b=e['V'],e['g'],e['idx'],e['inv'],e['b']
    xx=-lam*inv[:,None]*b
    for it in range(250):
        vx=V(xx); hb=lam*vx[idx,:]
        new=-inv[:,None]*(lam*b+lam*vx-xx@hb)
        err=norm(new-xx);xx=new
        if err<2e-15:break
    else:raise RuntimeError('graph iteration did not converge')
    hb=lam*V(xx)[idx,:]
    metric=np.eye(g.shape[1])+xx.conj().T@xx
    ev,u=la.eigh(metric);ss=(u*np.sqrt(ev))@u.conj().T
    si=(u/np.sqrt(ev))@u.conj().T
    heff=ss@hb@si
    physical_res=e['energy'][:,None]*xx+lam*V(g+xx)-(g+xx)@hb
    out=dict(lambda_=lam,iterations=it+1,graph_residual=norm(physical_res),hermiticity=norm(heff-heff.conj().T))
    for order in [2,3,4]:
        diff=heff-sum(lam**j*e['K'+str(j)] for j in range(2,order+1))
        out['error_after_'+str(order)]=norm(tl(diff))
    return out

def algebra(c):
    N,k=c['N'],c['k'];q=c['parity'];h=c['H'];g=c['G'];src=c['S']
    vals,u=la.eigh(h);spectral_error=0;root_error=0;cubic=0
    for i,m in enumerate(c['Ms']):
        ue=u[:,np.abs(vals-(i+1))<1e-8]
        spectral_error=max(spectral_error,norm(src.conj().T@ue@ue.conj().T@src-m))
        for off in [0,3*N]:
            block=m[off:off+3*N,off:off+3*N]
            dm=(block-np.eye(3*N)/k)/c['rho'];dm[N:2*N,N:2*N]=0
            radius=norm(dm)
            cubic=max(cubic,norm(dm@dm@dm-radius**2*dm)/max(1,radius**3))
    gs,_,p=sc.modes(c['n']);gam=[np.kron(x,np.eye(N)) for x in gs]+[np.kron(p,x) for x in gs]
    free=np.zeros((N*N,N*N),complex)
    for a in range(c['n']):
        ca=(gam[2*a]+1j*gam[2*a+1])/2;cb=(gam[2*c['n']+2*a]+1j*gam[2*c['n']+2*a+1])/2
        free-=np.exp(1j*c['theta'][a])*ca.conj().T@cb+np.exp(-1j*c['theta'][a])*cb.conj().T@ca
    return dict(D=c['D'],code_dimension=N,poles=k,rho2=c['rho']**2,code_error=norm(g.conj().T@g-np.eye(N)),
        kernel_error=norm(h@g),reflection_error=c['reflection_error'],zero_count=int(sum(abs(vals)<1e-9)),gap=float(vals[N]),
        first_order=max(norm(g.conj().T@f@g) for f in c['FL']+c['FR']),
        parity_error=norm(h*q[None,:]-q[:,None]*h),min_weight=float(min(la.eigvalsh(m)[0] for m in c['Ms'])),
        sum_weight_error=norm(sum(c['Ms'])-np.eye(6*N)),source_measure_error=spectral_error,cubic_relative_error=cubic,
        fermion_sign_error=norm(free-c['Kbond']),margins=margins(c['n']))

def main():
    out={'exact_filter_sizes':[],'runs':[]}
    for n in [1,2,4,8,12]:
        assert sc.exact_filter_check(n);out['exact_filter_sizes'].append(n)
    for n in [1,2]:
      for scale in [1.,.5,.25,-1.]:
        start=time.time();sc.RHO_SCALE=scale;c=sc.reflection_cell(n,True);rho2=c['rho']**2
        e=pair_engine(c);rec=dict(n=n,rho_multiplier=scale,dimension=e['dimension'],rho2=rho2,spectral_error=e['spectral_error'])
        for key in ['K2','K3','K4','irreduc','folded']:
            a=e[key]; rec[key]=dict(norm=norm(a),scalar=float(np.trace(a).real/len(a)),traceless_over_rho2=norm(tl(a))/rho2,
                nonGaussian_over_rho2=residual(a,n)/rho2,hermiticity=norm(a-a.conj().T))
        parity=sc.kron_all([c['logical_parity']]*2)
        rec['effective_parity_error']=max(norm(e[key]@parity-parity@e[key]) for key in ['K2','K3','K4'])
        rec['second_order_error']=norm(e['K2']-c['H2'])
        assert rec['second_order_error']<2e-12
        assert rec['effective_parity_error']<2e-12
        if scale==1:
            rec['algebra']=algebra(c)
            aa=rec['algebra']
            assert aa['zero_count']==c['N'] and abs(aa['gap']-1)<1e-11
            assert aa['min_weight']>=3/(4*c['k'])-1e-12
            for test in ['code_error','kernel_error','reflection_error','first_order',
                         'parity_error','sum_weight_error','source_measure_error',
                         'cubic_relative_error','fermion_sign_error']:
                assert aa[test]<2e-12, (test,aa[test])
            assert aa['margins']['root_error']<2e-12
            full=pair_engine(c,False)
            rec['full_vs_reduced']={key:norm(e[key]-full[key]) for key in ['K2','K3','K4']}
            assert max(rec['full_vs_reduced'].values())<3e-12
            rec['finite_coupling']=[exact_graph(e,lam) for lam in [.08,.04,.02]]
            np.savez('coefficients_n'+str(n)+'.npz',K2=e['K2'],K3=e['K3'],K4=e['K4'],rho2=rho2)
            positive={key:e[key].copy() for key in ['K2','K3','K4']}
        if scale==-1:
            rec['rho_sign_errors']={key:norm(e[key]-positive[key]) for key in positive}
        rec['seconds']=time.time()-start;out['runs'].append(rec)
        Path('checks_results.json').write_text(json.dumps(out,indent=2))
        print(json.dumps(rec),flush=True)
    print('PASS',flush=True)

if __name__=='__main__':
    from pathlib import Path
    main()
