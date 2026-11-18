"""Bounded, one-worker half-density contact response and clock retardation checks."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from itertools import combinations
import importlib.util
import json
import time
import numpy as np
import scipy.linalg as la
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import mpmath as mp

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('mobile_control',HERE.parent/'markerless_mobile_control_20261118'/'checks.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def uniform_fermions(N):
    L=2*N
    k=2*np.pi*(np.arange(L)-L/2+(.5 if N%2==0 else 0))/L
    e=-2*np.cos(k)
    occ=np.flatnonzero(e<0); emp=np.flatnonzero(e>0)
    assert len(occ)==N
    empty=set(map(int,emp))
    result=0.
    for h1,h2 in combinations(occ,2):
        for p1 in emp:
            p2=int((h1+h2-p1)%L)
            if p2 not in empty or p1>=p2:continue
            amp=2/L*(np.cos(k[p1]-k[h1])-np.cos(k[p1]-k[h2]))
            gap=e[p1]+e[p2]-e[h1]-e[h2]
            result+=amp**2/gap/N
    boost=8*np.sin(np.pi/L)**4/(L**3*(np.sin(np.pi/L)+np.sin(3*np.pi/L)))
    return dict(N=N,gamma_zero=result,one_current_branch=boost,
                scaled_current_branch=N**6*boost)


def response(N):
    H,C=old.bath(N,N); H=H.astype(float)
    es,us=sla.eigsh(H,k=1,which='SA',tol=2e-13,v0=np.ones(H.shape[0]))
    E=es[0];o=us[:,0]
    alpha=old.alpha(N,N)
    X=(C.astype(float)-alpha)*o[:,None]
    op=sla.LinearOperator(H.shape,matvec=lambda v:H@v-E*v+o*(o@v))
    y,info=sla.cg(op,X[:,0],rtol=2e-12,atol=1e-14,maxiter=3000)
    assert info==0
    col=X.T@y
    gam=np.fft.fft(col).real
    assert gam.min()>-2e-11
    assert abs(gam[0]-uniform_fermions(N)['gamma_zero'])<2e-10
    residual=la.norm(op@y-X[:,0])
    return dict(N=N,V=N,dimension=H.shape[0],gamma=list(gam),
                gamma_max=float(gam.max()),residual=float(residual),alpha=alpha)


def bethe_contacts(N):
    """Exact finite Bethe spectral representation, bounded checks N=4,6 only."""
    V=N; L=N+V
    integers=np.arange(L)-L//2+(.5 if V%2==0 else 0)
    ground=set(np.argsort(abs(integers))[:V])
    pg=2*np.pi*integers[sorted(ground)]/L
    E0=-2*sum(np.cos(pg))
    gam=np.zeros(N); beyond=np.zeros(N); counts=np.zeros(N,dtype=int)
    norm=N*L**(V-1)
    y=np.arange(1,L-1)
    for inds in combinations(range(L),V):
        momentum=int(round(sum(integers[list(inds)])))%L
        if momentum>=N:continue
        q=2*np.pi*momentum/N
        p=(2*np.pi*integers[list(inds)]+q)/L
        energy=-2*sum(np.cos(p))
        counts[momentum]+=1
        overlap=np.exp(-1j*np.outer(pg,y))@np.exp(1j*np.outer(y,p))
        form=np.linalg.det(overlap)/norm
        if abs(energy-E0)<1e-12:
            assert abs(form-old.alpha(N,N))<2e-13
            continue
        value=N*abs(form)**2/(energy-E0)
        gam[momentum]+=value
        if len(set(inds)-ground)>2:beyond[momentum]+=value
    from math import comb
    assert sum(counts)==comb(L-1,V)
    direct=response(N)
    assert np.max(abs(gam-direct['gamma']))<3e-11
    assert max(beyond)>1e-8
    return dict(N=N,states=int(sum(counts)),gamma=list(gam),
                gamma_from_more_than_two_changes_of_Bethe_integers=list(beyond),
                max_difference_from_gap_resolvent=float(np.max(abs(gam-direct['gamma']))))


def clocks():
    nu=old.selector(10)
    alpha=.5-2/np.pi**2
    forward=2/np.pi*(.25+1/np.pi**2)**2
    data={}
    for b in (9,10,11):
        c=old.cell(b,nu); h=c['sources'].sum(axis=0)
        es,U=la.eigh(h); v=U[:,0]
        D=np.einsum('s,sij->ij',np.arange(b+1)-1,c['sources'])
        X=np.diag(np.arange(b)-.5)
        assert la.norm(D-(X@h+h@X)/2)<1e-13
        dip=U[:,1:].T@D@v; gaps=es[1:]-es[0]
        spectral_moment=float(np.sum(gaps**3*abs(dip)**2))
        theta=np.pi/(b+1);epsilon=es[0]
        moment_formula=8*np.sin(theta)**2/(b+1)*(np.cos(theta)*(1+2*epsilon**2)+2*epsilon)
        assert abs(spectral_moment-moment_formula)<2e-13
        log_coefficient=forward*alpha**3/np.pi*spectral_moment
        data[b]=dict(b=b,gaps=list(gaps),dipole_sq=list(abs(dip)**2),
                     dipole_norm_sq=float(np.sum(abs(dip)**2)),
                     cubic_gap_moment=spectral_moment,
                     closed_cubic_gap_moment=moment_formula,
                     acoustic_J5_log_coefficient=log_coefficient)
    patches=[]
    for name,A,B in [('balanced',[10,10],[9,11]),('more_heads',[10]*9,[9]*10),('fewer_heads',[10]*11,[11]*10)]:
        delta=sum(data[b]['acoustic_J5_log_coefficient'] for b in B)-sum(data[b]['acoustic_J5_log_coefficient'] for b in A)
        patches.append(dict(name=name,J5_log_difference=delta))
    return dict(cells=list(data.values()),patches=patches)


def constants():
    mp.mp.dps=50
    i2=mp.quad(lambda y:y*y/mp.cosh(y)**2,[0,1,mp.inf])
    i4=mp.quad(lambda y:y*y/mp.cosh(y)**4,[0,1,mp.inf])
    assert abs(i2-mp.pi**2/12)<mp.mpf('1e-45')
    assert abs(i4-(mp.pi**2-6)/18)<mp.mpf('1e-45')
    gamma0=(24-mp.pi**2)/(18*mp.pi**3)
    C=2/mp.pi*(mp.mpf(1)/4+1/mp.pi**2)**2
    I=3*mp.log(2)/mp.pi
    I1=(9*mp.log(2)+3)/mp.pi**2
    I2=mp.mpf(9)/(2*mp.pi)+18*mp.log(2)/mp.pi**3
    e2bound=mp.mpf(4)/3*(mp.mpf(1)/4+I+2*I1+I2)+4/(3*mp.sqrt(3))*(mp.mpf(1)/4+I+I1)
    assert e2bound<8
    return dict(gamma_zero_infinite=float(gamma0),forward_jump=float(C),
                hydro_small_nonzero_q=float(gamma0+C),current_scaled_limit=float(mp.pi**3/32),
                homogeneous_XXZ_second_derivative_bound=float(e2bound))


def projected_zero_mode():
    nu=old.selector(10);eb=old.cell(10,nu)['eps']/10
    gamma=constants()['gamma_zero_infinite']
    out=[]
    for name,A,B in [('balanced',[10,10],[9,11]),('more_heads',[10]*9,[9]*10),('fewer_heads',[10]*11,[11]*10)]:
        delta=sum(old.cell(b,nu)['eps'] for b in B)-sum(old.cell(b,nu)['eps'] for b in A)
        out.append(dict(name=name,source_sum_difference=delta,
                        thermodynamic_uniform_cross_term=-2*eb*gamma*delta))
    return out


def projected_constants():
    """Local source norms; keep the two exterior contacts distinct."""
    nu=old.selector(10)
    cells={b:old.cell(b,nu) for b in (9,10,11)}
    means={}
    partial_norms={}
    for b,c in cells.items():
        es,U=la.eigh(c['sources'].sum(axis=0));v=U[:,0]
        a=np.einsum('i,sij,j->s',v,c['sources'],v)
        xi=np.einsum('sij,j->si',c['sources'],v)-a[:,None]*v
        assert la.norm(xi.sum(axis=0))<1e-13
        means[b]=a
        partial_norms[b]=float(la.norm(np.cumsum(xi,axis=0)[:-1],axis=1).sum())
    def unwrapped(word):
        out=np.zeros(sum(word)+1);start=0
        for b in word:
            out[start:start+b+1]+=means[b]
            start+=b
        return out
    periodic=means[10][:-1].copy();periodic[0]+=means[10][-1]
    harmonics=abs(np.fft.fft(periodic)/10)
    patches=[]
    for name,A,B in [('balanced',[10,10],[9,11]),('more_heads',[10]*9,[9]*10),('fewer_heads',[10]*11,[11]*10)]:
        d=unwrapped(B)-unwrapped(A)
        norm1=float(abs(d).sum())
        clock_norm=float(sum(partial_norms[b]**2 for b in A+B))
        patches.append(dict(name=name,source_l1=norm1,
                            source_sum=float(d.sum()),
                            changed_clock_B_squared_sum=clock_norm,
                            coefficient_of_P=2*norm1,
                            coefficient_of_T_after_M_le_4T=norm1**2+4*clock_norm))
    return dict(periodic_source_fourier_magnitudes=list(harmonics),
                clock_partial_sum_norms={str(b):value for b,value in partial_norms.items()},
                patches=patches)


def retarded_kernel_audit():
    N=6;word=[3,3];nu=old.selector(10);J=.2
    H,C=old.bath(N,N);H=H.astype(float)
    es,Ub=la.eigh(H.toarray());o=Ub[:,0]
    hs=old.sources(word,nu);Hm=sum(hs)
    em,Uc=la.eigh(Hm.toarray());v=Uc[:,0]
    alpha=old.alpha(N,N)
    residual=sum((sp.kron(sp.diags(C[:,j].astype(float)-alpha),hs[j],format='csr')
                  for j in range(N)),sp.csr_matrix((len(C)*len(v),)*2))
    psi=np.kron(o,v)
    amplitudes=Ub.T@(residual@psi).reshape(len(o),len(v))@Uc
    assert np.max(abs(amplitudes[0,:]))<1e-12
    weights=abs(amplitudes[1:,:])**2
    omega=(es[1:]-es[0])[:,None];d=alpha*(em-em[0])[None,:]
    d[:,0]=0
    second=-float(np.sum(J*J*weights/(omega+J*d)))
    curvature0=2*float(np.sum(weights/omega))
    curvatureJ=2*float(np.sum(weights*omega**2/(omega+J*d)**3))
    assert 0<=curvatureJ<=curvature0
    reference=sp.kron(H,sp.eye(len(v)),format='csr')+J*alpha*sp.kron(sp.eye(len(o)),Hm,format='csr')
    E0=es[0]+J*alpha*em[0]
    estimates=[]
    for eta in (.002,.004):
        e=[]
        for sign in (-1,1):
            e.append(float(sla.eigsh(reference+sign*eta*J*residual,k=1,which='SA',tol=1e-13,v0=psi)[0][0]))
        estimates.append((sum(e)/2-E0)/eta**2)
    extrapolated=(4*estimates[0]-estimates[1])/3
    assert abs(extrapolated-second)<2e-8
    return dict(N=N,V=N,word=word,J_over_t=J,dimension=len(psi),
                exact_second_eta_coefficient=second,finite_eta_estimates=estimates,
                extrapolated_eta_coefficient=extrapolated,
                absolute_error=abs(extrapolated-second),
                magnitude_J_curvature_at_zero=curvature0,
                magnitude_J_curvature_at_J=curvatureJ)


def main():
    started=time.perf_counter()
    out=dict(constants=constants(),fermion_zero_mode=[uniform_fermions(N) for N in (4,6,8,10,20,40,80)],
             gap_responses=[response(N) for N in (4,6,8,10)],
             bethe_contact_audit=[bethe_contacts(N) for N in (4,6)],
             projected_zero_mode=projected_zero_mode(),projected_constants=projected_constants(),
             retarded_kernel_audit=retarded_kernel_audit(),
             clock_retardation=clocks())
    out['versions']=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__)
    out['seconds']=time.perf_counter()-started
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
