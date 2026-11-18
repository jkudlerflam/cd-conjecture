#!/usr/bin/env python3
"""Bounded one-worker operator/continuum checks; never loads MPS or runs DMRG."""
import itertools,json,time,hashlib,os,importlib.util
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def covariance(h,metric):
    w,u=np.linalg.eigh(h);assert min(w)>0
    root=(u*np.sqrt(w))@u.T;inv=(u/np.sqrt(w))@u.T
    w,u=np.linalg.eigh(root@metric@root)
    return np.real_if_close(inv@(u*np.abs(w))@u.conj().T@inv)


def oscillator_check(h):
    # A single real spatial Fourier mode: theta_cos=-p_sin, theta_sin=p_cos.
    a,b,c=h[:2,:2],h[2:,2:],h[:2,2:]
    z=np.zeros((2,2));hq=np.block([[a,z],[z,a]]);hp=np.block([[b,z],[z,b]])
    cross=np.block([[z,-c],[c,z]])
    osc=np.block([[hq,cross],[cross.T,hp]])
    j=np.block([[np.zeros((4,4)),np.eye(4)],[-np.eye(4),np.zeros((4,4))]])
    g=covariance(osc,1j*j)
    # twice the oscillator covariance; same convention as logarithmic G.
    t=np.zeros((4,8));t[:2,:2]=np.eye(2);t[2:,6:]=-np.eye(2)
    return t@g@t.T


def gaussian_checks():
    z=np.zeros((2,2));eta=np.block([[z,np.eye(2)],[np.eye(2),z]])
    j=np.block([[z,np.eye(2)],[-np.eye(2),z]])
    h=np.eye(4);h[0,3]=h[3,0]=h[1,2]=h[2,1]=.6
    correct=covariance(h,eta);old=covariance(h,1j*j)
    assert np.max(abs(correct-np.eye(4)))<1e-12
    spec=importlib.util.spec_from_file_location('gauss191',ROOT/'codex/work_C191/gauss191.py')
    c191=importlib.util.module_from_spec(spec);spec.loader.exec_module(c191)
    transform=np.array([[-1,-1,1,1],[1,1,1,1],[-1,1,1,-1],[1,-1,1,-1]])/np.sqrt(2)
    ti=np.linalg.inv(transform)
    rng=np.random.default_rng(180);errors=[];chiral_errors=[]
    for _ in range(20):
        x=rng.normal(size=(4,4));h=x.T@x+np.eye(4)
        g=covariance(h,eta)
        errors.append(float(np.max(abs(g-oscillator_check(h)))))
        wf=2*ti.T@h@ti
        gc,_,_=c191.metric(c191.S4@wf@c191.S4)
        chiral_errors.append(float(np.max(abs(g-2*ti@gc@ti.T))))
        assert np.max(abs(g@eta@g-eta))<1e-11
    assert max(errors)<1e-11 and max(chiral_errors)<1e-11
    # A symmetry-allowed counterexample to checking only pure locking vertices.
    q=np.array([[1.1,.9],[.9,1.1]]) # coordinates (theta_c,phi_s)
    g=np.zeros((4,4));g[np.ix_([2,1],[2,1])]=q
    g[np.ix_([0,3],[0,3])]=np.linalg.inv(q)
    h=np.linalg.inv(g);assert np.max(abs(covariance(h,eta)-g))<1e-11
    dims={f'{r},{s}':float(2*np.array([r,-s])@q@np.array([r,-s]))
          for r,s in [(1,0),(0,1),(1,1),(1,-1)]}
    return dict(free_opposite_drift_correct_G=correct.tolist(),
        C180_oscillator_formula_G=old.tolist(),
        free_quartet_dimension_correct=float(2*correct[2,2]),
        free_quartet_dimension_C180=float(2*old[2,2]),
        independent_real_Fourier_oscillator_max_error=max(errors),
        agreement_with_existing_C191_chiral_code_max_error=max(chiral_errors),
        mixed_vertex_example=dict(Q=q.tolist(),H=h.tolist(),dimensions=dims,
                                  scope='Allowed Gaussian example, not fitted V4 parameters.'))


def vertex_checks():
    out=[];expected=[];points=[(np.sqrt(2)/5,np.sqrt(3)/9),(.071234,.293817),(-.19371,.17239)]
    for a,b,m,n in itertools.product(range(-4,5),repeat=4):
        if not (a or b or m or n) or (a+b)%2:continue
        phases=[np.exp(1j*np.pi*((a-b)*k+b+rho*(m+n))) for k,rho in points]
        if all(abs(z-1)<1e-10 for z in phases):out.append([a,b,m,n])
        if a==b and a%2==0 and m==-n:expected.append([a,b,m,n])
    assert out==expected
    # On the generic lattice A and Sigma phases are i**(a+b)=1.
    assert all((x[0]+x[1])%4==0 for x in out)
    return dict(range_each_integer=[-4,4],generic_vertices=out,phase_test_points_k_over_pi_rho=points,
                exact_rule='nA=nB=2r, mA=-mB=s; all integer r,s')


def microscopic_checks():
    L=8;D=1<<L;st=np.arange(D);bits=(st[:,None]>>np.arange(L))&1
    z=1-2*bits;N=bits.sum(1);H=np.zeros((D,D),complex)
    for j in range(L):
        k=(j+1)%L;to=st^(1<<j)^(1<<k)
        H[to,st]+=.5*(1+z[:,j]*z[:,k])+1j*(z[:,k]-z[:,j])/16
    H[st,st]+=(z*np.roll(z,-1,axis=1)).sum(1)-1.75*z.sum(1)
    R=np.exp(-1j*(bits@(np.pi/2+np.pi*np.arange(L))))
    inv=np.array([sum(int(bits[s,j])<<((-j)%L) for j in range(L)) for s in st])
    I=np.eye(D)[inv];Sigma=R[:,None]*I
    err=dict(theta=float(np.max(abs(I@H.conj()@I-H))),
             sigma=float(np.max(abs(Sigma@H@Sigma.conj().T-H))),
             sigma_squared_parity=float(np.max(abs(Sigma@Sigma-np.diag((-1.)**N)))))
    for r in (2,3):
        pair=np.zeros((D,D),complex)
        for j in range(L):
            k=(j+r)%L;mask=(bits[:,j]+bits[:,k])==0
            pair[(st[mask]^(1<<j)^(1<<k)),st[mask]]+=1
        err[f'uniform_pair_r{r}']=float(np.max(abs(Sigma@pair@Sigma.conj().T-(-1)**(r+1)*pair)))
    assert max(err.values())<1e-12
    return dict(L=L,errors=err)


def virtual_paths():
    L=10;D=1<<L
    adj=lambda s:sum(((s>>j)&1)*((s>>((j+1)%L))&1) for j in range(L))
    counts=[adj(s) for s in range(D)];P=[s for s in range(D) if counts[s]==0];ps=set(P)
    pairs=[];hops=[]
    for s in range(D):
        bp=[];wp=[]
        for j in range(L):
            k=(j+1)%L;x=(s>>j)&1;y=(s>>k)&1;t=s^(1<<j)^(1<<k)
            if x==y:bp.append(t)
            else:wp.append((t,1j if x else -1j))
        pairs.append(bp);hops.append(wp)
    g3={};mixed={}
    for s in P:
        for t in pairs[s]:
            if not counts[t]:continue
            for u in pairs[t]:
                if not counts[u]:continue
                for v in pairs[u]:
                    if v in ps and v.bit_count()==s.bit_count()+2:
                        key=(v,s);g3[key]=g3.get(key,0)+1/(counts[t]*counts[u])
            for v,w in hops[t]:
                if v in ps and v.bit_count()==s.bit_count()+2:
                    key=(v,s);mixed[key]=mixed.get(key,0)-w/counts[t]
        for t,w in hops[s]:
            if not counts[t]:continue
            for v in pairs[t]:
                if v in ps and v.bit_count()==s.bit_count()+2:
                    key=(v,s);mixed[key]=mixed.get(key,0)-w/counts[t]
    expected={}
    for s in P:
        for j in range(L):
            k=(j+3)%L
            if not (s&(1<<j) or s&(1<<k)):
                v=s|(1<<j)|(1<<k)
                if v in ps:expected[v,s]=expected.get((v,s),0)+2/3
    err=max(abs(g3.get(k,0)-expected.get(k,0)) for k in g3.keys()|expected.keys())
    # Retain any additional correlated paths rather than assuming their absence.
    extra=[dict(initial=s,final=v,coefficient=g3.get((v,s),0)-expected.get((v,s),0))
           for v,s in g3.keys()|expected.keys() if abs(g3.get((v,s),0)-expected.get((v,s),0))>1e-12]
    return dict(L=L,blockade_dimension=len(P),
        imaginary_hop_pair_second_order_max=max([abs(x) for x in mixed.values()]+[0]),
        third_order_vs_projected_separation3_max_error=err,extra_correlated_paths=extra,
        g3_units='Delta^3/V^2',scope='Exact finite-ring path sum about V sum nn, not a finite-V error bound.')


def main():
    start=time.monotonic()
    with threadpool_limits(limits=1):
        out=dict(gaussian=gaussian_checks(),vertices=vertex_checks(),
                 microscopic=microscopic_checks(),virtual_paths=virtual_paths())
    paths=['codex/work_C180/PROOF.md','codex/work_C180/optable.py','codex/work_C179/PROOF.md',
           'codex/work_C189/PROOF.md','codex/work_C191/PROOF.md','codex/work_C191/gauss191.py',
           'contrib/v4_interpretation_audit_20261117/audit.json']
    out['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    out['runtime_seconds']=time.monotonic()-start;out['workers']=out['BLAS_threads']=1
    (HERE/'checks.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(seconds=out['runtime_seconds'],gaussian=out['gaussian'],
                         virtual_paths=out['virtual_paths']),indent=2))


if __name__=='__main__':main()
