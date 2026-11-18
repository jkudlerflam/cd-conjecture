"""Small, one-worker exact-control calibrations; no target state is simulated."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys
import json
import time
import itertools
import numpy as np
from scipy import linalg
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits
from records import bare_band,carrier,wrap
from response import (ETA,SIGNS,gaussian,from_response,dimension,dangerous_vertices,
                      response_error_bound,gaussian_error_bounds)

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'codex/work_C191'))
from gauss191 import metric,S4


def chiral_crosscheck(H):
    C=np.array([[-.5,.5,0,0],[0,0,-.5,.5],[.5,.5,0,0],[0,0,.5,.5]])
    W=S4@(2*C.T@H@C)@S4
    old,_,_=metric(W)
    return float(np.max(abs(gaussian(H)['G']-2*C@old@C.T)))


def free_asymmetric():
    """Free spinless-fermion band with p153c's bare hopping and chemical potential.

    This intentionally removes the projector and interactions. Its endpoints,
    occupations and responses are controls, not estimates of the target.
    """
    b=bare_band(); mu=b['mu_eff']; phi=b['phi']
    energy=lambda k:2*np.cos(k-phi)-2*np.cos(2*k)
    velocity=lambda k:-2*np.sin(k-phi)+4*np.sin(2*k)
    ks=[]
    for k in b['k']:
        ks.extend((brentq(lambda p:energy(p)-mu,k-.8,k),
                   brentq(lambda p:energy(p)-mu,k,k+.8)))
    ks=np.array(ks);v=velocity(ks)
    old=np.array([np.ones(4),2*np.cos(ks),2*np.sin(ks),-2*np.cos(2*ks)])
    coefficients=np.array([-mu,np.cos(phi),np.sin(phi),1.])
    null_error=float(np.max(abs(coefficients@old)))
    new=old.copy();new[3]=2*np.sin(2*ks)
    M=np.array([[.5,0,-.5,0],[.5,0,.5,0],[0,.5,0,-.5],[0,.5,0,.5]])
    Z=new@M
    chi=new@np.diag(1/(2*np.pi*abs(v)))@new.T
    fit=from_response(chi,matching=Z)
    exact=np.zeros((4,4))
    for a in range(2):
        left,right=abs(v[2*a]),v[2*a+1]
        exact[a,a]=exact[a+2,a+2]=(right+left)/2
        exact[a,a+2]=exact[a+2,a]=(left-right)/2
    herror=float(np.max(abs(fit['H']-exact)))
    assert null_error<1e-11 and np.linalg.matrix_rank(old,tol=1e-10)==3
    assert np.linalg.matrix_rank(new,tol=1e-10)==4 and herror<1e-11
    assert np.max(abs(fit['G']-np.eye(4)))<1e-11
    return dict(scope='Exact free-band calibration only; not the interacting p153c medium.',
        endpoints=ks.tolist(),endpoint_velocities=v.tolist(),
        free_densities=[float((ks[2*i+1]-ks[2*i])/(2*np.pi)) for i in range(2)],
        old_probe_singular_values=np.linalg.svd(old,compute_uv=False).tolist(),
        null_relation_error=null_error,
        augmented_probe_singular_values=np.linalg.svd(new,compute_uv=False).tolist(),
        matching_Z=Z.tolist(),response=chi.tolist(),H=fit['H'].tolist(),
        H_error=herror,G_error=float(np.max(abs(fit['G']-np.eye(4)))),
        signed_velocities=fit['signed_velocities'].tolist(),
        C191_covariance_difference=chiral_crosscheck(exact))


def hessian(fun,h):
    E=np.eye(4)*h;out=np.empty((4,4));zero=fun(np.zeros(4))
    for i in range(4):
        out[i,i]=(fun(E[i])+fun(-E[i])-2*zero)/h**2
        for j in range(i):
            out[i,j]=out[j,i]=(fun(E[i]+E[j])-fun(E[i]-E[j])-
                               fun(-E[i]+E[j])+fun(-E[i]-E[j]))/(4*h*h)
    return out


def hard_rods():
    """Exact continuum hard-rod equation of state, then a compact canonical shear.

    The shear is a coordinate change, not a candidate Hamiltonian. It exercises
    mixed covariance entries with a known answer and unequal mode velocities.
    """
    rho=np.array([.16,.09]);a=np.array([.6,1.3]);mass=np.array([1.,1.4])
    K=(1-a*rho)**2;u=np.pi*rho/(mass*K)
    H0=np.diag(np.r_[u/K,u*K]);G0=np.diag(np.r_[K,1/K])
    F=np.array([[0.,2.],[-2.,0.]])
    T=np.block([[np.eye(2),np.zeros((2,2))],[F,np.eye(2)]])
    Ti=np.linalg.inv(T);H=Ti.T@H0@Ti;G=T@G0@T.T
    assert np.max(abs(T@ETA@T.T-ETA))<1e-14
    source_inverse=SIGNS@Ti@SIGNS
    e0=lambda r:np.pi**2*r**3/(6*mass*(1-a*r)**2)
    mu=np.pi**2*rho**2*(3-a*rho)/(6*mass*(1-a*rho)**3)
    def energy(X):
        old=source_inverse@X;r=rho+old[:2];w=old[2:]
        return float(np.sum(e0(r)-e0(rho)-mu*old[:2]+np.pi**2*r*w*w/(2*mass)))
    fits=[]
    for h in (1e-3,5e-4,2.5e-4):
        curvature=hessian(energy,h);chi=np.linalg.inv(curvature)
        fit=from_response(chi)
        fits.append(dict(step=h,H_max_error=float(np.max(abs(fit['H']-H))),
                         G_max_error=float(np.max(abs(fit['G']-G))),
                         signed_velocity_max_error=float(np.max(abs(
                             fit['signed_velocities']-gaussian(H)['signed_velocities'])))))
    assert fits[-1]['H_max_error']<2e-6 and fits[-1]['G_max_error']<2e-5
    exactchi=SIGNS@np.linalg.inv(H)@SIGNS/np.pi
    Z=np.array([[1.,1.,0,0],[1.,-1.,.2,0],[0,.3,1,0],[.2,0,0,1.]])
    contact=np.diag([.03,.02,.01,.04])
    obs=Z@exactchi@Z.T+contact
    matched=from_response(obs,matching=Z,contact=contact)
    matcherror=float(np.max(abs(matched['H']-H)))
    assert matcherror<1e-12
    # Perturb the *matched* response by a stated norm; check the analytic bounds.
    E=np.array([[1.,.2,0,0],[.2,-.3,.1,0],[0,.1,.4,.2],[0,0,.2,-.7]])
    E=E/np.linalg.norm(E,2)*1e-6
    moved=from_response(exactchi+E);dh=response_error_bound(exactchi,1e-6)
    bounds=gaussian_error_bounds(H,dh)
    actualG=float(np.linalg.norm(moved['G']-G,2))
    actualv=float(np.max(abs(moved['signed_velocities']-gaussian(H)['signed_velocities'])))
    assert np.linalg.norm(moved['H']-H,2)<=dh
    assert actualG<=bounds['covariance_error_norm'] and actualv<=bounds['signed_velocity_error']
    return dict(scope='Exact interacting control and integer field-coordinate change only.',
                rho=rho.tolist(),rod_length=a.tolist(),mass=mass.tolist(),K=K.tolist(),u=u.tolist(),
                H=H.tolist(),G=G.tolist(),finite_difference_checks=fits,
                matching_contact_recovery_error=matcherror,
                C191_covariance_difference=chiral_crosscheck(H),
                response_error_norm=1e-6,H_error_bound=dh,propagated_bounds=bounds,
                actual_covariance_error=actualG,actual_signed_velocity_error=actualv)


def lattice_checks():
    G=np.eye(4);ks=(-2.9,-.23);rho=(.11,.09)
    found=dangerous_vertices(G,ks,rho,dk=(.001,.002),drho=(.0001,.0002))
    keys={(tuple(r['n']),tuple(r['m'])) for r in found}
    # Independent enumeration on a larger box checks the proven radius code.
    brute=set()
    for v in itertools.product(range(-4,5),repeat=4):
        n=v[:2];m=v[2:];first=next((x for x in v if x),0)
        if first>0 and sum(n)%2==0 and dimension(G,n,m)<=2+1e-12:brute.add((n,m))
    assert keys==brute
    for n in ((2,0),(0,2),(1,1),(2,2)):
        assert (n,(0,0)) in keys
    # A mixed self-non-null vertex must be retained and flagged, not discarded.
    nonnull=[r for r in found if not r['self_null']]
    assert nonnull
    return dict(parity_even_vertices_with_x_at_most_two=len(found),
                matches_larger_box=True,non_null_count=len(nonnull),
                example_non_null=nonnull[0],scope='Enumeration calibration; no target covariance supplied.')


def xx_g(m):
    if m==0:return .5
    idx=np.arange(m);n=idx[:,None]-idx[None,:]+1
    G=np.zeros((m,m));odd=n%2!=0
    G[odd]=2*(-1.)**((n[odd]-1)//2)/(np.pi*n[odd])
    return float(linalg.det(G)/2)


def phase_control():
    # Exact XX x XX pair envelope, gauged to the declared pair momentum.
    r=np.arange(1,201);g=np.array([xx_g(m) for m in range(102)])
    assert abs(g[1]-1/np.pi)<1e-14 and abs(g[2]-2/np.pi**2)<1e-14
    Q=-np.pi+.0085;exact=np.exp(1j*Q*r)*g[r//2]*g[(r+1)//2]
    mean=.17+.09j;contaminant=.05*np.exp(.31j*r)/r**.3
    raw=exact+abs(mean)**2+contaminant
    connected=raw-abs(mean)**2
    fits=[]
    for window in ((12,32),(16,40),(32,96),(48,160)):
        f=carrier(r,connected,-np.pi,*window)
        f['momentum_error']=float(wrap(f['k']-Q));fits.append(f)
    assert max(abs(f['momentum_error']) for f in fits)<1e-4
    # A nearby harmonic is intentionally unresolved: a tiny regression residual
    # does not measure its bias. This must not be reported as an uncertainty bar.
    near=connected+.45*exact*np.exp(.012j*r)
    bad=carrier(r,near,-np.pi,16,40);bad['momentum_bias']=float(wrap(bad['k']-Q))
    assert abs(bad['momentum_bias'])>100*bad['formal_slope_error']
    return dict(exact_Q=Q,XX_two_leg_fits=fits,nearby_harmonic_failure=bad,
                scope='Known physical XX pair control plus explicit synthetic contamination; not target exponents.')


def main():
    start=time.monotonic()
    with threadpool_limits(1):
        out=dict(free_band=free_asymmetric(),hard_rods=hard_rods(),
                 lattice=lattice_checks(),phase=phase_control())
    out['runtime_seconds']=time.monotonic()-start
    (HERE/'calibration_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
