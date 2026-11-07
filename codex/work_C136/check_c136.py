"""C136 finite-sector checks. Floating point, not spectral enclosures.

Run small/moderate modes locally. The full grid is a SLURM task.
All writes are relative to this directory; inherited solvers are local copies.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

from bethe_scan import solve, functions
from check_c129 import energy as ed_energy, bond_M
from half_scan import constants as half_constants, root_checks, bulk

DELTAS = [-.9999, -.999, -.99, -.9, -.5, 0.]


def all_constant(gap):
    alpha=np.arccos(1-gap); s=np.sqrt(gap*(2-gap)); B=np.pi/(2*alpha)
    mb=np.sqrt(2)*np.sin(np.pi/4-alpha/2)
    C=4/(np.pi*s*(1-mb))
    H1=2/s**2; H2=8/s**4+4/s**2
    c=2/3*np.log(max(1.,16*C/3))
    t=s/(16*np.pi*np.cosh(c+1)**2)
    kap=min(.5,s/(np.cosh(2*(c+1))+1),t/2*np.exp(-2*c))
    P0=8*C/3; P1=P0*H1; P2=P0*H2
    W=(P2+4*P0+2*P1+H1*(P1+2*P0))/kap
    Vv=np.pi*(4+4*W/3+32*C*H1/(3*kap))
    D=1+3*B/(np.pi*s); D1=1.5*B*H2
    Vclip=2*Vv+8*np.pi*(H1+D1/D)+8*np.pi
    K=max((4*np.pi*D/kap+1)**2,Vclip/8+np.pi*(128*D**2/kap**2+72*D*B/kap+B**2/2))
    return dict(K=K,physical_C=5*K+np.pi**2/4,B=B,Mbeta=mb,C0=C,H1=H1,H2=H2,
                c=c,kappa=kap,P0=P0,P1=P1,P2=P2,W=W,Vv=Vv,D=D,D1=D1,Vclip=Vclip)


def compact_constant(gap, eta):
    s = np.sqrt(gap*(2-gap))
    B = np.pi/(2*np.arccos(1-gap))
    R = np.log(4*B/(np.pi*eta))
    m = s/(np.pi*(np.cosh(2*R)+1))
    S1 = 2/np.pi+1/(np.pi*s**3)
    U1, U2 = 2*B, 8+4*B/(np.pi*s)
    V = 2*R*(U2/m+U1*S1/m**2)
    return dict(K=max((2*B/eta+1)**2, V/8+np.pi*B**2/2),
                B=B, R=R, m=m, S1=S1, U1=U1, U2=U2, V=V)


def fd_checks(L, n, delta, roots, E, order=192):
    """Independent restricted dressing at the finite counting boundary.

    Returns exact-identity residuals and sampled inequality diagnostics.
    Quadrature convergence is separately checked at two orders by the caller.
    """
    assert 0 < n/L < .5
    nu = n/L
    p, th, dp, dt, bare, _ = functions(delta)
    gamma = np.arccos(delta)
    a, b = np.pi*np.sin(gamma), np.pi/gamma
    q, B = 2*gamma/np.pi-1, np.pi/(2*(np.pi-gamma))
    T = lambda x: -dt(x)/(2*np.pi)
    a1 = lambda x: dp(x)/(2*np.pi)
    def z(x):
        x = np.asarray(x)
        return (p(x)-np.sum(th(x[..., None]-roots), axis=-1)/L)/(2*np.pi)
    def sigma(x):
        x = np.asarray(x)
        return (dp(x)-np.sum(dt(x[..., None]-roots), axis=-1)/L)/(2*np.pi)
    upper = max(1., float(roots[-1])+1)
    while z(upper) <= nu/2:
        upper *= 2
        if upper > 30:
            raise RuntimeError('Finite boundary not bracketed')
    Q = brentq(lambda x: z(x)-nu/2, 0., upper, xtol=1e-14)
    nodes, weights = roots_legendre(order)
    y, w = Q*nodes, Q*weights
    A = np.eye(order)-T(y[:, None]-y[None, :])*w
    r, Z, u = np.linalg.solve(A, np.stack([a1(y), np.ones(order),
                                             -a*ap(y, gamma)], axis=1)).T
    def extend(x, vals, source):
        x = np.asarray(x)
        return source(x)+T(x[..., None]-y)@(w*vals)
    rQ = float(extend(Q, r, a1))
    ZQ = float(extend(Q, Z, lambda x: np.ones_like(x)))
    h = -a*rQ/ZQ
    eps_nodes = -a*r-h*Z
    m = float(w@r)
    eq = float(delta/4+w@(bare(y)*r))
    def eps(x):
        return extend(x, eps_nodes, lambda x: bare(x)-h)
    def uprime(x):
        return extend(x, u, lambda x: -a*ap(x, gamma))
    def epp(x):
        x = np.asarray(x)
        return -a*app(x, gamma)+Tp(x[..., None]-y, gamma)@(w*u)
    def sigprime(x):
        x = np.asarray(x)
        return ap(x, gamma)+np.sum(Tp(x[..., None]-roots, gamma), axis=-1)/L
    residual = float(np.sum(eps(roots))-L*(w@(eps_nodes*sigma(y))))
    direct = float(E)-L*(eq+h*(nu-m))
    mass_quad = float(np.sum(extend(roots, Z, lambda x: np.ones_like(x)))/L
                      -w@(Z*sigma(y)))
    grid = np.linspace(-Q,Q,1001)
    sg, up = sigma(grid), uprime(grid)
    variation = float(np.sum(np.abs(np.diff(up/sg))))
    eta = .5-nu
    c = compact_constant(1+delta, eta)
    # Actual fixed-Q constants, as distinct from compact majorants.
    s = np.sin(gamma)
    t0 = abs(delta)/(np.pi*s)
    low = s/(np.pi*(np.cosh(2*Q)+1))
    S1 = 2/np.pi+t0/s**2
    U1 = 2*a*B/np.pi
    U2 = 8*a/np.pi+2*t0*U1
    V = 2*Q*(U2/low+U1*S1/low**2)
    # Independent continuum of the same mass, for the new clipped CDF estimate.
    def cont_at(boundary):
        yy, ww = boundary*nodes, boundary*weights
        rr = np.linalg.solve(np.eye(order)-T(yy[:,None]-yy)*ww, a1(yy))
        return yy, ww, rr
    qtrue = brentq(lambda Q0: np.dot(*cont_at(Q0)[1:])-nu,
                  0., max(Q, 1.)*2, xtol=2e-13)
    yy, ww, rr = cont_at(qtrue)
    def ztrue(x):
        x = np.asarray(x)
        return (p(x)-th(x[...,None]-yy)@(ww*rr))/(2*np.pi)
    cgrid = np.unique(np.r_[np.linspace(-max(Q,qtrue)-1,max(Q,qtrue)+1,1201),roots])
    cc = np.clip(ztrue(cgrid),-nu/2,nu/2)
    cm_before = np.searchsorted(roots,cgrid,side='left')/L-nu/2
    cm_after = np.searchsorted(roots,cgrid,side='right')/L-nu/2
    cdf = max(np.max(abs(cm_before-cc)),np.max(abs(cm_after-cc)))
    ebulk = float(delta/4+ww@(bare(yy)*rr))
    tangent_gap = ebulk-eq-h*(nu-m)
    ca=all_constant(1+delta)
    rg=extend(grid,r,a1)
    r0=1/(2*gamma*np.cosh(b*grid))
    dr=ap(grid,gamma)+Tp(grid[:,None]-y,gamma)@(w*r)
    ddr=app(grid,gamma)+Tpp(grid[:,None]-y,gamma)@(w*r)
    r0p=-b*np.tanh(b*grid)*r0
    r0pp=b*b*(1-2/np.cosh(b*grid)**2)*r0
    decay=np.exp(-1.5*(Q-abs(grid)))
    density_comparison=float(np.min(rg/r0))
    boundary_ratios=[float(np.max(abs(v)/(r0*ca[key]*decay)))
                     for v,key in [(r0-rg,'P0'),(r0p-dr,'P1'),(r0pp-ddr,'P2')]]
    velocity_variation=float(np.sum(abs(np.diff(up/rg))))
    out = dict(Q=Q, Q_true=qtrue, continuum_mass=m, mass_error=nu-m,
               identity_error=abs(direct-residual), mass_identity_error=abs(nu-m-mass_quad),
               direct=direct, midpoint=residual, bulk=ebulk,
               cdf_ratio=float(2*L*cdf/B),
               density_ratio=float(2*L*np.max(abs(sigma(cgrid)-
                   (a1(cgrid)+T(cgrid[:,None]-yy)@(ww*rr))))/(B*2*t0)) if t0 else 0.,
               mass_ratio=float(L*abs(nu-m)/B),
               curvature=float(uprime(Q)/(2*rQ*ZQ**2)), curvature_bound=a,
               tangent_gap=tangent_gap, tangent_gap_bound=a*(nu-m)**2/2,
               sampled_variation=variation, variation_bound=V,
               midpoint_ratio=float(8*L*abs(residual)/V) if V else 0.,
               sampled_epsilon_derivative_ratio=float(np.max(abs(up))/U1),
               sampled_epsilon_second_ratio=float(np.max(abs(epp(grid)))/U2),
               sampled_sigma_derivative_ratio=float(np.max(abs(sigprime(grid)))/S1),
               sampled_min_sigma=float(np.min(sg)), min_sigma_bound=low,
               finite_bulk_scaled_error=float(L*abs(float(E)-L*ebulk)),
               finite_bulk_bound=c['K'], restricted_constant=c,
               all_filling_constant=ca,
               density_comparison=density_comparison,
               boundary_correction_ratios=boundary_ratios,
               velocity_variation=velocity_variation,
               Q_bound_active=bool(L>=2*c['B']/eta), order=order)
    if out['Q_bound_active']:
        assert Q<=c['R']+1e-8
    return out


def ap(x, gamma):
    c = np.cosh(2*x)
    return -2*np.sin(gamma)*np.sinh(2*x)/(np.pi*(c-np.cos(gamma))**2)


def app(x, gamma):
    c, u, d = np.cosh(2*x), np.sinh(2*x), np.cos(gamma)
    return -4*np.sin(gamma)/np.pi*(c/(c-d)**2-2*u*u/(c-d)**3)


def Tp(x, gamma):
    return np.sin(2*gamma)*np.sinh(2*x)/(2*np.pi*(np.sinh(x)**2+np.sin(gamma)**2)**2)


def Tpp(x,gamma):
    den=np.sinh(x)**2+np.sin(gamma)**2
    T=-np.sin(2*gamma)/(2*np.pi*den)
    return T*(2*(np.sinh(2*x)/den)**2-2*np.cosh(2*x)/den)


def scan(mode, index, output):
    if mode=='small':
        lengths, deltas = list(range(2,13)), DELTAS
    elif mode=='moderate':
        lengths, deltas = [16,17,32,33,64,65], DELTAS
    elif mode=='spot':
        lengths, deltas = [128,129,256,257], DELTAS
    else:
        # Full requested all-density grid is cluster-only.
        lengths, deltas = list(range(2,257)), [DELTAS[index]]
    rows=[]
    with Path(output).open('w') as f:
        for delta in deltas:
            for L in lengths:
                ns = (sorted(set([1,2,round(np.sqrt(L)),L//4,L//2-1,L//2]))
                      if mode=='spot' else range(1,L//2+1))
                for n in ns:
                    try:
                        E,diag,x=solve(L,n,delta,return_roots=True)
                        E2,diag2=solve(2*L,2*n,delta)
                        Q=float(2*L*(E2-2*E))
                        k,m=2*n,2*L
                        Ak=np.sin(np.pi*k/m)/np.sin(np.pi/m)
                        trial=2*L*(k-Ak)
                        band=4*(delta+1)*L*n*(2*n-1)/(2*L-1)
                        dilute=2*np.pi**2*n**3/L
                        assert Q<=min(trial,band,dilute)+2e-7
                        row=dict(status='ok',L=L,n=n,delta=delta,E=float(E),E2=float(E2),
                                 Q=Q,free_trial=trial,band=band,dilute=dilute,
                                 min_variational_slack=min(trial,band,dilute)-Q,
                                 diag=diag,diag2=diag2,
                                 correction_indicator=2*L*(diag2['newton_energy_correction']+
                                                          2*diag['newton_energy_correction']))
                        ca=all_constant(1+delta)
                        row['all_filling_constant_ratio']=Q/(5*ca['K'])
                        assert Q<=5*ca['K']+1e-7
                        if L<=12:
                            row['ed_error']=abs(float(E)-ed_energy(L,n,delta))
                            assert row['ed_error']<2e-9
                        if 2*L<=12:
                            row['ed2_error']=abs(float(E2)-ed_energy(2*L,2*n,delta))
                            assert row['ed2_error']<2e-9
                        if delta==0:
                            exact=2*L*np.sin(np.pi*n/L)*np.tan(np.pi/(4*L))/np.cos(np.pi/(2*L))
                            row['free_Q_error']=abs(Q-exact)
                            assert row['free_Q_error']<2e-7
                        if 2*n==L:
                            c=half_constants(delta)
                            row['half_constant_ratio']=Q/c['Q_bound']
                            if mode in ['small','spot']:
                                row['half_root_checks']=root_checks(L,delta,x)
                            assert Q<=c['Q_bound']+1e-7
                        else:
                            c=compact_constant(1+delta,.5-n/L)
                            row['restricted_constant_ratio']=Q/(5*c['K'])
                            assert Q<=5*c['K']+1e-7
                    except Exception as exc:
                        row=dict(status='failed',L=L,n=n,delta=delta,error=str(exc))
                    rows.append(row)
                    f.write(json.dumps(row)+'\n'); f.flush()
            print(f'{mode} delta={delta} finished',flush=True)
    ok=[r for r in rows if r['status']=='ok']
    summary=dict(mode=mode,cases=len(rows),failures=len(rows)-len(ok),
                 largest_Q=max((r['Q'] for r in ok),default=None),
                 largest_at=max(ok,key=lambda r:r['Q']) if ok else None,
                 max_ed_error=max((r.get('ed_error',0) for r in ok),default=0),
                 max_ed2_error=max((r.get('ed2_error',0) for r in ok),default=0),
                 max_free_Q_error=max((r.get('free_Q_error',0) for r in ok),default=0),
                 min_variational_slack=min((r['min_variational_slack'] for r in ok),default=None),
                 max_correction_indicator=max((r['correction_indicator'] for r in ok),default=0),
                 max_residual=max((max(r['diag']['residual'],r['diag2']['residual']) for r in ok),default=0),
                 maxL=max(lengths),max_paired_L=2*max(lengths),
                 per_delta=[dict(delta=d,max_Q=max((r['Q'] for r in ok if r['delta']==d),default=None))
                            for d in deltas])
    Path(output+'.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if summary['failures']:
        raise SystemExit(1)


def intermediate(output, max_order=432, index=None):
    rows=[]
    for d in (DELTAS if index is None else [DELTAS[index]]):
        for L in [8,17,32]:
            for n in sorted(set([1,L//4,(L-1)//2])):
                attempts=[]
                try:
                    E,diag,x=solve(L,n,d,return_roots=True)
                    low=None; convergence=float('inf'); attempts=[]
                    for order in [192,288,432,648,864,1296,1944,2916,4374]:
                        if order>max_order:
                            break
                        try:
                            high=fd_checks(L,n,d,x,E,order)
                            if low is not None:
                                convergence=max(abs(low[k]-high[k]) for k in
                                                ['Q_true','bulk','midpoint','continuum_mass'])
                            attempts.append(dict(order=order,change=convergence if np.isfinite(convergence) else None))
                            high['quadrature_change']=convergence
                            if convergence<2e-9:
                                break
                            low=high
                        except Exception as exc:
                            attempts.append(dict(order=order,error=str(exc)))
                            low=None; convergence=float('inf')
                    assert convergence<2e-8, f'Quadrature unresolved: {convergence}'
                    assert high['identity_error']<2e-8, high['identity_error']
                    assert high['mass_identity_error']<2e-9
                    for key in ['cdf_ratio','density_ratio','mass_ratio','midpoint_ratio',
                                'sampled_epsilon_derivative_ratio','sampled_epsilon_second_ratio',
                                'sampled_sigma_derivative_ratio']:
                        assert high[key]<=1+2e-7,(key,high[key])
                    assert -2e-8<=high['curvature']<=high['curvature_bound']+2e-8
                    assert -2e-8<=high['tangent_gap']<=high['tangent_gap_bound']+2e-8
                    assert high['sampled_variation']<=high['variation_bound']+2e-8
                    assert high['sampled_min_sigma']>=high['min_sigma_bound']-2e-8
                    assert high['finite_bulk_scaled_error']<=high['finite_bulk_bound']+2e-7
                    assert high['density_comparison']>=high['all_filling_constant']['kappa']-2e-8
                    assert max(high['boundary_correction_ratios'])<=1+2e-7
                    assert high['velocity_variation']<=high['all_filling_constant']['Vv']+2e-7
                    row=dict(status='ok',L=L,n=n,delta=d,checks=high,attempts=attempts)
                except Exception as exc:
                    row=dict(status='failed',L=L,n=n,delta=d,error=str(exc),attempts=attempts)
                rows.append(row)
        print(f'intermediate delta={d} finished',flush=True)
    Path(output).write_text(json.dumps(rows,indent=2)+'\n')
    ok=[r['checks'] for r in rows if r['status']=='ok']
    summary=dict(cases=len(rows),failures=len(rows)-len(ok),
                 maxima={k:max((r[k] for r in ok),default=None) for k in
                         ['identity_error','mass_identity_error','quadrature_change','cdf_ratio',
                          'density_ratio','mass_ratio','midpoint_ratio','sampled_epsilon_derivative_ratio',
                          'sampled_epsilon_second_ratio','sampled_sigma_derivative_ratio']})
    Path(output+'.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if summary['failures']:
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['small','moderate','spot','full','intermediate'],required=True)
    parser.add_argument('--index',type=int,choices=range(6))
    parser.add_argument('--max-order',type=int,default=432)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    if args.mode=='intermediate':
        intermediate(args.output,args.max_order,args.index)
    else:
        if args.mode=='full' and args.index is None:
            parser.error('Full grid needs a SLURM array index')
        scan(args.mode,args.index,args.output)
