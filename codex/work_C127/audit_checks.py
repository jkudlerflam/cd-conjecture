"""Assertions and independent analytic checks for C127. Run after checks.py."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from checks import boson_log, boson_sectors

ROOT = Path(__file__).resolve().parent
x = json.loads((ROOT/'checks.json').read_text())
tol = 2e-8
for r in x['XX_exact']:
    assert r['sharp_margin'] >= -tol
    assert r['boson_margin'] >= -tol
    assert -tol <= r['product_error'] <= r['product_tail']+tol
for r in x['XX_matsubara']:
    assert r['LR_margin'] >= -1e-12
assert any(r['global_linear_margin'] < -1 for r in x['XX_matsubara'])
for r in x['QTM']:
    assert r['checks']['commutator'] < 1e-12
    assert abs(r['checks']['contraction_error']) < tol
    assert r['checks']['min_Ke_relative'] > -1e-11
    assert r['checks']['min_U_relative'] > -1e-11
    for t in r['traces']:
        assert t['candidate_margin'] > t['discarded_trace_bound']
for key in ['boson_inequalities', 'theta_inequality', 'Gibbs_comparison', 'nonlinear_product']:
    for r in x[key]:
        assert all(v >= -1e-10 for k,v in r.items() if k.endswith('margin'))

resolvent = []
for nu in [-3., -.2, .2, 3.]:
    for spatial in [1, 2, 3, 5]:
        f = lambda k: np.cos(spatial*k)/(1j*nu-np.cos(k))
        val = (quad(lambda k:f(k).real, 0, 2*np.pi, epsabs=1e-12)[0]
               +1j*quad(lambda k:f(k).imag, 0, 2*np.pi, epsabs=1e-12)[0])/(2*np.pi)
        exact = -1j*np.sign(nu)/np.sqrt(1+nu*nu)*(-1j*np.sign(nu)*np.exp(-np.arcsinh(abs(nu))))**spatial
        residual = float(abs(val-exact))
        assert residual < 1e-11
        resolvent.append(dict(nu=nu, x=spatial, residual=residual))

bubble = []
V = 2*np.e
for beta in [2., 10., 100.]:
    for k in [1, 2, 5]:
        omega = 2*np.pi*k/beta
        if omega>=V: continue
        nu = (2*np.arange(-10000,10001)+1)*np.pi/beta
        other = nu-omega
        energy = np.arcsinh(abs(nu))+np.arcsinh(abs(other))
        den = beta*np.sqrt((1+nu*nu)*(1+other*other))
        for spatial in [1, 3, 5, 11]:
            positive = float(np.sum(np.exp(-spatial*energy)/den))
            G = -1j*np.sign(nu)/np.sqrt(1+nu*nu)*(-1j*np.sign(nu)*np.exp(-np.arcsinh(abs(nu))))**spatial
            Gp = -1j*np.sign(other)/np.sqrt(1+other*other)*(-1j*np.sign(other)*np.exp(-np.arcsinh(abs(other))))**spatial
            reflected = np.sum(G*Gp).real/beta
            bound = .5*V/(omega*(V-omega))*np.exp(-omega*spatial/V)
            tail = 2*beta**(2*spatial+1)/(2.**(2*spatial)*np.pi**(2*spatial+2)
                    *(2*spatial+1)*10000.**(2*spatial+1))
            assert abs(reflected-positive) < 1e-12
            assert positive+tail<=bound
            bubble.append(dict(beta=beta,k=k,x=spatial,reflected_value=positive,
                               tail_bound=float(tail),LR_bound=float(bound),margin=float(bound-positive-tail)))

# Convex midpoint sum and the integral used in the exact XX proof.
integrals = []
for L in [2., 3., 4., 8., 16.]:
    val = quad(lambda t: np.cosh(t)*np.log1p(np.exp(-L*t)),0,50,epsabs=1e-12)[0]
    exact = (np.pi/np.sin(np.pi/L)-L)/2
    assert abs(val-exact)<1e-11
    for beta in [.1,1.,10.,100.]:
        n=np.arange(20000)
        midpoint=np.log1p(np.exp(-L*np.arcsinh((2*n+1)*np.pi/beta))).sum()
        right=np.log1p(np.exp(-L*np.arcsinh((2*n+2)*np.pi/beta))).sum()
        cap=beta/(2*np.pi)*val
        assert max(midpoint,right)<=cap+1e-10
    integrals.append(dict(L=L,value=float(val),exact=float(exact)))

# Independent partition coefficients check the square-root Fock count.
N=100; S=6
coef=np.zeros(N+1,dtype=object);coef[0]=1
for k in range(1,N+1):
    for species in range(S):
        for n in range(k,N+1):coef[n]+=coef[n-k]
cum=0;count=[]
for n,c in enumerate(coef):
    cum+=int(c)
    upper=2*np.pi*np.sqrt(S*n/6)
    margin=float(upper-np.log(float(cum)))
    assert margin>=-1e-12
    if n in [0,1,2,10,50,100]:count.append(dict(level=n,log_count=float(np.log(float(cum))),margin=margin))

sector_consistency=[]
for a in [.1,1.,10.,30.]:
    for m in [2,3,4]:
        sectors=boson_sectors(a,m)
        top=max(s for s in sectors if s is not None)
        total=top+np.log(sum(np.exp(s-top) for s in sectors if s is not None))
        expected,_=boson_log(a,6)
        err=abs(total-expected)
        assert err<1e-10
        sector_consistency.append(dict(a=a,m=m,residual=float(err)))

failures=dict(
    worst_sector=min((v,r['delta'],r['beta'],r['m'],t['L'],q)
                    for r in x['QTM'] for t in r['traces']
                    for q,v in enumerate(t['sector_candidate_margins']) if v is not None),
    worst_finite_step_low_cone=min((s['finite_m_cone_margin'],r['delta'],r['beta'],r['m'],s['q'])
                    for r in x['QTM'] for s in r['spectrum']
                    if s['q'] and s['omega']<4*np.e*r['M']),
    K_zero_classical=x['classical_K_necessity'][-1])
out=dict(resolvent=resolvent,bubble=bubble,integrals=integrals,Fock_count=count,
         sector_consistency=sector_consistency,expected_failures=failures)
(ROOT/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status='PASS for stated assertions; expected failures retained',
    max_resolvent_error=max(r['residual'] for r in resolvent),
    min_bubble_margin=min(r['margin'] for r in bubble),expected_failures=failures),indent=2))
