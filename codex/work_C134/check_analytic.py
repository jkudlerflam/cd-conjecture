"""Independent numerical checks of identities in PROOF.md, not proofs."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from bethe_scan import solve, functions
from half_scan import bulk, constants, compact_constant, sech, root_checks, DELTAS
from check_c129 import energy

out = dict(dressed_identity_error=0., dressing_error=0., phase_error=0.,
           special_twist_energy_error=0., special_odd_energy_error=0.,
           special_bound_violation=0., cases=0, cutoff_cases=0)
records=[]
for delta in [-.9, -.5, -.1, 0.]:
    gamma=np.arccos(delta)
    p,t,dp,dt,bare,_=functions(delta)
    rho=lambda x: sech(np.pi*x/gamma)/(2*gamma)
    eps=lambda x: -np.pi*np.sin(gamma)*rho(x)
    kernel=lambda x: dt(x)/(2*np.pi)
    for x in [-2.,-.3,0.,1.7]:
        conv=quad(lambda y: kernel(x-y)*rho(y),-40,40,epsabs=1e-12)[0]
        out['dressing_error']=max(out['dressing_error'],abs(rho(x)+conv-dp(x)/(2*np.pi)))
        y=.23
        z,w=np.exp(1j*p(x)),np.exp(1j*p(y))
        ratio=(1+z*w+2*delta*z)/(1+z*w+2*delta*w)
        out['phase_error']=max(out['phase_error'],float(abs(ratio-np.exp(1j*t(x-y)))))
    for L in [4,5,16,17,64,65]:
        E,diag,roots=solve(L,L//2,delta,return_roots=True)
        sigma=lambda x: (dp(x)-np.sum(dt(x-roots))/L)/(2*np.pi)
        integral=quad(lambda x: eps(x)*sigma(x),-40,40,epsabs=2e-13,epsrel=2e-13)[0]
        direct=float(E-L*bulk(delta))
        dressed=float(np.sum(eps(roots))-L*integral)
        out['dressed_identity_error']=max(out['dressed_identity_error'],abs(direct-dressed))
        rc=root_checks(L,delta,roots)
        out['cutoff_cases']+=int('cutoff' in rc)
        records.append(dict(L=L,delta=delta,direct=direct,dressed=dressed,checks=rc))
        out['cases']+=1
for L in range(2,13):
    if L%2==0:
        val=energy(L,L//2,.5,np.pi+2*np.pi/(3*L))
        out['special_twist_energy_error']=max(out['special_twist_energy_error'],abs(val+3*L/8))
    else:
        val=energy(L,L//2,.5)
        out['special_odd_energy_error']=max(out['special_odd_energy_error'],abs(val+3*L/8))
    E,_=solve(L,L//2,.5)
    E2,_=solve(2*L,L,.5)
    out['special_bound_violation']=max(out['special_bound_violation'],float(L*(E2-2*E)-2*np.pi**2/9))
assert out['dressed_identity_error']<1e-10
assert out['dressing_error']<1e-11
assert out['phase_error']<1e-10
assert out['special_twist_energy_error']<1e-10
assert out['special_odd_energy_error']<1e-10
assert out['special_bound_violation']<1e-10
out['iterative_solver_error']=0.
for delta in DELTAS:
    E,_=solve(12,6,delta)
    Ei,_=solve(12,6,delta,iterative_threshold=0)
    out['iterative_solver_error']=max(out['iterative_solver_error'],float(abs(E-Ei)))
assert out['iterative_solver_error']<1e-10
out['point_constants']=[dict(delta=d,**constants(d)) for d in [-.999,-.99,-.9,-.5,0.]]
out['compact_constants']=[compact_constant(d) for d in [.001,.01,.1,.5,1.]]
Path('analytic_checks.json').write_text(json.dumps(out,indent=2)+'\n')
Path('analytic_detail.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(out,indent=2))
