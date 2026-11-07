"""Small independent field, kernel, constant, and crossover checks for C136."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from check_c129 import energy, bond_M
from check_c136 import DELTAS, compact_constant, all_constant
from half_scan import constants

out=dict(field_cases=0,twist_cases=0,min_field_slack=float('inf'),
         min_combined_slack=float('inf'),max_twist_ratio=0.,
         max_physical_Q=-float('inf'),max_kernel_norm_error=0.,
         endpoint_error=0.,constant_failures=0)
for d in DELTAS:
    gamma=np.arccos(d)
    ker=lambda x: np.sin(2*gamma)/(2*np.pi*(np.sinh(x)**2+np.sin(gamma)**2))
    q=2*gamma/np.pi-1
    integ=2*quad(lambda x: abs(ker(x)),0,30,epsabs=3e-13,limit=200)[0]
    out['max_kernel_norm_error']=max(out['max_kernel_norm_error'],abs(integ-q))
    for L in range(2,7):
        opt=np.array([energy(L,n,d) for n in range(L//2+1)])
        opt2=np.array([energy(2*L,n,d) for n in range(L+1)])
        for theta in [0.,.37,np.pi/2]:
            small=np.array([energy(L,n,d,theta) for n in range(L//2+1)])
            large=np.array([energy(2*L,n,d,theta) for n in range(L+1)])
            for field in [0.,(1+d)/2,15*(1+d)/16,1+d,2.]:
                n=int(np.argmin(small-field*(L/2-np.arange(len(small)))))
                EL=small[n]-field*(L/2-n)
                E2=np.min(large-field*(L-np.arange(len(large))))
                defect=E2-2*EL
                raw=large[2*n]-2*small[n]
                sector=opt2[2*n]-2*opt[n]
                M=bond_M(1,d,field)
                fs=raw-defect
                cs=sector+np.pi**2*M/(4*L)-defect
                out['min_field_slack']=min(out['min_field_slack'],float(fs))
                out['min_combined_slack']=min(out['min_combined_slack'],float(cs))
                out['max_physical_Q']=max(out['max_physical_Q'],float(L*defect/M))
                assert fs>=-2e-9 and cs>=-2e-9
                if field>=1+d:
                    assert abs(defect)<2e-9
                out['field_cases']+=1
        for n in range(1,L//2+1):
            for phi in [0.,np.pi/3,2*np.pi/3,np.pi]:
                actual=energy(L,n,d,np.pi+phi/L)-opt[n]
                exact=L*(1-np.cos(phi/L))/2
                assert -2e-9<=actual<=exact+2e-9
                out['max_twist_ratio']=max(out['max_twist_ratio'],actual/(np.pi**2/(4*L)))
                out['twist_cases']+=1

for L in range(2,13):
    for n in range(L//2+1):
        out['endpoint_error']=max(out['endpoint_error'],abs(energy(L,n,-1)+L/4))
assert out['endpoint_error']<2e-9

out['high_field_placement_cases']=0
out['largest_high_field_density']=0.
out['minimum_convex_combination_slack']=float('inf')
for d in DELTAS:
    for L in range(8,13):
        nvals=np.arange(L//2+1)
        for n in nvals:
            lower=d*L/4-(1+d)*np.sin(np.pi*n/L)/np.sin(np.pi/L)
            slack=energy(L,int(n),d)-lower
            out['minimum_convex_combination_slack']=min(out['minimum_convex_combination_slack'],slack)
            assert slack>=-2e-9
        for theta in [0.,.37,np.pi/2]:
            vals=np.array([energy(L,int(n),d,theta) for n in nvals])
            field=15*(1+d)/16
            n=int(np.argmin(vals-field*(L/2-nvals)))
            assert n/L<.25
            out['largest_high_field_density']=max(out['largest_high_field_density'],n/L)
            out['high_field_placement_cases']+=1

out['crossover_constants']=[]
for eps in [1e-2,1e-3,1e-4,1e-5,1e-6]:
    c=constants(-1+eps)
    alpha=np.pi-c['gamma']; a=np.pi*np.sin(alpha); b=np.pi/c['gamma']
    large=a*(1.5*b+c['D1']/(2*c['D'])+32*c['D']**2/b+8*c['D']*c['B'])
    out['crossover_constants'].append(dict(epsilon=eps,band_length=eps**(-.5),
        cutoff_length=c['cutoff_length'],cutoff_ratio=c['cutoff_length']*alpha**2/(4*np.pi),
        K_small_ratio=(c['cutoff_length']+1)**2*alpha**4/(16*np.pi**2),
        K_large_ratio=large*alpha**3/(32*np.pi),K_large=large))

out['example_constants']=[dict(delta=gap,eta=eta,**compact_constant(gap,eta))
                         for gap in [.1,.5,1.] for eta in [.1,.25]]
out['all_filling_constants']=[dict(delta=1+d,**all_constant(1+d)) for d in DELTAS]
out['moment_errors']=[]
for d in DELTAS:
    g=np.arccos(d)
    T=lambda x: -np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
    m=2*quad(lambda x: np.cosh(x/2)*T(x),0,30,epsabs=3e-13,limit=200)[0]
    out['moment_errors'].append(abs(m-np.sqrt(2)*np.sin(g/2-np.pi/4)))
assert max(out['moment_errors'])<2e-10
out['max_simple_constant_ratio']=0.
for gap in np.geomspace(1e-6,1,301):
    cc=all_constant(gap)
    ratio=cc['K']*gap**(25/3)/1e15
    out['max_simple_constant_ratio']=max(out['max_simple_constant_ratio'],ratio)
    assert ratio<=1

# Constant-cost free-fermion checks in the very large cutoff branch.
# No array with L entries, Bethe solve, or extensive matrix is constructed.
out['free_clipping']=[]
cc=all_constant(1.)
for L in [2**20,2**22]:
    for n in [L//2,L//2-1]:
        nu=np.longdouble(n)/L
        A=np.arccosh(cc['kappa']*L/(2*np.pi*cc['D']))/2
        pA=np.longdouble(2*np.arctan(np.tanh(A)))
        x=np.longdouble(np.pi)/L
        radius=pA/(2*x)
        lo=max(1,int(np.ceil(-radius+(np.longdouble(n)+1)/2)))
        hi=min(n,int(np.floor(radius+(np.longdouble(n)+1)/2)))
        count=max(0,hi-lo+1)
        assert count%2==n%2
        mid=(np.sin(pA)/x-np.sin(count*x)/np.sin(x)
             +(count-pA/x)*np.cos(pA))
        exact=-np.sin(np.longdouble(np.pi)*nu)*(x/6+7*x**3/360+31*x**5/15120)
        tail=exact-mid
        var=8*np.pi*np.sin(float(pA))
        bound_tail=np.pi*(128*cc['D']**2/cc['kappa']**2+72*cc['D']*cc['B']/cc['kappa'])
        assert abs(mid)*L<=cc['Vclip']/8
        assert abs(tail)*L<=bound_tail
        assert var<=cc['Vclip']
        out['free_clipping'].append(dict(L=L,n=n,cutoff=A,core_roots=count,
             scaled_midpoint=float(mid*L),scaled_tail=float(tail*L),variation=var,
             variation_bound=cc['Vclip'],tail_bound=bound_tail))
Path('audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
