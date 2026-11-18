"""Conditional coupled-flow comparison, with no fitted p153c inputs.

The constants C and D and all example dimensions are illustrative inputs
to a differential-inequality theorem, not microscopic estimates.
"""
from medium import np, threadpool_limits
from scipy.integrate import solve_ivp
from pathlib import Path
import json


def budget(a0,q_lower,x_lower,C,D,rho=0.1):
    """Sufficient total initial coupling for Y(l)<=rho until detuning.

    Assumptions: Y' <= (2-min x)Y+C Y^2, |x_i'|<=D Y^2,
    a maintained physical mismatch >=q_lower for every retained channel,
    and a controlled omitted-channel budget. See RG.md for scope.
    """
    assert a0>0 and q_lower>0 and C>=0 and D>=0 and rho>0
    stop=max(0.0,float(np.log(1/(a0*q_lower))))
    growth=max(0.0,2-x_lower+D*rho*rho*stop)
    expm1=float(np.expm1(growth*stop))
    phi=stop if growth==0 else expm1/growth
    ycrit=rho/(1+expm1+C*rho*phi)
    return dict(a0=a0,q_lower=q_lower,x_lower=x_lower,C=C,D=D,rho=rho,
                l_detuning=stop,dimension_drift_cap=D*rho*rho*stop,
                comparison_growth=growth,Y0_max=ycrit,
                lambda4_over_t_per_velocity_max=ycrit/(4*np.pi*a0*a0))


def check_comparison(row):
    """Integrate a saturating coupled system, not a p153c RG calculation."""
    C,D=row['C'],row['D']
    def ode(t,v):
        Y,x=v
        return [(2-x)*Y+C*Y*Y,-D*Y*Y]
    end=row['l_detuning']; y0=row['Y0_max']
    sol=solve_ivp(ode,(0,end),[y0,row['x_lower']],rtol=2e-12,atol=2e-14,
                  dense_output=True,max_step=max(end/100,1e-3))
    samples=sol.sol(np.linspace(0,end,201))
    assert sol.success
    assert samples[0].max()<=row['rho']+1e-10
    drift=row['x_lower']-samples[1].min()
    assert drift<=row['dimension_drift_cap']+1e-10
    return dict(final_Y=float(samples[0,-1]),max_Y=float(samples[0].max()),
                actual_dimension_drift=float(drift))


def main():
    with threadpool_limits(limits=1):
        rows=[]
        for a0 in (4,8):
            for x in (1.2,1.5,1.8,2.0):
                row=budget(a0,.01,x,1,1)
                row['comparison_check']=check_comparison(row)
                rows.append(row)
        zero=budget(4,.01,2,1,0)
        exact=.1/(1+.1*zero['l_detuning'])
        assert abs(zero['Y0_max']-exact)<1e-15
        decaying=budget(4,.01,2.5,0,0)
        assert abs(decaying['Y0_max']-.1)<1e-15
    out=dict(status='Conditional calibration only; actual p153c inputs are unknown.',
             convention='H4=lambda4 V4+h.c.; S4=y4/(2pi a0^2) integral dx d(u tau) cos; y4=4pi a0^2 abs(lambda4)/u.',
             velocity='u/t in physical sites per effective time',
             actual_inputs=dict(x_lower=None,q_lower=None,a0=None,velocity=None,
                 lambda4=None,other_vertex_budget=None,C=None,D=None),
             examples=rows,
             zero_growth_check=dict(calculated=zero['Y0_max'],analytic=exact))
    Path(__file__).with_name('flow.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in rows:
        print(row['a0'],row['x_lower'],f"Y0 <= {row['Y0_max']:.9g}",
              f"abs(lambda4)/u <= {row['lambda4_over_t_per_velocity_max']:.9g}")


if __name__=='__main__': main()
