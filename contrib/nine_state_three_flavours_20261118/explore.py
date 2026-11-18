"""Optional bounded search history; the certificate does not depend on rerunning it."""
import os,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
from physical import onsite
from mechanism import scalar,pair_coefficients,normal_order,bloch
from flow import match,beta
from audit import encode


def matrix(z):return np.diag(z[:3])+np.diag(z[3:],1)+np.diag(z[3:],-1)


def controls(z):
    mu=2*z[2]-z[1];eta=3*(z[1]-z[2]);b=12/11*(z[0]-mu-11/12*eta)
    return np.array([mu,eta,b,(12*z[3]-3*z[4])/7,(-12*z[3]+24*z[4])/7])


def main():
    assert os.environ.get('OPENBLAS_NUM_THREADS')=='1'
    started=time.perf_counter();T=scalar(3)[2];q4=normal_order(pair_coefficients(3)[4],6)
    nodes=np.array([-1.,-.3,.3,1.]);Vi=np.linalg.inv(np.polynomial.polynomial.polyvander(nodes,3))
    def coeff(z):return Vi@np.array([np.linalg.det(bloch(T/.01,matrix(z),np.arccos(x))).real for x in nodes])
    c3=coeff(np.zeros(5))[3]
    plans=[([[-.75,-.1,.65]],1500),([[-.9,-.3,.4],[-.9,0,.9],[-.8,.2,.8],[-.6,0,.6],[-.95,-.7,.5],[-.5,.7,.95],[-.95,-.2,.95],[-.7,.1,.6]],500)]
    seeds=[];attempts=0
    for targets,budget in plans:
        rng=np.random.default_rng(918)
        for roots in targets:
            target=np.polynomial.polynomial.polyfromroots(roots)*c3
            for repeat in range(12):
                attempts+=1
                result=least_squares(lambda z:coeff(z)[:3]-target[:3],rng.normal(0,.5,5),max_nfev=budget,gtol=1e-13,ftol=1e-13,xtol=1e-13)
                z=result.x;rr=np.polynomial.polynomial.polyroots(coeff(z))
                if np.max(abs(rr.imag))>1e-7 or max(abs(rr.real))>=1:continue
                r=match(T,.01*matrix(z),q4,require_generic=False);initial=np.r_[r['D'],r['c']];norm=max(abs(initial))
                def stop(t,y):return 20*norm-max(abs(y))
                stop.terminal=True
                generic=r['generic_momentum_channels'] and r['other_pair_sum_chord_margin']>1e-4
                sol=solve_ivp(lambda t,y:beta(t,y,r['a']),(0,100/norm),initial,rtol=1e-8,atol=1e-10,events=stop) if generic else None
                seeds.append(dict(target_roots=roots,repeat=repeat,max_nfev=budget,nfev=result.nfev,z=z,controls=.01*controls(z),
                    roots=r['cos_roots'],velocities=r['velocities'],D=r['D'],c=r['c'],score=min(r['D'])/np.linalg.norm(r['c']),
                    flow_status=('twentyfold_growth' if len(sol.t_events[0]) else 'finite_endpoint') if generic else 'excluded_additional_or_near_umklapp',
                    flow_end=sol.y[:,-1] if generic else None,flow_time=sol.t[-1] if generic else None,solver_success=sol.success if generic else None))
    rounding=[];selected=[]
    for denominator in (10,20,50,100,200,500,1000,2000,5000,10000,100000):
        good=[]
        for index,p in enumerate(seeds):
            numerators=np.rint(controls(p['z'])*denominator).astype(int);physical=numerators/(100*denominator)
            z=onsite(physical)/.01;rr=np.polynomial.polynomial.polyroots(Vi@np.array([np.linalg.det(bloch(T/.01,z,np.arccos(x))).real for x in nodes]))
            if np.max(abs(rr.imag))>1e-7 or max(abs(rr.real))>=.98 or min(np.diff(rr.real))<.12:continue
            r=match(T,onsite(physical),q4,require_generic=False)
            if not r['generic_momentum_channels'] or r['other_pair_sum_chord_margin']<1e-4:continue
            if np.prod(r['c'])>=0 or min(r['D'])>=0:continue
            good.append(dict(seed=index,numerators=numerators,denominator=100*denominator,roots=r['cos_roots'],D=r['D'],c=r['c'],minimum_velocity=min(r['velocities'])))
        rounding.append(dict(relative_denominator=denominator,accepted=len(good)))
        if good:selected=sorted(good,key=lambda x:-x['minimum_velocity']);break
    out=dict(attempts=attempts,accepted=len(seeds),seconds=time.perf_counter()-started,seeds=seeds,rounding=rounding,selected=selected,
        scope='Finite prescribed-root local solves and rounding checks; no exclusion theorem for the five-dimensional family.')
    (Path(__file__).parent/'exploration.json').write_text(json.dumps(out,indent=2,default=encode)+'\n')
    print(json.dumps(dict(attempts=attempts,accepted=len(seeds),seconds=out['seconds'],
        flow_counts={s:sum(p['flow_status']==s for p in seeds) for s in ('twentyfold_growth','finite_endpoint','excluded_additional_or_near_umklapp')},rounding=rounding,selected=selected),indent=2,default=encode))


if __name__=='__main__':main()
