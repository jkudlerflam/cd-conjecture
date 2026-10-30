"""Independent small checks of the new bounds and SDP calibrations."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh, svdvals
from scipy.special import gammaln
from ring_tools import embed, local_bond
from sdp import MinSDP

ROOT=Path(__file__).resolve().parent


def main():
    rng=np.random.default_rng(122)
    out=dict(diagonal=[],order=[],majorant=[],tail=[])
    # Classical diagonal W has mu=sum_D max_C W(C,D), independent of optimization.
    for d,n,seam in [(2,4,[0,3]),(3,3,[0,2])]:
        opt=MinSDP(d,n,seam)
        diagonal=.1+rng.random(d**n)
        target=float(diagonal[opt.perm].reshape(opt.q,opt.m).max(axis=0).sum())
        sol=opt.solve(np.diag(diagonal))
        assert sol['lower']-1e-8<=target<=sol['upper']+1e-8
        out['diagonal'].append(dict(d=d,n=n,exact=target,SDP=sol))
    for model in ['XX','Heisenberg','SU3']:
        h,d=local_bond(model)
        L=4 if d==2 else 3
        h0=embed(h,d,2*L,[(a*L+i,a*L+(i+1)%L,1.) for a in range(2) for i in range(L)])
        h1=embed(h,d,2*L,[(i,(i+1)%(2*L),1.) for i in range(2*L)])
        e0,u0=eigh(h0); e1,u1=eigh(h1)
        overlap=u1.T@u0
        M=.5 if d==2 else 1.
        for beta in [.01,.1,.5,1.]:
            f=overlap*np.exp(beta*(e1[:,None]-e0[None,:])/2)
            measured=float(2*np.log(svdvals(f)[0]))
            b=beta*M
            ceiling=12*b+4*np.expm1(2*b)
            assert measured<=ceiling+1e-10
            if model=='XX':
                lower=2*beta/(L*np.sin(np.pi/L))-2*np.log(2)
                assert measured>=lower-1e-10
            out['order'].append(dict(model=model,beta=beta,log_c=measured,ceiling=float(ceiling)))
    # Check the positive coefficient recurrence against its closed generating sum.
    v,J,s=4.,2.,2
    for t in [.01,.1,.3]:
        k=np.arange(100)
        a=2*np.expm1(2*J*t)
        terms=np.exp((v+2*J*s)*t+k*np.log(a)-gammaln(k+1))
        total=np.exp((v+2*J*s)*t+a)
        residual=abs(terms.sum()/total-1)
        assert residual<1e-12
        derivative=terms*((v+2*J*s)+k*(2*J*np.exp(2*J*t)/np.expm1(2*J*t)))
        rhs=(v+2*J*(s+k))*terms+4*J*np.r_[0.,terms[:-1]]
        err=float(np.max(np.abs(derivative-rhs))/max(derivative))
        assert err<1e-12
        out['majorant'].append(dict(t=t,sum_relative_error=float(residual),recurrence_error=err))
    # Width inequality used after multiplying by q=d^(4+4r).
    for d in [2,3]:
        for b in [.01,.1]:
            a=2*np.expm1(2*b)
            r=int(np.ceil(2*np.e*a*d**4))+1
            k=np.arange(r,r+500)
            logtail=np.logaddexp.reduce(k*np.log(a)-gammaln(k+1))
            log_actual=(4+4*r)*np.log(d)+np.log(2)+6*b+logtail
            log_bound=np.log(4)+4*np.log(d)+6*b-r*np.log(2)
            assert log_actual<=log_bound+1e-10
            out['tail'].append(dict(d=d,b=b,r=r,log_q_tail=float(log_actual),log_bound=float(log_bound)))
    (ROOT/'analytic_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Passed diagonal SDP, full-order, recurrence, and width checks.')


if __name__=='__main__':
    main()
