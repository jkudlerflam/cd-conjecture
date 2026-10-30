"""C122 small full Gibbs SDPs. All computations and outputs stay in this folder."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import eigh, eigvalsh, svdvals
from scipy.special import logsumexp
from ring_tools import local_bond, embed
from sdp import MinSDP

ROOT = Path(__file__).resolve().parent


def run(model,L,seam=None,taus=(0.,.1,1.,4.,16.)):
    h,d = local_bond(model)
    n = 2*L
    if seam is None:
        seam = [0,L-1,L,n-1]
    sdp = MinSDP(d,n,seam)
    h0 = embed(h,d,n,[(a*L+i,a*L+(i+1)%L,1.) for a in range(2) for i in range(L)])
    h1 = embed(h,d,n,[(i,(i+1)%n,1.) for i in range(n)])
    v = h1-h0
    es,us = zip(eigh(h0),eigh(h1))
    ground = []
    for e,u in zip(es,us):
        indices = np.flatnonzero(e-e[0]<1e-9)
        p = u[:,indices]@u[:,indices].T
        sol = sdp.solve(p)
        if len(indices)==1:
            psi = u[sdp.perm,0].reshape(sdp.q,sdp.m)
            exact = float(svdvals(psi).sum()**2)
            assert sol['lower']-1e-7 <= exact <= sol['upper']+1e-7
            sol['pure_exact'] = exact
        sol['degeneracy'] = len(indices)
        ground.append(sol)
    rows = []
    for tau in taus:
        beta = float(L*tau)
        sols,logzs,logs,wmat = [],[],[],[]
        for e,u in zip(es,us):
            weights = np.exp(-beta*(e-e[0]))
            w = (u*weights)@u.T
            sol = sdp.solve(w)
            if beta == 0:
                assert sol['lower']-1e-7 <= sdp.m <= sol['upper']+1e-7
            sols.append(sol)
            logzs.append(float(-beta*e[0]+logsumexp(-beta*(e-e[0]))))
            logs.append(float(-beta*e[0]+np.log(sol['value'])))
            wmat.append(w)
        delta = es[1][0]-es[0][0]
        r = logzs[0]-logzs[1]
        ratio = logs[0]-logs[1]
        lo = beta*delta+np.log(sols[0]['lower']/sols[1]['upper'])
        hi = beta*delta+np.log(sols[0]['upper']/sols[1]['lower'])
        entropy_difference = logzs[0]-logs[0]-logzs[1]+logs[1]
        assert abs(r-ratio-entropy_difference)<1e-10
        assert abs(r-ratio)<=2*np.log(sdp.q)+1e-8
        vs = float(np.sum(wmat[0]*v)/np.trace(wmat[0]))
        row=dict(tau=tau,beta=beta,R=r,log_mu_ratio=ratio,
                 ratio_interval=[float(lo),float(hi)],log_mu_Y=logs[0],log_mu_X=logs[1],
                 energy_part=float(beta*delta),entropy_difference=entropy_difference,
                 centered_ratio=float(ratio-beta*delta),SDP_Y=sols[0],SDP_X=sols[1],
                 order_log_lower_bound=beta*vs)
        if beta*max(np.ptp(e) for e in es)<=25:
            # F = X^(-1/2) Y^(1/2), so Y <= ||F||^2 X.
            overlap = us[1].T@us[0]
            weighted = overlap*np.exp(beta*(es[1][:,None]-es[0][None,:])/2)
            logc = float(2*np.log(svdvals(weighted)[0]))
            assert logc+1e-8 >= beta*vs
            assert ratio<=logc+3e-5
            row['full_order_log_c'] = logc
        rows.append(row)
        print(model,L,len(seam),tau,round(ratio,8),flush=True)
    return dict(model=model,L=L,d=d,seam=seam,dimension=d**n,
                complement_dimension=sdp.m,SDP_variables=len(sdp.cost),
                energy_defect=float(es[1][0]-es[0][0]),ground=ground,rows=rows)


def main():
    start = time.time()
    results=[]
    jobs=[('XX',4,None),('Heisenberg',4,None),('SU3',3,None),
          ('XX',2,None),('Heisenberg',2,None),
          ('XX',4,[0,1,3,4,5,7]),('Heisenberg',4,[0,1,3,4,5,7])]
    for model,L,seam in jobs:
        results.append(run(model,L,seam))
        (ROOT/'checks.json').write_text(json.dumps(dict(results=results,seconds=time.time()-start),indent=2)+'\n')
    print('Completed',len(results),'cases in',time.time()-start,'seconds',flush=True)


if __name__=='__main__':
    main()
