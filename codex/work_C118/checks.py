"""C118 checks. Full traces, no spectral truncation. Run in this folder.

Large interacting calculations require --cluster and a SLURM allocation.
"""
from pathlib import Path
import argparse
import json
import os
import time
import numpy as np
from scipy.special import logsumexp
from scipy.linalg import eigh
from scipy.integrate import quad
from ring_tools import spectrum

ROOT = Path(__file__).resolve().parent


def xx_logz(n, beta):
    # Long double limits cancellation of the bulk term in large-ring ratios.
    answers = []
    for shift, parity in [(0.5, 0), (0.0, 1)]:
        eps = np.cos(2*np.longdouble(np.pi)*(np.arange(n, dtype=np.longdouble)+shift)/n)
        even, odd = np.longdouble(0), np.longdouble(-np.inf)
        for energy in eps:
            even, odd = np.logaddexp(even, odd-beta*energy), np.logaddexp(odd, even-beta*energy)
        answers.append((even, odd)[parity])
    return np.logaddexp(*answers)


def stats(ev, beta):
    shifted = ev-ev[0]
    lw = logsumexp(-beta*shifted)
    weights = np.exp(-beta*shifted-lw)
    mean = np.dot(weights, shifted)
    variance = np.dot(weights, (shifted-mean)**2)
    return float(lw-beta*ev[0]), float(lw), float(variance)


def local_bond(model):
    if model == 'SU3':
        return np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9), 3
    x = np.array([[0.,1.],[1.,0.]])/2
    y = np.array([[0.,-1j],[1j,0.]])/2
    z = np.diag([.5,-.5])
    h = np.kron(x,x)+np.kron(y,y)
    if model == 'Heisenberg':
        h += np.kron(z,z)
    return h.real, 2


def embed(h, d, n, edges):
    powers = d**np.arange(n, dtype=np.int64)
    result = np.zeros((d**n, d**n))
    for state in range(d**n):
        for i,j,weight in edges:
            a,b = (state//powers[i])%d, (state//powers[j])%d
            for aa in range(d):
                for bb in range(d):
                    value = weight*h[d*aa+bb,d*a+b]
                    if value:
                        target = state+(aa-a)*powers[i]+(bb-b)*powers[j]
                        result[target,state] += value
    return result


def entropy(rho):
    ev = np.linalg.eigvalsh((rho+rho.T)/2)
    ev = ev[ev>1e-15]
    return float(-np.dot(ev,np.log(ev)))


def seam_checks(model, length):
    h,d = local_bond(model)
    edges0 = [(a*length+i, a*length+(i+1)%length, 1.) for a in range(2) for i in range(length)]
    h0 = embed(h,d,2*length,edges0)
    vs=[]
    for i in range(length):
        j=(i+1)%length
        edges=[(i,length+j,1.), (length+i,j,1.), (i,j,-1.), (length+i,length+j,-1.)]
        vs.append(embed(h,d,2*length,edges))
    v=vs[-1]
    h1=h0+v
    vbar=sum(vs)/length
    e0,u0=eigh(h0)
    e1,u1=eigh(h1)
    ea=np.linalg.eigvalsh(h0+vbar)
    vnorm=np.max(np.abs(np.linalg.eigvalsh(vbar)))
    rows=[]
    for beta in [.25,1.,4.,16.]:
        z0=logsumexp(-beta*e0)
        z1=logsumexp(-beta*e1)
        za=logsumexp(-beta*ea)
        sigma=(u0*np.exp(-beta*e0-z0))@u0.T
        rho=(u1*np.exp(-beta*e1-z1))@u1.T
        r=float(z0-z1)
        v0=float(np.sum(sigma*v.T))
        v1=float(np.sum(rho*v.T))
        # Relative entropies evaluated using exact Gibbs logarithms.
        dsr=float(-entropy(sigma)+beta*np.sum(sigma*h1.T)+z1)
        drs=float(-entropy(rho)+beta*np.sum(rho*h0.T)+z0)
        n=d**length
        a=np.trace(rho.reshape(n,n,n,n),axis1=1,axis2=3)
        b=np.trace(rho.reshape(n,n,n,n),axis1=0,axis2=2)
        mi=entropy(a)+entropy(b)-entropy(rho)
        ha=embed(h,d,length,[(i,(i+1)%length,1.) for i in range(length)])
        zl=z0/2
        da=-entropy(a)+beta*np.sum(a*ha.T)+zl
        db=-entropy(b)+beta*np.sum(b*ha.T)+zl
        err=max(abs(r-beta*v0+dsr),abs(r-beta*v1-drs),abs(drs-mi-da-db))
        assert err<2e-10
        assert za<=z1+2e-10
        assert dsr>=-1e-10 and drs>=-1e-10
        rows.append(dict(beta=beta,R=r,D_split_joined=dsr,D_joined_split=drs,
                         MI=mi,seam_split=beta*v0,seam_joined=beta*v1,
                         marginal_D=da+db,identity_error=err,
                         Jensen_ceiling=float(z0-za)))
    return dict(model=model,L=length,averaged_seam_norm=float(vnorm),
                averaged_energy_cost=float(ea[0]-e0[0]),rows=rows)


def constant_swap_check(model):
    h,d=local_bond(model)
    htwo=embed(h,d,4,[(0,1,1.),(2,3,1.)])
    swap=np.eye(d*d).reshape(d,d,d,d).transpose(1,0,2,3).reshape(d*d,d*d)
    gen=embed(swap,d,4,[(0,2,1.),(1,3,1.)])
    comm=htwo@gen-gen@htwo
    commnorm=float(np.linalg.norm(comm,2))
    assert commnorm>1e-3
    # The same angle at both sites has zero spatial gradient.
    es,us=eigh(gen)
    unitary=(us*np.exp(1j*np.pi*es/4))@us.T
    mismatch=unitary.conj().T@htwo@unitary-htwo
    norm=float(np.linalg.norm(mismatch,2))
    assert norm>1e-3
    return dict(model=model,commutator_norm=commnorm,constant_angle_mismatch=norm)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cluster',action='store_true')
    args=parser.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Large ED requires a SLURM allocation.')
    start=time.time()
    out=dict(mode='cluster' if args.cluster else 'local',models={},xx=[],seams=[])
    sizes={'XX':[4,8],'Heisenberg':[4,8],'SU3':[3,6]}
    if args.cluster:
        sizes={'XX':[4,6,8,12],'Heisenberg':[4,6,8,12],'SU3':[3,6,12]}
    for model, ns in sizes.items():
        spec={n:spectrum(model,n)[0] for n in ns}
        rows=[]
        for n in ns:
            if 2*n not in spec:
                continue
            e,e2=spec[n],spec[2*n]
            delta=float(e2[0]-2*e[0])
            m=.5 if model!='SU3' else 1.
            for tau in np.geomspace(.025,32,121):
                beta=float(tau*n)
                z,w,var=stats(e,beta)
                z2,w2,var2=stats(e2,beta)
                r=2*z-z2
                assert abs(r)<=4*beta*m+1e-8
                d=3 if model=='SU3' else 2
                # A sampled extension, not a theorem for interacting chains.
                trial=2*np.log(d)+(10*np.pi/3)*m*np.log2(d)*beta/n
                assert r<=trial+1e-8
                rows.append(dict(L=n,tau=float(tau),R=r,curvature=2*var-var2,
                                 thermal_remainder=2*w-w2,
                                 delta=delta,scaled_defect=n*delta/m,
                                 trial_ceiling=float(trial),trial_slack=float(trial-r)))
        out['models'][model]=rows
    err=0.
    for n in [4,6,8]:
        ev=spectrum('XX',n)[0]
        for beta in [.1,1.,4.,16.]:
            err=max(err,abs(float(xx_logz(n,beta))-stats(ev,beta)[0]))
    assert err<1e-10
    out['xx_ED_error']=err
    out['xx_gaps']=[]
    for n in [4,8]:
        ev=spectrum('XX',n)[0]
        gap=float(ev[np.flatnonzero(ev>ev[0]+1e-10)[0]]-ev[0])
        mult=int(np.count_nonzero(np.abs(ev-ev[0]-gap)<1e-10))
        out['xx_gaps'].append(dict(L=n,gap=gap,multiplicity=mult))
    for n in [4,6,8,16,64,256,1024,4096]:
        for tau in [.01,.05,.1,.25,.5,1.,2.,4.,8.,32.]:
            beta=np.longdouble(tau*n)
            r=float(2*xx_logz(n,beta)-xx_logz(2*n,beta))
            ceiling=2*np.log(2)+5*np.pi*float(beta)/(3*n)
            assert r<=ceiling+2e-7
            out['xx'].append(dict(L=n,tau=tau,R=r,proved_ceiling=ceiling))
    # Check the quadrature precursor independently against integration.
    out['quadrature']=[]
    for n in [4,6,16,64]:
        for beta in [.1,1.,8.,64.]:
            pressure=quad(lambda k:np.logaddexp(0.,-beta*np.cos(k)),0,2*np.pi,
                          epsabs=1e-10)[0]/(2*np.pi)
            xi=float(xx_logz(n,beta))-n*pressure
            c=2*np.pi*beta/(3*n)
            assert -c-1e-8<=xi<=np.log(2)+c+1e-8
            out['quadrature'].append(dict(L=n,beta=beta,Xi=xi,lower=-c,upper=np.log(2)+c))
    # Small dense matrices only. SU3 L=3 uses dimension 729.
    for model,n in [('XX',4),('Heisenberg',4),('SU3',3)]:
        out['seams'].append(seam_checks(model,n))
    out['constant_swaps']=[constant_swap_check(m) for m in ['XX','Heisenberg','SU3']]
    # Classical nonnegative transfer comparison, oriented and nonsymmetric.
    h=np.array([[.3,-.4,.8],[.2,.1,-.7],[-.5,.6,-.2]])
    classical=[]
    for beta in [.1,1.,4.,16.]:
        t=np.exp(-beta*(h-h.min()))
        for n in [3,4,8,16]:
            a=np.linalg.matrix_power(t,n)
            r=float(2*np.log(np.trace(a))-np.log(np.trace(a@a)))
            assert r<=np.log(3)+1e-12
            classical.append(dict(beta=beta,L=n,R=r))
    out['classical']=classical
    out['seconds']=time.time()-start
    filename='checks_cluster.json' if args.cluster else 'checks_local.json'
    (ROOT/filename).write_text(json.dumps(out,indent=2)+'\n')
    summary={model:dict(min_curvature=min(x['curvature'] for x in rows),
                        max_curvature=max(x['curvature'] for x in rows),
                        max_thermal_remainder=max(x['thermal_remainder'] for x in rows))
             for model,rows in out['models'].items()}
    print(json.dumps(dict(file=filename,seconds=out['seconds'],xx_ED_error=err,summary=summary,
                         seams=out['seams']),indent=2))


if __name__=='__main__':
    main()
