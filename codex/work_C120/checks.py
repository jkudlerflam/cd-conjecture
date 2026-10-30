"""C120. Small full ED, exact trace identities, and finite Trotter histories.

Run with single-threaded BLAS. No network access or large ED.
"""
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import eigh, eigvalsh, svdvals
from scipy.special import logsumexp
from ring_tools import local_bond, embed, xx_logz

ROOT = Path(__file__).resolve().parent


def gibbs(ev, u, beta):
    logw = -beta*(ev-ev[0])
    logw -= logsumexp(logw)
    return logw, (u*np.exp(logw))@u.T, float(logsumexp(-beta*ev))


def maxdiv(logp, logq, overlap):
    # Largest singular value of sigma^(-1/2) rho^(1/2).
    a = (logp[:, None]-logq[None, :])/2
    shift = float(np.max(a))
    matrix = np.exp(a-shift)*overlap
    return float(2*(np.log(svdvals(matrix, check_finite=False)[0])+shift))


def petz(logp, logq, overlap2, alpha):
    with np.errstate(divide='ignore'):
        lq = logsumexp(alpha*logp[:, None]+(1-alpha)*logq[None, :]
                       + np.log(overlap2))
    return float(lq/(alpha-1))


def sandwich(logp, logq, overlap, alpha):
    a = (logp[:, None]+((1-alpha)/alpha)*logq[None, :])/2
    shift = float(np.max(a))
    singular = svdvals(np.exp(a-shift)*overlap, check_finite=False)
    with np.errstate(divide='ignore'):
        value = logsumexp(2*alpha*(np.log(singular)+shift))
    return float(value/(alpha-1))


def partial_site(rho, d):
    n = rho.shape[0]//d
    r = rho.reshape(n,d,n,d)
    return np.trace(r,axis1=0,axis2=2), np.trace(r,axis1=1,axis2=3)


def min_witness(rho,d):
    # The fixed marginal is feasible for the dimension estimate, not necessarily optimal.
    _, env = partial_site(rho,d)
    ev,u = eigh(env)
    assert ev[0]>0
    inv = (u*(ev**-.5))@u.T
    inv = np.kron(inv,np.eye(d))
    c = inv@rho@inv
    return float(np.log(eigvalsh(c,subset_by_index=[len(c)-1,len(c)-1])[0]))


def ed_checks(model, L):
    h,d=local_bond(model)
    edges0=[(a*L+i,a*L+(i+1)%L,1.) for a in range(2) for i in range(L)]
    h0=embed(h,d,2*L,edges0)
    v=embed(h,d,2*L,[(L-1,L,1.),(2*L-1,0,1.),(L-1,0,-1.),(2*L-1,L,-1.)])
    h1=h0+v
    e0,u0=eigh(h0); e1,u1=eigh(h1)
    small=embed(h,d,L,[(i,(i+1)%L,1.) for i in range(L)])
    es=eigvalsh(small)
    assert np.max(np.abs(np.sort((es[:,None]+es[None,:]).ravel())-e0))<1e-11
    overlap=u1.T@u0
    overlap2=overlap**2
    rows=[]
    for beta in [0.,.1,.5,1.,2.,4.,8.,16.,32.]:
        lp,rho,z1=gibbs(e1,u1,beta); lq,sigma,z0=gibbs(e0,u0,beta)
        r=z0-z1
        dsr=float(np.dot(np.exp(lq),lq)+beta*np.sum(sigma*h1)+z1)
        drs=float(np.dot(np.exp(lp),lp)+beta*np.sum(rho*h0)+z0)
        vr=float(np.sum(rho*v)); vs=float(np.sum(sigma*v))
        assert abs(r-beta*vr-drs)<1e-10
        assert abs(r-beta*vs+dsr)<1e-10
        # Root fidelity from singular values, no square root of a nearly singular density matrix.
        rootf=float(np.sum(svdvals(np.exp((lp[:,None]+lq[None,:])/2)*overlap,
                                      check_finite=False)))
        td=float(np.sum(np.abs(eigvalsh(rho-sigma)))/2)
        assert 1-rootf<=td+1e-10
        assert td<=np.sqrt(max(0.,1-rootf**2))+2e-7
        row=dict(beta=beta,R=r,relative_forward=drs,relative_reverse=dsr,
                 seam_forward=beta*vr,seam_reverse=beta*vs,root_fidelity=rootf,
                 trace_distance=td,fidelity_divergence=-2*np.log(rootf),support_D0=0.,
                 cap_one_cut=2*np.log(d),cap_four_sites=4*np.log(d),
                 fidelity_mixed_term=r+2*np.log(rootf),
                 logZ_ratio_identity_error=abs(r-beta*vr-drs),renyi=[])
        # Exact powers of Gibbs spectra are used in Petz traces, including very small eigenvalues.
        for alpha in [.05,.1,.25,.5,.75,.9]:
            da=petz(lp,lq,overlap2,alpha)
            rev=petz(lq,lp,overlap2.T,alpha)
            mixed=logsumexp(-alpha*beta*e1[:,None]-(1-alpha)*beta*e0[None,:]
                               +np.log(np.maximum(overlap2,np.finfo(float).tiny)))
            identity=da+(mixed-z1)/(1-alpha)
            assert abs(identity-r)<2e-9
            assert beta*vr+da<=r+1e-9
            assert r<=beta*vs-rev+1e-9
            entry=dict(alpha=alpha,Petz=da,Petz_reverse=rev,
                       mixed_term=float((mixed-z1)/(1-alpha)),identity_error=abs(identity-r),
                       forward_replacement_slack=beta*vr+da-r,
                       reverse_bound_slack=beta*vs-rev-r)
            # Low orders amplify singular-value roundoff. Restrict the numerical condition number.
            if beta*max(np.ptp(e0),np.ptp(e1))*max(1.,(1-alpha)/alpha)<=25:
                sa=sandwich(lp,lq,overlap,alpha)
                assert sa<=drs+1e-8
                assert sa<=da+1e-8
                entry['sandwiched']=sa
                if alpha==.5:
                    assert abs(sa+2*np.log(rootf))<1e-9
            row['renyi'].append(entry)
        if beta*max(np.ptp(e0),np.ptp(e1))<=25:
            dm=maxdiv(lp,lq,overlap)
            dmr=maxdiv(lq,lp,overlap.T)
            assert dm>=drs-1e-9 and dmr>=dsr-1e-9
            row.update(Dmax=dm,Dmax_reverse=dmr,hilbert_metric=dm+dmr,
                       max_upper_slack=beta*vr+dm-r,
                       unnormalized_log_order=dm-r)
        if beta<=2.:
            witness=min_witness(rho,d)
            assert -np.log(d)-1e-8<=witness<=np.log(d)+1e-8
            row['conditional_min_fixed_reference_Dmax']=witness
        rows.append(row)
    # Unique ground states, where conditional min-entropy has a closed expression.
    assert e0[1]-e0[0]>1e-8 and e1[1]-e1[0]>1e-8
    ground=[]
    for u in [u0,u1]:
        c,_=partial_site(np.outer(u[:,0],u[:,0]),d)
        hm=-2*np.log(np.sum(np.sqrt(np.maximum(eigvalsh(c),0))))
        assert -np.log(d)-1e-9<=hm<=np.log(d)+1e-9
        ground.append(float(hm))
    return dict(model=model,L=L,d=d,dimension=d**(2*L),
                energy_defect=float(e1[0]-e0[0]),rows=rows,
                ground_conditional_min=ground,
                ground_root_fidelity=float(abs(np.dot(u0[:,0],u1[:,0]))))


def transfer(model,L,period,n,beta):
    h,d=local_bond(model)
    ev,u=eigh(h)
    # Subtract a scalar to keep local gates bounded, which leaves the ratio invariant.
    h=h-ev[0]*np.eye(d*d)
    gate=(u*np.exp(-beta*(ev-ev[0])/n))@u.T
    depth=period*n
    D=d**depth
    hist=(np.arange(D)[:,None]//(d**np.arange(depth)))%d
    ws=[]
    for color in range(period):
        w=np.ones((D,D))
        for t in range(depth):
            nxt=(t+1)%depth
            if t%period==color:
                out=d*hist[:,nxt,None]+hist[None,:,nxt]
                ins=d*hist[:,t,None]+hist[None,:,t]
                w*=gate[out,ins]
            elif t%period!=(color-1)%period:
                w*=(hist[:,nxt]==hist[:,t])[:,None]
        ws.append(w)
    t=np.eye(D)
    for w in ws: t=t@w
    b=np.linalg.matrix_power(t,L//period)
    z=float(np.trace(b)); z2=float(np.trace(b@b))
    assert z>0 and z2>0
    # Independent temporal contraction on the smaller physical ring.
    product=np.eye(d**L)
    for color in range(period):
        layer=embed(h,d,L,[(i,(i+1)%L,1.) for i in range(color,L,period)])
        e,u=eigh(layer)
        product=product@((u*np.exp(-beta*e/n))@u.T)
    ztime=float(np.trace(np.linalg.matrix_power(product,n)))
    assert abs(ztime-z)<2e-9*max(1.,z)
    product2=np.eye(d**(2*L))
    for color in range(period):
        layer=embed(h,d,2*L,[(i,(i+1)%(2*L),1.) for i in range(color,2*L,period)])
        e,u=eigh(layer)
        product2=product2@((u*np.exp(-beta*e/n))@u.T)
    ztime2=float(np.trace(np.linalg.matrix_power(product2,n)))
    assert abs(ztime2-z2)<2e-9*max(1.,z2)
    neg=float(np.sum(np.maximum(-(b*b.T),0.)))
    penalty=float(np.log1p(neg/z2))
    r=2*np.log(z)-np.log(z2)
    assert r<=np.log(D)+penalty+1e-9
    exact=eigvalsh(embed(h,d,L,[(i,(i+1)%L,1.) for i in range(L)]))
    doubled=eigvalsh(embed(h,d,2*L,[(i,(i+1)%(2*L),1.) for i in range(2*L)]))
    rexact=float(2*logsumexp(-beta*exact)-logsumexp(-beta*doubled))
    return dict(model=model,L=L,period=period,slices=n,beta=beta,history_dimension=D,
                R_trotter=float(r),R_ED=rexact,temporal_residual=abs(ztime-z),
                temporal_doubled_residual=abs(ztime2-z2),
                log_history_dimension=float(np.log(D)),negative_pair_penalty=penalty,
                negative_entries=int(np.count_nonzero(t < -1e-14)),
                negative_pairs=int(np.count_nonzero(b*b.T < -1e-14)),
                signed_bound_slack=float(np.log(D)+penalty-r),
                seam_dimension_bound_slack=float(2*np.log(d)-r))


def channel_examples():
    rows=[]
    for m in [1,2,4,8]:
        angle=np.pi/(2*m)
        rot=np.array([[np.cos(angle),-np.sin(angle),0.],
                      [np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])
        for a in [.9,.99,.9999]:
            t=np.zeros((4,4));t[0,0]=1;t[1:,1:]=a*rot
            b=np.linalg.matrix_power(t,m)
            z=np.trace(b);z2=np.trace(b@b)
            r=2*np.log(z)-np.log(z2)
            formula=np.log((1+a**m)/(1-a**m))
            assert abs(r-formula)<1e-10
            if 3*a**(2*m)<1:
                ceiling=2*np.log1p(3*a**m)-np.log1p(-3*a**(2*m))
                assert r<=ceiling+1e-10
            rows.append(dict(power=m,contraction=a,ratio_log=float(r),log_super_dimension=np.log(4),
                             projective_diameter=float(2*np.log((1+a)/(1-a)))))
    return rows


def qubit_examples():
    rows=[]
    for t in [1.,4.,16.,64.]:
        lq=np.array([-t,0.])-np.logaddexp(0.,-t)
        lp=lq[::-1]
        ds=[dict(alpha=a,Petz=petz(lp,lq,np.eye(2),a),
                 sandwiched=sandwich(lp,lq,np.eye(2),a)) for a in [.05,.25,.5,.75,.9]]
        rows.append(dict(t=t,Dmax=t,support_D0=0.,root_fidelity=float(1/np.cosh(t/2)),
                         trace_distance=float(np.tanh(t/2)),renyi=ds))
    return rows


def xx_large():
    rows=[]
    def logw(n,beta):
        answers=[]
        for shift,parity in [(.5,0),(0.,1)]:
            eps=np.cos(2*np.pi*(np.arange(n)+shift)/n)
            eps[np.abs(eps)<1e-13]=0.
            target=parity^(int(np.count_nonzero(eps<0))%2)
            even,odd=0.,-np.inf
            for energy in np.abs(eps):
                even,odd=np.logaddexp(even,odd-beta*energy),np.logaddexp(odd,even-beta*energy)
            good=(shift==.5 and n%4==0) or (shift==0. and n%4==2)
            answers.append((even,odd)[target]-(0. if good else beta*np.tan(np.pi/(2*n))))
        return float(np.logaddexp(*answers))
    calibration=0.
    for L in [4,8,32,128,512,2048,8192]:
        for tau in [.1,1.,4.,16.,64.]:
            beta=L*tau
            angle=np.pi/L
            defect=2*np.sin(angle/4)**2/(np.sin(angle/2)*np.cos(angle/2))
            r=float(beta*defect+2*logw(L,beta)-logw(2*L,beta))
            if L<=32:
                calibration=max(calibration,abs(r-float(2*xx_logz(L,beta)-xx_logz(2*L,beta))))
            cap=4*np.log(2)
            ceiling=2*np.log(2)+5*np.pi*tau/3
            assert r<=ceiling+2e-6
            # Any normalized-state constant-capacity substitution is <= cap (or 1 for distance).
            rows.append(dict(L=L,tau=tau,R=r,cap_four_sites=cap,cap_slack=cap-r,
                             distance_ceiling_slack=1-r,proved_XX_slack=ceiling-r))
    assert calibration<2e-9
    return rows,calibration


def xx_ground():
    rows=[]
    for L in [4,8,16,32,64,128,256,512]:
        def occupied(n):
            shift=.5 if (n//2)%2==0 else 0.
            k=2*np.pi*(np.arange(n)+shift)/n
            k=k[np.cos(k)<-1e-13]
            return np.exp(1j*np.arange(n)[:,None]*k[None,:])/np.sqrt(n)
        little=occupied(L)
        split=np.zeros((2*L,L),complex)
        split[:L,:L//2]=little;split[L:,L//2:]=little
        joined=occupied(2*L)
        phase,logf=np.linalg.slogdet(joined.conj().T@split)
        assert abs(abs(phase)-1)<1e-10
        f=float(np.exp(logf))
        conditional=[]
        for modes in [split,joined]:
            corr=modes[:4]@modes[:4].conj().T
            nu=np.clip(eigvalsh(corr),0,1)
            nu[nu<1e-12]=0.;nu[nu>1-1e-12]=1.
            conditional.append(float(-2*np.sum(np.log(np.sqrt(nu)+np.sqrt(1-nu)))))
        rows.append(dict(L=L,ground_root_fidelity=f,ground_trace_distance=float(np.sqrt(1-f*f)),
                         support_D0=-2*float(logf),Dmax='infinite for distinct pure states',
                         ground_conditional_min_four_sites=conditional,
                         renyi=[dict(alpha=a,Petz=-2*float(logf)/(1-a),
                                     sandwiched=-2*a*float(logf)/(1-a)) for a in [.05,.25,.5,.75,.9]]))
    return rows


def main():
    start=time.time()
    large,calibration=xx_large()
    out=dict(ed=[ed_checks(m,L) for m,L in [('XX',4),('Heisenberg',4),('SU3',3)]],
             channels=channel_examples(),qubits=qubit_examples(),xx_large=large,
             xx_ground=xx_ground(),xx_stable_calibration=calibration,trotter=[])
    assert abs(out['ed'][0]['ground_root_fidelity']-out['xx_ground'][0]['ground_root_fidelity'])<1e-12
    for model,L,period in [('XX',4,2),('Heisenberg',4,2),('SU3',3,3)]:
        ns=[1,2,3] if period==2 else [1,2]
        for n in ns:
            for beta in [.5,2.,8.]:
                out['trotter'].append(transfer(model,L,period,n,beta))
    out['seconds']=time.time()-start
    (ROOT/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    summary=dict(seconds=out['seconds'],ed_rows=sum(len(x['rows']) for x in out['ed']),
                 max_identity_error=max(q['identity_error'] for x in out['ed'] for row in x['rows']
                                        for q in row['renyi']),
                 trotter_rows=len(out['trotter']),xx_large_rows=len(out['xx_large']),
                 models={x['model']:dict(beta16=next(r for r in x['rows'] if r['beta']==16.),
                                             ground_conditional_min=x['ground_conditional_min']) for x in out['ed']})
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='models'},indent=2))


if __name__=='__main__': main()
