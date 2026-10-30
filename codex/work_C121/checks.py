"""C121. Small full ED and linear-cost XX checks. All output stays here."""
from pathlib import Path
from itertools import product
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parent

def bond(model):
    if model == 'SU3':
        return np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9),3,1.
    x=np.array([[0.,1.],[1.,0.]])/2
    y=np.array([[0.,-1j],[1j,0.]])/2
    z=np.diag([.5,-.5])
    h=np.kron(x,x)+np.kron(y,y)
    if model=='Heisenberg': h+=np.kron(z,z)
    return h.real,2,.5

def ring(model,n):
    h,d,m=bond(model)
    powers=d**np.arange(n,dtype=np.int64)
    dim=d**n
    a=np.arange(dim)
    out=np.zeros((dim,dim))
    for i in range(n):
        j=(i+1)%n
        u=(a//powers[i])%d
        v=(a//powers[j])%d
        for aa in range(d):
            for bb in range(d):
                target=a+(aa-u)*powers[i]+(bb-v)*powers[j]
                out[target,a]+=h[d*aa+bb,d*u+v]
    ev,vec=eigh(out)
    return dict(ev=ev,vec=vec,d=d,M=m,n=n,H=out)

def propagator(spec,t):
    return (spec['vec']*np.exp(-t*(spec['ev']-spec['ev'][0])))@spec['vec'].conj().T

def measured(spec,beta,count,sites,histories=False):
    d,n=spec['d'],spec['n']
    labels=np.zeros(d**n,dtype=int)
    for j,site in enumerate(sites):
        labels+=d**j*((np.arange(d**n)//d**site)%d)
    mask=labels[:,None]==labels[None,:]
    step=beta/count
    f=propagator(spec,step/4)
    e=f@f
    rho=e.copy()
    for j in range(count):
        rho*=mask
        if j<count-1: rho=e@rho@e
    value=float(np.real(np.sum(rho*e.T)))
    result=dict(logQ=float(np.log(value)))
    if histories:
        vals=[]
        for history in product(range(d**len(sites)),repeat=count):
            w=f.copy()
            for j,s in enumerate(history):
                w=(labels==s)[:,None]*w
                w=(e if j<count-1 else f)@w
            vals.append(float(np.sum(np.abs(w)**2)))
        assert abs(sum(vals)-value)<2e-10
        result.update(logS=float(np.log(np.dot(vals,vals))),
                      history_residual=float(abs(sum(vals)-value)),
                      min_history=min(vals))
    return result

def coherent(spec,beta,count):
    d,n=spec['d'],spec['n']
    labels=np.arange(d**n)%d
    e=propagator(spec,beta/count)
    vals=[]
    for history in product(range(d),repeat=count):
        w=np.eye(d**n)
        for label in history:
            w=(labels==label)[:,None]*(e@w)
        vals.append(float(np.trace(w)))
    target=float(np.sum(np.exp(-beta*(spec['ev']-spec['ev'][0]))))
    assert abs(sum(vals)-target)<2e-10
    return float(np.log(np.dot(vals,vals))),float(abs(sum(vals)-target)),min(vals)

def xx_parts(n,beta):
    """Ground-shifted full spin trace and one-pinch echo. Even rings only."""
    pieces=[]
    for shift,wanted in [(0.5,0),(0.,1)]:
        eps=np.cos(2*np.pi*(np.arange(n)+shift)/n)
        eps[np.abs(eps)<1e-13]=0.
        occupied=eps<0
        parity=int(np.count_nonzero(occupied)%2)
        even,odd=0.,-np.inf
        grand=0.
        for energy in np.abs(eps):
            x=-beta*energy
            even,odd=np.logaddexp(even,odd+x),np.logaddexp(odd,even+x)
            grand+=np.logaddexp(0.,x)
        # Which grid has the lower filled sea depends on n modulo four.
        lower=(shift==.5) if n%4==0 else (shift==0.)
        offset=0. if lower else np.tan(np.pi/(2*n))
        physical=(even,odd)[wanted^parity]-beta*offset
        # grand-canonical midpoint density echo, with a stable sech.
        x=beta*np.abs(eps)/2
        g=np.sum(np.exp(-x)/(1+np.exp(-2*x)))
        logecho=-np.inf if g==0 else grand-beta*offset+np.log(4*g*g/(n*n))
        pieces.append((physical,logecho,grand-beta*offset))
    logw=float(np.logaddexp(pieces[0][0],pieces[1][0]))
    logecho=float(np.logaddexp(pieces[0][1],pieces[1][1])-np.log(2)-logw)
    echo=float(np.exp(logecho))
    logfraction=float(np.log1p(echo)-np.log(2))
    return logw,logfraction,echo

def xx_ratio(n,beta):
    w,c,echo=xx_parts(n,beta)
    w2,_,_=xx_parts(2*n,beta)
    x=np.pi/(2*n)
    delta=2*np.sin(x/2)**2/(np.sin(x)*np.cos(x))
    return float(beta*delta+2*w-w2),c,echo

def zeno_xx(n,beta):
    # The seam site is isolated and the remaining n-1 sites form an open chain.
    eps=np.cos(np.pi*np.arange(1,n)/n)
    eps[np.abs(eps)<1e-13]=0.
    wz=np.log(2)+np.sum(np.logaddexp(0.,-beta*np.abs(eps)))
    w,_,_=xx_parts(n,beta)
    # E_zeno - E_ring = 1/2 + (1/2) tan(pi/(2n)).
    defect=.5+.5*np.tan(np.pi/(2*n))
    return float(beta*defect+w-wz),float(defect)

def main():
    start=time.time()
    out=dict(ed=[],xx=[],zeno=[],calibration=[],algebra=[])
    max_echo_error=0.
    for model,n in [('XX',4),('Heisenberg',4),('SU3',3)]:
        small,big=ring(model,n),ring(model,2*n)
        d,m=small['d'],small['M']
        delta=float(big['ev'][0]-2*small['ev'][0])
        for a in [.25,1.,2.,4.]:
            for count in [1,2,4,8]:
                beta=count*a*n/m
                w=float(logsumexp(-beta*(small['ev']-small['ev'][0])))
                w2=float(logsumexp(-beta*(big['ev']-big['ev'][0])))
                q=measured(small,beta,count,[0],histories=count<=4)
                q2=measured(big,beta,count,[0,n])
                loss=w-q['logQ']
                loss2=w2-q2['logQ']
                assert -2e-9<=loss<=count*np.log(d)+2e-9
                assert -2e-9<=loss2<=2*count*np.log(d)+2e-9
                r=beta*delta+2*w-w2
                rm=beta*delta+2*q['logQ']-q2['logQ']
                row=dict(model=model,L=n,N=count,a=a,beta=beta,R=r,
                         loss_L=loss,loss_2L=loss2,measured_R=rm,
                         G_lower=rm-count*np.log(d))
                if 'logS' in q:
                    g=beta*delta+q['logS']-q2['logQ']
                    cs,cr,cm=coherent(small,beta,count)
                    row.update(coherent_G=beta*delta+cs-w2,coherent_residual=cr,
                               min_coherent_history=cm,G=g,history_residual=q['history_residual'],
                               G_unmeasured=beta*delta+q['logS']-w2)
                    assert rm<=count*np.log(d)+g+2e-9
                    assert r<=3*count*np.log(d)+g+2e-9
                out['ed'].append(row)
                if model=='XX' and count==1:
                    wr,lc,_=xx_parts(n,beta)
                    err=abs(lc+loss)
                    max_echo_error=max(max_echo_error,err,abs(wr-w))
        # Zeno convergence on the small ring, modest dimensions only.
        for count in [1,4,16,64,256]:
            beta=4.
            q=measured(small,beta,count,[0])
            labels=np.arange(d**n)%d
            hz=small['H']*(labels[:,None]==labels[None,:])
            ez=np.linalg.eigvalsh(hz)
            logz=float(logsumexp(-beta*small['ev']))
            logzz=float(logsumexp(-beta*ez))
            if model=='XX':
                analytic_loss,_=zeno_xx(n,beta)
                assert abs(logz-logzz-analytic_loss)<1e-11
            out['zeno'].append(dict(model=model,L=n,N=count,beta=beta,
                                   loss=logz+beta*small['ev'][0]-q['logQ'],
                                   zeno_loss=logz-logzz))
    # Independent even-length and finite-temperature calibration of the echo formula.
    for n in [4,6,8]:
        spec=ring('XX',n)
        for beta in [.01,.2,1.,4.,16.]:
            q=measured(spec,beta,1,[0])
            w=float(logsumexp(-beta*(spec['ev']-spec['ev'][0])))
            wf,cf,_=xx_parts(n,beta)
            err=max(abs(w-wf),abs(q['logQ']-w-cf))
            max_echo_error=max(max_echo_error,err)
            out['calibration'].append(dict(L=n,beta=beta,error=err))
    assert max_echo_error<2e-10
    out['max_echo_error']=max_echo_error
    for n in [4,8,32,128,512,2048,8192]:
        for a in [.25,1.,2.,4.,16.]:
            beta=2*a*n
            r,c,echo=xx_ratio(n,beta)
            lossz,defect=zeno_xx(n,beta)
            assert -np.log(2)-1e-12<=c<=1e-12
            assert r<=2*np.log(2)+5*np.pi*beta/(3*n)+2e-8
            out['xx'].append(dict(L=n,a=a,beta=beta,R=r,pinch_loss=-c,
                                  echo=echo,G_unmeasured=r+2*c-np.log(2),
                                  zeno_loss=lossz,zeno_defect=defect))
    # General complex matrices test the channel theorem without locality or signs.
    rng=np.random.default_rng(121)
    raw=rng.normal(size=(9,9))+1j*rng.normal(size=(9,9))
    for name,h in [('complex',(raw+raw.conj().T)/10),
                   ('commuting',np.diag(np.arange(9)/8))]:
        ev,vec=eigh(h)
        spec=dict(ev=ev,vec=vec,d=3,n=2)
        for beta in [.2,2.,20.]:
            for count in [1,2,4]:
                q=measured(spec,beta,count,[0],histories=True)
                w=float(logsumexp(-beta*(ev-ev[0])))
                loss=w-q['logQ']
                assert -1e-10<=loss<=count*np.log(3)+1e-10
                if name=='commuting': assert abs(loss)<1e-10
                out['algebra'].append(dict(case=name,beta=beta,N=count,loss=loss,
                                           history_residual=q['history_residual']))
    out['seconds']=time.time()-start
    (ROOT/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    summary=dict(seconds=out['seconds'],max_echo_error=max_echo_error,
                 ed_cases=len(out['ed']),large_xx_cases=len(out['xx']),
                 positivity_failures=int(sum(x.get('G',x['G_lower'])>1e-8 for x in out['ed'])),
                 largest_xx=[x for x in out['xx'] if x['L']==8192],
                 examples=[x for x in out['ed'] if x['N']==1 and x['a']==1.])
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
