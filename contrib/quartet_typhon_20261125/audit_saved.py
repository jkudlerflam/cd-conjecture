"""Deterministic audit of historical local harmonic filtering.

No regression standard error is a statistical error bar: adjacent filtered
windows overlap. Signed fits reject zero crossings. Forward fits apply the
same filter to the model power law, explicitly accounting for window bias.
These fits do not identify the emergent charge of a lattice operator.
"""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.linalg import toeplitz, det
from common import BASE,HERE,sha256,write_json

NAMES=['T1_L256_chi384','T1_L256','T1_L384_chi512','T1_L384','T3_L256','T3_L384']


def historical():
    argv=sys.argv[:];sys.argv=sys.argv[:1]
    spec=importlib.util.spec_from_file_location('historical_quartet',BASE/'scripts/ana_quartet.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);sys.argv=argv
    return m


def qs(kappa,nmax):
    return np.array([(m*np.pi+n*kappa)%(2*np.pi) for m in (0,1) for n in range(-nmax,nmax+1)])


def filter_matrix(r,kappa,width=48,nmax=2,degree=0,step=4):
    r=np.asarray(r,float);freq=qs(kappa,nmax);zero=int(np.argmin(abs(freq)))
    centers=[];kernels=[];conds=[]
    # Preserve historical end-exclusion rule for a strict regression.
    for a in range(0,len(r)-width,step):
        rr=r[a:a+width];x=(rr-rr.mean())/width
        w=np.sqrt(np.hanning(width+2)[1:-1])
        X=np.column_stack([x**d*np.exp(1j*q*(rr-rr.mean())) for d in range(degree+1) for q in freq])
        X*=w[:,None]
        kernel=np.linalg.pinv(X,rcond=1e-13)[zero]*w
        row=np.zeros(len(r),complex);row[a:a+width]=kernel
        centers.append(rr.mean());kernels.append(row);conds.append(np.linalg.cond(X))
    return np.array(centers),np.array(kernels),max(conds,default=0.)


def signed_fit(r,g,lo,hi):
    s=(r>=lo)&(r<=hi);x=r[s];y=np.asarray(g)[s].real
    result={'n':int(s.sum()),'sign_changes':0,'log_p':None,'relative_rms':None,'amplitude40':None}
    if not len(y): return result
    result['sign_changes']=int(np.sum(np.sign(y[1:])!=np.sign(y[:-1])))
    result['min_abs']=float(np.min(abs(y)));result['max_abs']=float(np.max(abs(y)))
    if len(y)<4 or np.any(y==0) or np.any(np.sign(y)!=np.sign(y[0])):return result
    p=-np.polyfit(np.log(x),np.log(abs(y)),1)[0]
    f=(x/40)**(-p);A=float(np.dot(f,y)/np.dot(f,f))
    result.update(log_p=float(p),amplitude40=A,relative_rms=float(np.linalg.norm(y-A*f)/np.linalg.norm(y)))
    return result


def forward_fit(r,centers,K,g,lo,hi):
    s=(centers>=lo)&(centers<=hi);y=np.asarray(g).real[s]
    if len(y)<4 or not np.linalg.norm(y):return None
    def err(p,ret=False):
        f=(K@(r/40.)**(-p)).real[s]
        a=np.dot(f,y)/np.dot(f,f)
        return (float(a),float(np.linalg.norm(y-a*f)/np.linalg.norm(y))) if ret else np.linalg.norm(y-a*f)**2/np.linalg.norm(y)**2
    fit=minimize_scalar(err,bounds=(.05,8),method='bounded',options={'xatol':1e-7})
    a,e=err(fit.x,True)
    return {'p':float(fit.x),'amplitude40':a,'relative_rms':e,'at_boundary':bool(fit.x<.06 or fit.x>7.99)}


def raw_matrix(block,ref):
    b=block[str(ref)];vals={}
    for lab in ('md','dm','mm','dd'):
        a=np.array(b[lab]);vals[lab]=a[:,0]+1j*a[:,1]
    md,dm,mm,dd=(vals[k] for k in ('md','dm','mm','dd'))
    C=np.empty((len(md),2,2),complex)
    C[:,0,0]=(mm+md+dm+dd)/4
    C[:,1,1]=-(mm-md-dm+dd)/4
    C[:,0,1]=(mm-md+dm-dd)/4j
    C[:,1,0]=(mm+md-dm-dd)/4j
    return np.arange(len(md))+b['j0']-ref,C


def eigenvalues(C):
    sym=(C.real+np.swapaxes(C.real,-1,-2))/2
    val=np.linalg.eigvalsh(sym)
    small=np.where(abs(val[:,0])<abs(val[:,1]),val[:,0],val[:,1])
    big=np.where(abs(val[:,0])>=abs(val[:,1]),val[:,0],val[:,1])
    return big,small


def audit(out):
    t0=time.time();A=historical();sources={};rows=[];curves={};baseline=[]
    with tempfile.TemporaryDirectory() as tmp:
        for name in NAMES:
            p=BASE/'meas/quartet'/f'{name}.json.gz'
            Path(tmp,f'{name}.json').write_bytes(gzip.decompress(p.read_bytes()))
            sources[str(p.relative_to(BASE.parent.parent))]=sha256(p)
        tab=subprocess.check_output([sys.executable,str(BASE/'scripts/quartet_table.py'),tmp],text=True)
    expected=(BASE/'tables/quartet_k0.txt').read_text();expected=expected[expected.index('== '):]
    assert tab==expected,'Historical table differs'
    (out/'historical_table_reproduced.txt').write_text(tab)
    for name in NAMES:
        m=json.loads(gzip.decompress((BASE/'meas/quartet'/f'{name}.json.gz').read_bytes()))
        kappa=A.fit_kappa(m['P']['2'],m['refs'][0]);A.KAPPA=kappa
        krefs={str(i):float(A.fit_kappa(m['P']['2'],i)/np.pi) for i in m['refs']}
        for kind in ('P','Q'):
            for key,blk in m[kind].items():
                for width in (24,32,40,48,56,64):
                    for nh in (1,2,3):
                        refdata=[]
                        for ref in m['refs']:
                            r,C=raw_matrix(blk,ref)
                            rc,K,condition=filter_matrix(r,kappa,width,nh)
                            if len(K)==0:continue
                            F=np.einsum('ar,rij->aij',K,C)
                            refdata.append((ref,r,rc,K,F,condition))
                        for choice in ('all',*[str(i) for i in m['refs']]):
                            use=refdata if choice=='all' else [x for x in refdata if str(x[0])==choice]
                            if not use:continue
                            n=min(len(x[2]) for x in use);rc=use[0][2][:n]
                            F=np.mean([x[4][:n] for x in use],axis=0)
                            big,small=eigenvalues(F);g=big if kind=='P' else small
                            for lo,hi in ((28,48),(32,56),(40,64),(48,80),(64,96)):
                                fit=signed_fit(rc,g,lo,hi)
                                if fit['n']<4:continue
                                row=dict(dataset=name,L=m['L'],chi=m['chi'],kind=kind,operator=key,
                                         width=width,n_harmonics=nh,refs=choice,kappa_pi=float(kappa/np.pi),
                                         lo=lo,hi=hi,condition=max(x[5] for x in use),**fit)
                                rows.append(row)
                            if width==48 and nh==2 and choice=='all':
                                hist_r,hist_big,hist_small,_=A.quad(blk,m['refs'])
                                assert np.max(abs(g-(hist_big if kind=='P' else hist_small)))<1e-14
                                f=forward_fit(use[0][1],rc,use[0][3][:n],g,32,56)
                                baseline.append({'dataset':name,'kind':kind,'operator':key,
                                                 'signed_fit':signed_fit(rc,g,32,56),'forward_power_fit':f,
                                                 'kappa_by_ref_pi':krefs})
                                curves[f'{name}/{kind}/{key}']={'r':rc.tolist(),'big':big.tolist(),'small':small.tolist(),
                                    'antisymmetric_max':float(np.max(abs(F[:,0,1]-F[:,1,0]))),
                                    'imaginary_max':float(np.max(abs(F.imag)))}
        print('AUDITED',name,flush=True)
    with (out/'sensitivity.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_json(out/'baseline_signed_fits.json',baseline)
    with gzip.open(out/'baseline_curves.json.gz','wt') as f:json.dump(curves,f)
    write_json(out/'audit_provenance.json',{'historical_table_exact_match':True,'input_sha256':sources,
               'rows':len(rows),'seconds':time.time()-t0,'error_bar_policy':'No independent-sample regression errors; overlapping windows are correlated.'})
    return baseline


def controls(out):
    """Known p signals plus exact thermodynamic XX x XX quartet control.

    With two independent zero-field XX chains (physical d=4 per rung),
    C(r)=sin(pi*r/2)/(pi*r). Wick gives one-chain adjacent-pair correlation
    C(r)^2-C(r-1)C(r+1); the two-chain quartet is its square, p=4 (x=2).
    The two free gapless fermion bands give two compact bosons. This is a
    control only, not a new microscopic d=2 candidate.
    """
    r=np.arange(10,201,dtype=float);kappa=.2077*np.pi;results=[];raw={}
    for p in (3.68,4.,4.32):
        charge=1e-8*(r/40)**(-p)
        neutral=3e-7*(r/40)**(-2)
        u=np.array([np.cos(.31),np.sin(.31)]);v=np.array([-u[1],u[0]])
        oscill=3e-7*(r/40)**(-.95)*(np.cos(kappa*r+.4)+.5*np.cos(np.pi*r))
        for label,amp,multi in (('pure_power',0,0),('rank1_plus_oscillation',1,0),('two_neutral_components',1,1)):
            C=charge[:,None,None]*np.eye(2)+neutral[:,None,None]*np.outer(u,u)
            C+=amp*oscill[:,None,None]*np.outer(v,v)
            C+=multi*(2e-8*(r/40)**(-4))[:,None,None]*np.outer(v,v)
            raw[f'{label}_p{p}']=C.tolist()
            for width in (24,32,40,48,56,64):
                rc,K,cond=filter_matrix(r,kappa,width,2)
                F=np.einsum('ar,rij->aij',K,C);_,small=eigenvalues(F)
                results.append({'control':label,'true_charge_p':p,'width':width,'condition':cond,
                                'signed_fit':signed_fit(rc,small,32,56),
                                'forward_power_fit':forward_fit(r,rc,K,small,32,56)})
    # An explicit nonidentifiability construction: a charged isotropic r^-4
    # covariance and two independent neutral r^-4 covariances are identical.
    charge=1e-8*(r/40)**-4
    nonid=float(np.max(abs(charge[:,None,None]*np.eye(2)-charge[:,None,None]*(np.outer(u,u)+np.outer(v,v)))))
    cr=lambda x:np.sinc(x/2)/2
    G=(cr(r)**2-cr(r-1)*cr(r+1))**2
    # P_i=Sm_A_i Sm_B_i. The XX transverse correlator is half the
    # Jordan-Wigner Toeplitz determinant, and independent chains factor.
    pair=[]
    for distance in r.astype(int):
        g=lambda x:2*cr(x)-(np.asarray(x)==0)
        mat=toeplitz(g(1+np.arange(distance)),g(1-np.arange(distance)))
        pair.append((.5*det(mat))**2)
    pair=np.array(pair)
    for width in (24,32,40,48,56,64):
        rc,K,cond=filter_matrix(r,kappa,width,2)
        filtered=(K@G).real
        results.append({'control':'exact_XX_times_XX_quartet_d4','true_charge_p':4.,'width':width,
                        'signed_fit':signed_fit(rc,filtered,32,56),'forward_power_fit':forward_fit(r,rc,K,filtered,32,56)})
        filtered_pair=(K@pair).real
        results.append({'control':'exact_XX_times_XX_pair_d4','true_charge_p':1.,'width':width,
                        'signed_fit':signed_fit(rc,filtered_pair,32,56),
                        'forward_power_fit':forward_fit(r,rc,K,filtered_pair,32,56)})
    write_json(out/'synthetic_controls.json',{'r':r.tolist(),'results':results,'nonidentifiability_max_abs':nonid,
        'two_boson_control':{'model':'two decoupled zero-field XX chains; microscopic d=4','G':G.tolist(),
            'operator':'(Sm_A_i Sm_A_i+1)(Sm_B_i Sm_B_i+1)','asymptotic_p':4,'asymptotic_x':2,
            'pair_operator':'Sm_A_i Sm_B_i','pair_G':pair.tolist(),'pair_asymptotic_p':1},
        'note':'Forward fitting a pure power repairs smoothing bias only; neutral contamination and oscillatory leakage remain hypotheses.'})
    with gzip.open(out/'synthetic_raw.json.gz','wt') as f:json.dump({'r':r.tolist(),'matrices':raw},f)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default=str(HERE/'analysis'));a=p.parse_args()
    out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    audit(out);controls(out)
