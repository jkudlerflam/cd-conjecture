#!/usr/bin/env python3
"""Closed-book D1 checks. All model identifications are qualified in D1_report.md.
No external input files are read. Analysis reads only outputs generated here.
"""
import argparse
import json
import math
import time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, eigs


def qnum(x, gamma):
    return np.sin(x * gamma) / np.sin(gamma)


def spin1_projectors(gamma):
    """U_q(sl2) coproduct projectors; q=exp(i gamma), spin one."""
    z = np.diag([1., 0., -1.])
    eye = np.eye(3)
    e = np.sqrt(qnum(2, gamma)) * np.diag(np.ones(2), 1)
    f = e.T
    k = np.diag(np.exp(1j * gamma * np.diag(z)))
    et = np.kron(e, k) + np.kron(np.linalg.inv(k), e)
    ft = np.kron(f, k) + np.kron(np.linalg.inv(k), f)
    mt = np.diag(np.kron(z, eye) + np.kron(eye, z))
    cas = ft @ et + np.diag(qnum(mt, gamma) * qnum(mt + 1, gamma))
    vals = [qnum(j, gamma) * qnum(j + 1, gamma) for j in range(3)]
    ps = []
    for j in range(3):
        p = np.eye(9, dtype=complex)
        for k0 in range(3):
            if j != k0:
                p = p @ (cas - vals[k0] * np.eye(9)) / (vals[j] - vals[k0])
        ps.append(p)
    return ps


def baxter_r(x, gamma, a=4, b=6, product=False):
    """Candidate twisted spin-one baxterization, checked by YBE below."""
    p0, p1, p2 = spin1_projectors(gamma)
    q = np.exp(1j * gamma)
    r1 = (q**a * x - 1) / (q**a - x)
    r0 = (q**b * x + 1) / (q**b + x)
    if product:
        r0 *= r1
    return p2 + r1 * p1 + r0 * p0


def ybe_residual(gamma, a=4, b=6, product=False):
    x, y = 1.17, .83
    r = lambda t: baxter_r(t, gamma, a, b, product)
    eye = np.eye(3)
    lhs = np.kron(r(x), eye) @ np.kron(eye, r(x*y)) @ np.kron(r(y), eye)
    rhs = np.kron(eye, r(y)) @ np.kron(r(x*y), eye) @ np.kron(eye, r(x))
    return float(np.linalg.norm(lhs-rhs) / np.linalg.norm(lhs))


def local_sum(local, sites, d, length):
    """Embed local matrices at ordered sites into a periodic spin Hilbert space."""
    n = d ** length
    powers = d ** np.arange(length-1, -1, -1, dtype=np.int64)
    states = np.arange(n, dtype=np.int64)
    digits = (states[:, None] // powers[None, :]) % d
    rows, cols, data = [], [], []
    k = len(sites)
    lp = d ** np.arange(k-1, -1, -1, dtype=np.int64)
    ins = digits[:, sites] @ lp
    outdigits = (np.arange(d**k)[:, None] // lp[None, :]) % d
    for out, inp in zip(*np.where(np.abs(local) > 1e-13)):
        col = states[ins == inp]
        if len(col):
            change = (outdigits[out] - outdigits[inp]) @ powers[list(sites)]
            rows.extend(col + change)
            cols.extend(col)
            data.extend(np.full(len(col), local[out, inp], dtype=complex))
    return sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def periodic_nn(h, d, length):
    out = sparse.csr_matrix((d**length, d**length), dtype=complex)
    for j in range(length):
        out += local_sum(h, (j, (j+1)%length), d, length)
    return out


def antihermitian_ratio(h):
    if sparse.issparse(h):
        return float(sparse.linalg.norm(h-h.getH()) / max(sparse.linalg.norm(h), 1e-30))
    return float(np.linalg.norm(h-h.conj().T) / max(np.linalg.norm(h), 1e-30))


def ik_density(gamma):
    """i d/dlog(x) Rcheck(x)|x=1 for exponents (4,6)."""
    p0, p1, _ = spin1_projectors(gamma)
    return 1 / np.tan(2*gamma) * p1 - np.tan(3*gamma) * p0


def tl_e(gamma):
    e = np.zeros((4,4), dtype=complex)
    e[1:3,1:3] = [[np.exp(-1j*gamma), -1], [-1, np.exp(1j*gamma)]]
    return e


def tl_three_density(gamma):
    """H=sum [2 cos(gamma)e_j - {e_j,e_{j+1}}]."""
    e = tl_e(gamma)
    a, b = np.kron(e, np.eye(2)), np.kron(np.eye(2), e)
    return 2*np.cos(gamma)*a - a@b - b@a


def periodic_tl(gamma, length):
    h = tl_three_density(gamma)
    out = sparse.csr_matrix((2**length, 2**length), dtype=complex)
    for j in range(length):
        out += local_sum(h, (j, (j+1)%length, (j+2)%length), 2, length)
    return out


def graded_swap(m,n):
    d = m+n
    p = np.zeros((d*d,d*d))
    for a in range(d):
        for b in range(d):
            p[b*d+a,a*d+b] = -1 if a>=m and b>=m else 1
    return p


def osp_density(m,n):
    d = m+2*n
    delta = m-2*n
    if delta == 2:
        raise ValueError('Singular standard rational normalization at m-2n=2')
    vb = np.zeros(d*d)
    vf = np.zeros(d*d)
    for a in range(m):
        vb[a*d+a] = 1
    for r in range(n):
        a,b = m+2*r,m+2*r+1
        vf[a*d+b],vf[b*d+a] = 1,-1
    e = np.outer(vb+vf,vb-vf)
    return graded_swap(m,2*n) - 2/(delta-2)*e, e


def algebra_checks():
    g=.19
    ps=spin1_projectors(g)
    report = {'gamma':g,
        'projectors': {'sum':float(np.linalg.norm(sum(ps)-np.eye(9))),
                       'idempotence_max':float(max(np.linalg.norm(p@p-p) for p in ps))},
        'baxterization_trials': []}
    for a in [2,4,-2,-4]:
        for b in [3,6,-3,-6]:
            for prod in [False,True]:
                report['baxterization_trials'].append({'a':a,'b':b,'product':prod,
                    'relative_ybe_residual':ybe_residual(g,a,b,prod)})
    report['baxterization_trials'].sort(key=lambda z:z['relative_ybe_residual'])
    h=ik_density(g)
    report['ik_local_antihermitian_ratio']=antihermitian_ratio(h)
    report['ik_periodic_antihermitian_ratio']={str(L):antihermitian_ratio(periodic_nn(h,3,L)) for L in [3,4]}
    e=tl_e(g)
    a,b=np.kron(e,np.eye(2)),np.kron(np.eye(2),e)
    report['tl_relations']={'e_squared':float(np.linalg.norm(e@e-2*np.cos(g)*e)),
        'e_e_e':float(np.linalg.norm(a@b@a-a))}
    report['tl_periodic_antihermitian_ratio']={str(L):antihermitian_ratio(periodic_tl(g,L)) for L in [4,6]}
    report['osp']={}
    for m,n in [(3,0),(0,2),(1,1),(2,1),(3,1)]:
        h,e=osp_density(m,n)
        report['osp'][f'{m}|{2*n}']={'local_antihermitian_ratio':antihermitian_ratio(h),
          'contraction_relation':float(np.linalg.norm(e@e-(m-2*n)*e)),
          'periodic_L3_antihermitian_ratio':antihermitian_ratio(periodic_nn(h,m+2*n,3))}
    return report



# Callable spectrum routines, also used by d1_run_ed.py.
def sector_basis(d, length, charge):
    powers = d**np.arange(length-1,-1,-1,dtype=np.int64)
    allstates = np.arange(d**length,dtype=np.int64)
    digits = (allstates[:,None]//powers[None,:])%d
    values = (d-1)-2*digits
    keep = values.sum(axis=1)==charge
    return allstates[keep], digits[keep], powers


def sector_hamiltonian(local, d, length, charge=0):
    k = int(round(np.log(len(local))/np.log(d)))
    states,digits,powers = sector_basis(d,length,charge)
    lp=d**np.arange(k-1,-1,-1,dtype=np.int64)
    od=(np.arange(d**k)[:,None]//lp[None,:])%d
    rows,cols,data=[],[],[]
    transitions=list(zip(*np.where(abs(local)>1e-12)))
    for j in range(length):
        sites=tuple((j+a)%length for a in range(k))
        ins=digits[:,sites]@lp
        for out,inp in transitions:
            col=np.where(ins==inp)[0]
            if not len(col): continue
            new=states[col]+(od[out]-od[inp])@powers[list(sites)]
            row=np.searchsorted(states,new)
            if np.any(row>=len(states)) or np.any(states[row]!=new):
                raise ValueError('Local operator does not conserve specified charge')
            rows.append(row);cols.append(col);data.append(np.full(len(col),local[out,inp]))
    if not data:
        return sparse.csr_matrix((len(states),len(states)),dtype=complex),states,powers
    h=sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(len(states),len(states))).tocsr()
    return h,states,powers


def spectrum(model, length, param, charge=0, count=12):
    d=3 if model=='ik' else 2
    local=ik_density(param) if model=='ik' else (1 if model=='tl' else -1)*tl_three_density(param)
    start=time.time()
    h,states,powers=sector_hamiltonian(local,d,length,charge)
    if h.shape[0]<300:
        values,vectors=np.linalg.eig(h.toarray())
    else:
        rng=np.random.default_rng(1907)
        v0=rng.normal(size=h.shape[0])+1j*rng.normal(size=h.shape[0])
        values,vectors=eigs(h,k=min(count,h.shape[0]-2),which='SR',tol=1e-10,maxiter=6000,v0=v0)
    rotated=(states % powers[0])*d+states//powers[0]
    shift=np.searchsorted(states,rotated)
    # Resolve translations within every returned energy-degenerate subspace.
    # A truncated degenerate multiplet is flagged by its translation residual.
    unused=set(range(len(values)))
    while unused:
        i=min(unused)
        group=[j for j in unused if abs(values[j]-values[i])<1e-7]
        unused.difference_update(group)
        q,_=np.linalg.qr(vectors[:,group])
        tq=np.empty_like(q);tq[shift,:]=q
        phases,w=np.linalg.eig(q.conj().T@tq)
        vectors[:,group]=q@w
    order=np.lexsort((values.imag,values.real))[:count]
    values,vectors=values[order],vectors[:,order]
    momenta=[]
    variances=[]
    residuals=[]
    for v,val in zip(vectors.T,values):
        tv=np.empty_like(v);tv[shift]=v
        ov=np.vdot(v,tv)/np.vdot(v,v)
        momenta.append(float(np.angle(ov)))
        variances.append(float(np.linalg.norm(tv-ov*v)/np.linalg.norm(v)))
        residuals.append(float(np.linalg.norm(h@v-val*v)/np.linalg.norm(v)))
    return {'model':model,'L':length,'parameter':param,'charge_2Sz':charge,'dimension':h.shape[0],
      'eigenvalues':[[float(v.real),float(v.imag)] for v in values],
      'translation_phases_unresolved_if_degenerate':momenta,
      'translation_residuals':variances,'eigen_residual_max':max(residuals),
      'antihermitian_ratio':antihermitian_ratio(h),'seconds':time.time()-start}


def momentum_block(h, states, powers, d, length, momentum_integer):
    rotated=(states % powers[0])*d+states//powers[0]
    shift=np.searchsorted(states,rotated)
    seen=np.zeros(len(states),dtype=bool)
    rows=[];cols=[];data=[];column=0
    for i in range(len(states)):
        if seen[i]: continue
        orbit=[];j=i
        while not seen[j]:
            orbit.append(j);seen[j]=True;j=shift[j]
        ell=len(orbit)
        if (momentum_integer*ell)%length: continue
        rows.extend(orbit);cols.extend([column]*ell)
        data.extend(np.exp(-2j*np.pi*momentum_integer*np.arange(ell)/length)/np.sqrt(ell))
        column+=1
    u=sparse.coo_matrix((data,(rows,cols)),shape=(len(states),column)).tocsr()
    return (u.getH()@h@u).tocsr(),u


def momentum_spectrum(model,length,param,charge=0,momenta=(0,1),count=8):
    d=3 if model=='ik' else 2
    local=ik_density(param) if model=='ik' else tl_three_density(param)
    h,states,powers=sector_hamiltonian(local,d,length,charge)
    out=[]
    for r in momenta:
        hp,u=momentum_block(h,states,powers,d,length,r)
        if hp.shape[0]<100:
            vals,vecs=np.linalg.eig(hp.toarray())
        else:
            rng=np.random.default_rng(1907)
            v0=rng.normal(size=hp.shape[0])+1j*rng.normal(size=hp.shape[0])
            vals,vecs=eigs(hp,k=min(count,hp.shape[0]-2),which='SR',tol=1e-11,maxiter=10000,v0=v0)
        order=np.argsort(vals.real)[:count];vals=vals[order];vecs=vecs[:,order]
        residuals=[float(np.linalg.norm(h@(u@v)-val*(u@v))/np.linalg.norm(u@v)) for val,v in zip(vals,vecs.T)]
        out.append({'model':model,'L':length,'parameter':param,'charge_2Sz':charge,
          'momentum_integer':r,'dimension':hp.shape[0],
          'eigenvalues':[[float(v.real),float(v.imag)] for v in vals],
          'full_sector_residual_max':max(residuals)})
    return out



def run_ed():
    """Small-sector ED, largest sector dimension 184756."""
    import json
    from pathlib import Path
    out=[]
    for model,param,lengths in [('ik',.19,[4,6,8,10,12]),('tl',.4,[4,6,8,10,12,14,16,18,20]),('tl_af',0.,[4,8,12,16])]:
     for L in lengths:
      for charge in ([0,2] if model in ['ik','tl'] else [0]):
       result=spectrum(model,L,param,charge,count=18)
       out.append(result)
       Path('d1_spectra.json').write_text(json.dumps(out,indent=2)+'\n')
       print(model,L,charge,'dim',result['dimension'],'E0',result['eigenvalues'][0],
             'gap',result['eigenvalues'][1][0]-result['eigenvalues'][0][0],
             'residual',result['eigen_residual_max'],'seconds',round(result['seconds'],2),flush=True)


def run_momentum():
    import json
    from pathlib import Path
    out=[]
    for model,param,lengths in [('ik',.19,[4,6,8,10,12]),('tl',.4,[4,8,12,16,20])]:
     for L in lengths:
      result=momentum_spectrum(model,L,param,count=10)
      out.extend(result)
      Path('d1_momenta.json').write_text(json.dumps(out,indent=2)+'\n')
      print(model,L,[(r['momentum_integer'],r['eigenvalues'][0],r['full_sector_residual_max']) for r in result],flush=True)


def run_audit():
    """Check the remaining conserved-charge sectors at the largest ED sizes."""
    import json,time
    import numpy as np
    from scipy.sparse.linalg import eigs
    from pathlib import Path
    out=[]
    for model,d,L,g,local in [('ik',3,12,.19,ik_density(.19)),('tl',2,20,.4,tl_three_density(.4))]:
     for charge in range(0,(d-1)*L+1,2):
      start=time.time();h,_,_=sector_hamiltonian(local,d,L,charge)
      if h.shape[0]<=20: vals,vecs=np.linalg.eig(h.toarray())
      else:
       rng=np.random.default_rng(1907)
       v0=rng.normal(size=h.shape[0])+1j*rng.normal(size=h.shape[0])
       vals,vecs=eigs(h,k=2,which='SR',tol=1e-10,maxiter=10000,v0=v0)
      a=int(np.argmin(vals.real));v=vecs[:,a];z=vals[a]
      r={'model':model,'L':L,'charge_2Sz':charge,'dimension':h.shape[0],
         'lowest_real_part_eigenvalue':[float(z.real),float(z.imag)],
         'residual':float(np.linalg.norm(h@v-z*v)/np.linalg.norm(v)), 'seconds':time.time()-start}
      out.append(r);Path('d1_ground_audit.json').write_text(json.dumps(out,indent=2)+'\n')
      print(r,flush=True)


def run_analyze():
    """Finite-size diagnostics. Velocity-based c estimates are explicitly conditional."""
    import json, math
    from pathlib import Path
    import numpy as np
    r=json.loads(Path('d1_momenta.json').read_text())
    results={}
    for model in ['ik','tl']:
     rows=[]
     for a in r:
      if a['model']!=model or a['momentum_integer']!=0:continue
      L=a['L'];b=next(t for t in r if t['model']==model and t['L']==L and t['momentum_integer']==1)
      ev=np.array(a['eigenvalues'])[:,0]; e0=ev[0]
      gaps=L*(ev-e0)
      nonzero=gaps[gaps>1e-7]
      vp=L*(b['eigenvalues'][0][0]-e0)/(2*math.pi)
      rows.append({'L':L,'E0':e0,'L_neutral_gap':float(nonzero[0]),
        'neutral_gap_ratios':list(nonzero[:4]/nonzero[0]),'velocity_descendant_assumption':vp})
     rows.sort(key=lambda a:a['L'])
     v_recalled=math.pi/(math.pi-6*.19) if model=='ik' else 2*math.pi*math.sin(.8)/(math.pi-.8)
     for i,a in enumerate(rows):
      if i==0:continue
      b=rows[i-1]
      A=-(a['E0']/a['L']-b['E0']/b['L'])/(1/a['L']**2-1/b['L']**2)
      a['pair_cv_eff']=6*A/math.pi
      a['pair_c_eff_using_descendant_velocity']=6*A/(math.pi*a['velocity_descendant_assumption'])
      a['pair_c_eff_using_unverified_recalled_velocity']=6*A/(math.pi*v_recalled)
     # Test linearity of 1/sqrt(Lgap) against log L; it cannot establish asymptotics.
     selected=rows[1:]
     x=np.log([a['L'] for a in selected]);y=1/np.sqrt([a['L_neutral_gap'] for a in selected])
     coef=np.polyfit(x,y,1);pred=np.polyval(coef,x)
     logfit={'slope':coef[0],'intercept':coef[1],'max_y_residual':max(abs(y-pred)),
       'warning':'Short-size consistency test only; power-law alternatives are not excluded.'}
     results[model]={'velocity_formula_status':'recalled/guessed, not derived or externally verified',
                      'velocity_candidate':v_recalled,'rows':rows,'log_gap_test':logfit}
    Path('d1_diagnostics.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


def run_transfer():
    """Direct small-size alternating six-vertex transfer-matrix check."""
    import numpy as np

    def transfer(u,g,L):
        e=tl_e(g);p=graded_swap(2,0)
        mon=np.eye(2**(L+1),dtype=complex)
        for j in range(L):
            w=u-(np.pi/2 if j%2 else 0)
            r=p@(np.sin(g-w)*np.eye(4)+np.sin(w)*e)
            mon=local_sum(r,(0,j+1),2,L+1)@mon
        return np.trace(mon.reshape(2,2**L,2,2**L),axis1=0,axis2=2)

    g=.4;L=4;eps=1e-6
    H=np.zeros((2**L,2**L),dtype=complex)
    for u in [0,np.pi/2]:
        t=transfer(u,g,L)
        dt=(transfer(u+eps,g,L)-transfer(u-eps,g,L))/(2*eps)
        H+=np.linalg.solve(t,dt)
    A=periodic_tl(g,L).toarray();eye=np.eye(2**L)
    c=np.linalg.lstsq(np.column_stack([A.ravel(),eye.ravel()]),H.ravel(),rcond=None)[0]
    print('fit Htransfer = alpha Htl + beta I:',c)
    print('relative residual',np.linalg.norm(H-c[0]*A-c[1]*eye)/np.linalg.norm(H))
    print('relative antihermitian part',np.linalg.norm(H-H.conj().T)/np.linalg.norm(H))
    print('expected alpha',1/(np.sin(g)*np.cos(g)))
    print('expected beta',-2*L/np.tan(2*g))


def run_symbolic():
    """Symbolic checks for D1; no data inputs."""
    import sympy as s
    q=s.symbols('q',nonzero=True)
    qn=lambda j: (q**j-q**(-j))/(q-q**(-1))
    k=s.diag(q,1,1/q); ki=k.inv()
    e=s.Matrix([[0,1,0],[0,0,1],[0,0,0]]); f=e.T
    et=s.kronecker_product(e,k)+s.kronecker_product(ki,e)
    ft=s.kronecker_product(f,k)+s.kronecker_product(ki,f)
    ms=[2,1,0,1,0,-1,0,-1,-2]
    cas=(q+1/q)*ft*et+s.diag(*[qn(m)*qn(m+1) for m in ms])
    vals=[qn(j)*qn(j+1) for j in range(3)]
    p=[]
    for j in range(3):
        a=s.eye(9)
        for l in range(3):
            if l!=j: a=a*(cas-vals[l]*s.eye(9))/(vals[j]-vals[l])
        p.append(a.applyfunc(s.factor))
    minimal=(cas-vals[0]*s.eye(9))*(cas-vals[1]*s.eye(9))*(cas-vals[2]*s.eye(9))
    print("Casimir minimal polynomial identically zero:", all(s.factor(v)==0 for v in minimal))
    print("Projector traces:", [s.factor(s.trace(a)) for a in p])
    h=s.I*(q**4+1)/(q**4-1)*p[1]+s.I*(q**6-1)/(q**6+1)*p[0]
    for i,j in [(2,4),(4,2),(4,6),(1,1),(2,2)]:
        v=s.factor(h[i,j]); anti=s.factor(v+v.subs(q,1/q))
        print('IK',i,j,'h=',v,'h-hdag=',anti)
    c,t=s.symbols('c t',real=True)
    E=s.zeros(4); E[1,1]=c-s.I*t;E[2,2]=c+s.I*t;E[1,2]=E[2,1]=-1
    A=s.kronecker_product(E,s.eye(2));B=s.kronecker_product(s.eye(2),E)
    H=2*c*A-A*B-B*A
    paulis=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    for a in range(4):
     for b in range(4):
      for d in range(4):
       coeff=s.expand(s.trace(H*s.kronecker_product(paulis[a],paulis[b],paulis[d]))/8)
       if coeff!=0: print('TL','IXYZ'[a]+'IXYZ'[b]+'IXYZ'[d],s.factor(coeff))



def run_validate():
    """Independent algebra, blocking, symmetry, and complex-spectrum checks."""
    report={}
    theta=.19
    h=ik_density(theta)
    z=np.fliplr(np.eye(3))
    z2=np.kron(z,z)
    report['ik_spin_reversal_conjugation_residual']=float(np.linalg.norm(z2@h@z2-h.conj()))
    entry=h[2,4]-h.conj().T[2,4]
    report['ik_analytic_antihermitian_entry_residual']=float(abs(entry-1j*np.cos(2*theta)/np.cos(3*theta)))
    step=1e-6
    derivative=1j*(baxter_r(1+step,theta)-baxter_r(1-step,theta))/(2*step)
    report['ik_derivative_residual']=float(np.linalg.norm(h-derivative))
    gamma=.4
    g=tl_three_density(gamma)
    flip=np.array([[0,1],[1,0]])
    flip3=np.kron(np.kron(flip,flip),flip)
    report['tl_spin_reversal_conjugation_residual']=float(np.linalg.norm(flip3@g@flip3-g.conj()))
    blocked=np.kron(g,np.eye(2))+np.kron(np.eye(2),g)
    hb=periodic_nn(blocked,4,3)
    ht=periodic_tl(gamma,6)
    report['blocking_residual']=float(sparse.linalg.norm(hb-ht))
    ident=np.eye(2);x=np.array([[0,1],[1,0]])
    y=np.array([[0,-1j],[1j,0]]);zz=np.diag([1,-1])
    kron3=lambda a,b,c:np.kron(np.kron(a,b),c)
    pauli=(np.cos(2*gamma)/2*np.eye(8)+np.sin(gamma)**2*kron3(zz,zz,ident)
           -.5*sum(kron3(a,ident,a) for a in [x,y,zz])
           +.5j*np.sin(gamma)*sum(kron3(a,a,zz)-kron3(zz,a,a) for a in [x,y]))
    hp=sparse.csr_matrix(ht.shape,dtype=complex)
    for j in range(6):hp+=local_sum(pauli,(j,(j+1)%6,(j+2)%6),2,6)
    report['tl_periodic_pauli_expansion_residual']=float(sparse.linalg.norm(ht-hp))
    zero=periodic_tl(0.,6)
    expected=sparse.csr_matrix(zero.shape,dtype=complex)
    for j in range(6):expected+=local_sum(np.eye(4)-graded_swap(2,0),(j,(j+2)%6),2,6)
    report['tl_zero_gamma_permutation_limit_residual']=float(sparse.linalg.norm(zero-expected))
    impurity=[]
    for m,n,base,other,sign,hop in [(1,1,1,0,1,1),(2,1,2,0,1,1),
                                    (0,2,0,1,1,-1),(1,1,0,1,-1,-1)]:
        d=m+n;L=6
        H=sign*periodic_nn(graded_swap(m,n),d,L)
        powers=d**np.arange(L-1,-1,-1)
        states=base*sum(powers)+(other-base)*powers
        block=H[states,:][:,states].toarray()+L*np.eye(L)
        exact=2*np.eye(L)
        for j in range(L):exact[j,(j+1)%L]+=hop;exact[j,(j-1)%L]+=hop
        impurity.append({'m':m,'n':n,'sign':sign,'residual':float(np.linalg.norm(block-exact))})
    report['graded_impurity_checks']=impurity
    complex_checks=[]
    for model,L,H in [('ik',3,periodic_nn(h,3,3)),('ik',4,periodic_nn(h,3,4)),('tl',6,ht)]:
        vals,vecs=np.linalg.eig(H.toarray());i=int(np.argmax(abs(vals.imag)))
        val=vals[i];v=vecs[:,i]
        complex_checks.append({'model':model,'L':L,'eigenvalue':[float(val.real),float(val.imag)],
          'residual':float(np.linalg.norm(H@v-val*v)/np.linalg.norm(v)),
          'antihermitian_ratio':antihermitian_ratio(H)})
    report['complex_spectrum']=complex_checks
    Path('d1_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    Path('d1_complex_spectrum.json').write_text(json.dumps(complex_checks,indent=2)+'\n')
    print(json.dumps(report,indent=2))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['algebra','ed','momentum','audit','analyze','transfer','symbolic','validate'])
    args=parser.parse_args()
    if args.command=='algebra':
        report=algebra_checks()
        Path('d1_algebra.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2))
    else:
        actions={'ed':run_ed,'momentum':run_momentum,'audit':run_audit,
                 'analyze':run_analyze,'transfer':run_transfer,'symbolic':run_symbolic,'validate':run_validate}
        actions[args.command]()


if __name__=='__main__':
    main()
