"""Modest literal eleven-state checks, with one BLAS thread.

No long certificate scan and no inference of a microscopic CFT.
The original C84 files are read only. Outputs are local JSON.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
from itertools import product, combinations_with_replacement
import hashlib
import json
import sys
import time
import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, splu
import model as m

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT/'codex/work_C84'))
import bond11 as old
import bond18 as old_parent
from parent_checks import gaussian_residual
mp.mp.dps = 80


def lowest(H, count=1):
    if H.shape[0] < 150:
        e, v = la.eigh(H.toarray(), subset_by_index=[0, count-1])
    else:
        e, v = eigsh(H, k=count, which='SA', tol=1e-12,
                     v0=np.random.default_rng(1117).normal(size=H.shape[0]))
        order = np.argsort(e); e, v = e[order], v[:, order]
    residual = float(np.max(la.norm(H@v-v*e, axis=0)))
    assert residual < 2e-8
    return e, v, residual


def literal_checks():
    p = old_parent.params(3, g=.071, td=.003, flux=.34, theta=.23,
                          field=.17, musel=.4, K=100.)
    ours = dict(p, n=2, q=1., w=1., bias=0., field_side='both')
    error = float(np.max(abs(m.bond(ours, .13)-old.bond(p, .13))))
    assert error == 0
    p = m.parameters(2, td=.03, field=.007)
    H = m.bond(p, .013)
    q = np.array([m.charge(x) for x in m.LETTERS]); qq=(q[:, None]+q[None, :]).ravel()
    marks=np.array([m.marker(x) for x in m.LETTERS])
    ml=np.repeat(marks, 11); mr=np.tile(marks, 11)
    errors = dict(hermiticity=float(np.max(abs(H-H.conj().T))),
                  charge=float(np.max(abs(H*(qq[:, None]-qq[None, :])))),
                  left_marker=float(np.max(abs(H*(ml[:, None]!=ml[None, :])))),
                  right_marker=float(np.max(abs(H*(mr[:, None]!=mr[None, :])))))
    assert max(errors.values()) < 1e-14
    return dict(alphabet=[list(x) for x in m.LETTERS], d=11, original_bond_difference=error,
                original_nonzero_entries=int(np.count_nonzero(m.bond(ours, .13))),
                variant_nonzero_entries=int(np.count_nonzero(H)),
                variant_undirected_offdiagonal=int(np.count_nonzero(np.triu(H, 1))), **errors)


def simplicity_audit():
    def canonical(h):
        # Remove the constant and every one-site telescoping bond term.
        tensor=h.reshape(11,11,11,11)
        A=(np.einsum('ijkj->ik',tensor)-np.einsum('ijil->jl',tensor))/22
        out=h-np.trace(h)/121*np.eye(121)-np.kron(A,np.eye(11))+np.kron(np.eye(11),A)
        assert abs(np.trace(out))<1e-10
        t=out.reshape(11,11,11,11)
        assert np.max(abs(np.einsum('ijkj->ik',t)-np.einsum('ijil->jl',t)))<1e-10
        return out

    ranks={}
    for name,keys,extra in [
        ('variant',['q','w','D','K','musel','td','g','field'],{}),
        ('original',['Delta','D','K','musel','td','g','field'],
         dict(q=1.,w=1.,bias=0.,field_side='both'))]:
        p=m.parameters(2,**dict({key:0. for key in keys},**extra))
        base=m.bond(p); ops=[canonical(base)]
        ops.extend(canonical(m.bond(dict(p,**{key:1.}))-base) for key in keys)
        ops.append(canonical(m.bond(p,1.)-base))
        array=np.array([np.r_[a.real.ravel(),a.imag.ravel()] for a in ops]).T
        singular=la.svdvals(array);rank=int(np.linalg.matrix_rank(array))
        assert rank==len(ops)
        ranks[name]=dict(rule_basis_rank=rank,independent_ratios_after_common_scale=rank-1,
                         smallest_singular_value=float(singular[-1]))
    p=m.parameters(2,td=.03,field=.007)
    counts={}
    for name,td in [('loop',.03),('tree',0.)]:
        h=m.bond(dict(p,td=td),.013)
        counts[name]=dict(diagonal_matrix_units=int(np.count_nonzero(np.diag(h))),
                          directed_offdiagonal_matrix_units=int(np.count_nonzero(h-np.diag(np.diag(h)))),
                          undirected_transitions=int(np.count_nonzero(np.triu(h,1))))
    d0,d1=mp.mpf(1),mp.mpf(1004)
    for j in range(12):d0,d1=d1,1004*d1-d0
    beta=1/d1;alpha=mp.sqrt(mp.mpf(2)/14)*mp.sin(mp.pi/14)
    field=mp.power(2,-960-432-101)*alpha**4*beta**2
    return dict(quotient_basis=ranks,raw_matrix_unit_counts=counts,
                reference_n11_log2_K_over_field=str(mp.log(mp.power(2,2048)/field,2)),
                reference_n11_log10_K_over_field=str(mp.log10(mp.power(2,2048)/field)),
                note='Ranks concern independent rule coefficients, not the number of literal monomials; counts are before quotienting.')


def clock_and_readout():
    reports=[]
    for n in (1, 2, 3):
        p=m.parameters(n, D=3., td=.01)
        states,H,_=m.cell_matrix(p,n+1,full=False); index={w:i for i,w in enumerate(states)}
        ev,vec=la.eigh(m.qclock(n+1,p['q'])); v=vec[:,0]*np.sign(vec[0,0])
        target=m.qclock(n+1,p['q'])+(p['musel']-1-p['q']**2)*np.eye(n+2)
        maxerr=0.;readout=0.
        for bits in product((0,1),repeat=n):
            ws=m.vacuum_words(bits); ids=[index[w] for w in ws]
            maxerr=max(maxerr,float(np.max(abs(H[ids][:,ids].toarray()-target))))
            literal=sum(v[j]**2*sum(x==('L',1) for x in word) for j,word in enumerate(ws))
            analytic=sum(bit*sum(v[i+1:]**2) for i,bit in enumerate(bits))
            readout=max(readout,abs(literal-analytic))
        assert maxerr<1e-12 and readout<1e-14
        reports.append(dict(n=n,data_sites=n+1,cell_length=n+2,code_dimension=2**n,
                            legal_dimension=len(states),matrix_error=maxerr,
                            readout_error=readout,endpoint_visibility=float(v[0]**2*v[-1]**2),
                            E=float(ev[0]),clock_gap=float(ev[1]-ev[0])))
    return reports


def naive_projector_obstruction():
    q=.25; p=2; local=np.zeros((121,121))
    moves=[((('R',0),('R',a)),(('L',a),('R',0))) for a in (0,1)]
    moves += [((('R',0),('e',0)),(('L',0),('e',0)))]
    for source,target in moves:
        v=np.zeros(121);v[11*m.INDEX[source[0]]+m.INDEX[source[1]]]=q
        v[11*m.INDEX[target[0]]+m.INDEX[target[1]]]=-1
        local+=np.outer(v,v)
    diagonal=[];energies=[]
    for bit in (0,1):
        ws=m.vacuum_words((bit,));index={w:i for i,w in enumerate(ws)}; H=np.zeros((3,3))
        for col,w in enumerate(ws):
            for j in range(p):
                source=11*m.INDEX[w[j]]+m.INDEX[w[j+1]]
                for row in np.flatnonzero(local[:,source]):
                    z=list(w);z[j],z[j+1]=m.LETTERS[row//11],m.LETTERS[row%11]
                    if tuple(z) in index:H[index[tuple(z)],col]+=local[row,source]
        diagonal.append(np.diag(H).tolist()); energies.append(float(la.eigvalsh(H)[0]))
    assert energies[0]>0 and abs(energies[1])<1e-14
    return dict(q=q,n=1,diagonals_by_bit=diagonal,lowest_energies_by_bit=energies,
                status='payload degeneracy already fails at infinite K')


def finite_K_checks():
    reports=[]
    params=m.parameters(2,D=3.,td=.01,g=.05,K=100.)
    for p in (1,2,3,4):
        states,H,bad=m.cell_matrix(params,p)
        coo=H.tocoo();keep=bad[coo.row]==bad[coo.col]
        good=sp.coo_matrix((coo.data[keep],(coo.row[keep],coo.col[keep])),shape=H.shape).tocsr()
        eclock=la.eigvalsh(m.qclock(p,params['q']))[0]-1-params['q']**2
        ideal=params['musel']+eclock
        form=good-sp.diags((params['K']-4*params['q'])*bad)-ideal*sp.eye(len(states))
        formmin,_,rf=lowest(form);actual,_,ra=lowest(H)
        lower=ideal-max(p-1,0)*params['q']**2/(params['K']-4*params['q'])
        assert formmin[0]>-1e-9 and actual[0]>=lower-1e-9
        reports.append(dict(data_sites=p,dimension=len(states),good_form_min=float(formmin[0]),
                            energy=float(actual[0]),uniform_lower=lower,
                            eigenpair_residual=max(rf,ra)))
    splitting=[]
    for n in (1,2):
        for K in (100.,1000.):
            params=m.parameters(n,D=3.,K=K,td=.01)
            _,H,_=m.cell_matrix(params,n+1)
            energies,_,res=lowest(H,2**n+1)
            splitting.append(dict(n=n,K=K,payload_bandwidth=float(np.ptp(energies[:2**n])),
                                  gap=float(energies[-1]-energies[-2]),residual=res))
    for n in (1,2):
        values=[x['payload_bandwidth'] for x in splitting if x['n']==n]
        assert 0<values[1]<values[0]/5
    return dict(form_and_energy=reports,payload_splitting=splitting,
                note='K=100,1000 diagnose splitting; not the full selector/band schedule')


def scalar(n,p,D,w,a,td=0.):
    assert p==n+1
    S=2*p+1
    H=np.diag([a]*p+[a+D]*(p+1)).astype(complex)
    for j in range(p-1):H[j,j+1]=H[j+1,j]=-w
    for j in range(p,2*p):H[j,j+1]=H[j+1,j]=-1
    H[0,p]=H[p,0]=-1
    H[-1,p-1]=-1j*td;H[p-1,-1]=1j*td
    R=splu(sp.kronsum(sp.csc_matrix(H),sp.csc_matrix(H),format='csc'))
    source=np.zeros(S*S);source[(S-1)*S]=1
    x=R.solve(source).reshape(S,S)
    U=np.zeros((S*S,n),complex);V=U.copy()
    for j in range(n):U[j*S:(j+1)*S,j]=x[j,:];V[j::S,j]=x[:,j]
    C=U.conj().T@R.solve(V)
    A=U.conj().T@R.solve(U)+V.conj().T@R.solve(V)
    return -2*C.real,-2j*A.imag


def pair_mapping():
    reports=[]
    for n in (1,2):
        p=m.parameters(n,D=4.,w=.6,g=0.,td=.003)
        states,H0,_=m.cell_matrix(p,n+1,full=False); index={w:i for i,w in enumerate(states)}
        _,H1,_=m.cell_matrix(dict(p,g=1.),n+1,full=False);C=H1-H0
        size=len(states);dc=2**n
        ev,vec=la.eigh(m.qclock(n+1,p['q']));v=vec[:,0]*np.sign(vec[0,0])
        G=np.zeros((size,dc))
        for bits in range(dc):
            ws=m.vacuum_words(tuple((bits>>i)&1 for i in range(n)))
            for j,word in enumerate(ws):G[index[word],bits]=v[j]
        E0=p['musel']+ev[0]-1-p['q']**2
        _,_,_,port=m.parts(p)
        rows=[];cols=[];data=[]
        for il,left in enumerate(states):
            for ir,right in enumerate(states):
                col=il*size+ir;source=11*m.INDEX[left[-1]]+m.INDEX[right[0]]
                for target in np.flatnonzero(port[:,source]):
                    ll=left[:-1]+(m.LETTERS[target//11],)
                    rr=(m.LETTERS[target%11],)+right[1:]
                    if ll in index and rr in index:
                        rows.append(index[ll]*size+index[rr]);cols.append(col);data.append(port[target,source])
        V=sp.coo_matrix((data,(rows,cols)),shape=(size*size,size*size)).tocsr()
        ex=np.array([any(m.excited(x) for x in w) for w in states])
        src=np.flatnonzero(np.kron(ex,ex))
        pair=sp.kronsum(H0,H0,format='csc')[src][:,src]-2*E0*sp.eye(len(src),format='csc')
        conv=sp.kronsum(C,C,format='csc')[src][:,src]
        code=np.kron(G,G)[:,[a*dc+b for b in range(dc) for a in range(dc)]]
        S=(V@code)[src];lu=splu(pair);Y=conv@lu.solve(S)
        H22=-Y.conj().T@lu.solve(Y)
        inds=[1<<j for j in range(2*n)]
        h=H22[np.ix_(inds,inds)]-H22[0,0]*np.eye(2*n)
        a=p['Delta']-(ev[0]-1-p['q']**2)
        T,M=scalar(n,n+1,p['D'],p['w'],a,p['td'])
        visibility=v[0]**2*v[-1]**2
        scale=visibility*max(la.norm(T),la.norm(M))
        terr=float(la.norm(h[:n,n:]-visibility*T)/scale)
        merr=float(la.norm(h[:n,:n]+h[n:,n:]-visibility*M)/scale)
        gaussian=gaussian_residual(H22,n)
        assert max(terr,merr,gaussian)<1e-8
        reports.append(dict(n=n,active_conversions=n,inactive_endpoint=n,
                            scalar_vertices=2*n+3,pair_legal_dimension=size**2,
                            scalar_T_error=terr,scalar_onsite_error=merr,
                            full_Fock_gaussian_residual=gaussian))
    return reports


def low_gram(p,w,a=4.,precision=False,active=None):
    # Two time orderings of the Gram integral, retaining all p low vertices.
    A=mp.eye(p)*mp.mpf(a)
    for j in range(p-1):A[j,j+1]=A[j+1,j]=-mp.mpf(w)
    R=mp.matrix(p*p)
    for i,j,k,l in product(range(p),repeat=4):
        R[i*p+j,k*p+l]=A[i,k]*(j==l)+(i==k)*A[j,l]
    R=R**-1;x=R[:,0];n=p-1 if active is None else active
    U=mp.matrix(p*p,n);V=mp.matrix(p*p,n)
    for j in range(n):
        for b in range(p):U[j*p+b,j]=x[j*p+b];V[b*p+j,j]=x[b*p+j]
    Gamma=2*U.T*R*V
    return Gamma if precision else np.array(Gamma.tolist(),float)


def gram_limit_checks():
    reports=[]
    for n in (1,2,3):
        p=n+1;w=.5;G=low_gram(p,w)
        truncated=low_gram(n,w,active=n)
        inactive_effect=float(la.norm(G-truncated)/la.norm(G))
        assert inactive_effect>1e-9
        errors=[]
        for D in (100.,500.,2500.):
            T,_=scalar(n,p,D,w,4.)
            error=float(la.norm(-D**(2*p+2)*T-G)/la.norm(G));errors.append(error)
            reports.append(dict(n=n,D=D,relative_Gram_error=error,
                                incorrect_endpoint_deletion_relative_change=inactive_effect))
        assert errors[-1]<errors[0]/10
    return reports


def simple_tree_bands():
    reports=[]
    for n in (1,2,3):
        p=n+1;w=mp.mpf(1)/4;u=mp.mpf(1)/4
        for attempt in range(10):
            q=(w*u)**2;Q=mp.matrix(p+1)
            for j in range(p+1):Q[j,j]=q*q+(j>0)
            for j in range(p):Q[j,j+1]=Q[j+1,j]=-q
            energies,vectors=mp.eigsy(Q);v=vectors[:,0]
            a=4+q*q-energies[0]
            Gamma=low_gram(p,w,a,True)
            A=mp.matrix(n);E=mp.matrix(n)
            for i,j in product(range(n),repeat=2):A[i,j]=Gamma[i,j]/w**(2*i+2*j)
            for i in range(n):E[i,i]=sum(v[j]**2 for j in range(i+1,p+1))/(q*q*w**(4*i))
            rootE=mp.diag([mp.sqrt(E[i,i]) for i in range(n)])
            invE=rootE**-1
            vals=mp.eigsy(invE*A*invE,eigvals_only=True)
            r=sorted([1/x for x in vals],reverse=True)
            if n==1 or max(r[j+1]/r[j] for j in range(n-1))<mp.mpf(1)/32:break
            u/=2
        else:raise AssertionError('small check did not terminate')
        nu=min(mp.mpf('.2')/r[0],mp.sqrt(r[-1]/(4*r[0]**3))/2)
        k=[mp.acos(nu*x/2) for x in r]
        roots=[sign*x for x in k for sign in (-1,1)]
        wrap=lambda x:abs((x+mp.pi)%(2*mp.pi)-mp.pi)
        d2=min(wrap(x-y) for i,x in enumerate(roots) for y in roots[i+1:])
        pairs=list(combinations_with_replacement(range(2*n),2))
        cooper=lambda pair:pair[0]//2==pair[1]//2 and pair[0]%2!=pair[1]%2
        nonpair=min(wrap(roots[a]+roots[b]-roots[c]-roots[d])
                    for j,(a,b) in enumerate(pairs) for c,d in pairs[j+1:]
                    if not (cooper((a,b)) and cooper((c,d))))
        margin=nu*r[-1]/4
        assert d2>=margin and nonpair>=margin
        reports.append(dict(n=n,w=str(w),u=str(u),q=str(q),nu=str(nu),
                            analytic_nonpair_margin=str(margin),
                            delta2=str(d2),delta4_excluding_Cooper_equalities=str(nonpair),
                            full_delta4='0' if n>1 else str(nonpair),
                            reciprocal_pair_sums=[str(roots[2*j]+roots[2*j+1]) for j in range(n)]))
    return reports


def source_certificates():
    base=ROOT/'codex/work_C84'
    files=['MEMBER.json','rejected_parameters.json','root_certificate_exact_n11.json','finite_budget.json']
    hashes={f:hashlib.sha256((base/f).read_bytes()).hexdigest() for f in files}
    roots=json.loads((base/files[2]).read_text());budget=json.loads((base/files[3]).read_text())
    return dict(status='read and consistency-audited; original certificate not rerun',
                hashes=hashes,root_intervals=len(roots['intervals']),
                delta2=roots['delta2'],delta4=roots['delta4'],
                finite_budget_passes=budget['sum_below_one_eighth_budget'],
                rejected=json.loads((base/files[1]).read_text()))


def main():
    start=time.time();out=dict(status='literal and finite checks, no all-n CFT claim',
                              versions=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
                              BLAS_threads=1)
    for name,run in [('source_certificates',source_certificates),('literal',literal_checks),
                     ('simplicity',simplicity_audit),
                     ('clock_readout',clock_and_readout),('naive_bias_obstruction',naive_projector_obstruction),
                     ('finite_K',finite_K_checks),('pair_mapping',pair_mapping),
                     ('Gram_limit',gram_limit_checks),('simple_tree_bands',simple_tree_bands)]:
        out[name]=run();print(name,'passed',flush=True)
        (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    out['runtime_seconds']=time.time()-start
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('all checks passed in',round(out['runtime_seconds'],2),'seconds')


if __name__=='__main__':main()
