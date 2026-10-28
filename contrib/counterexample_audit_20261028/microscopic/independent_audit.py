"""Independent literal-cell to C33 leading-matrix audit; writes only here.

Uses bond23's literal word matrices, but constructs the two-cell resolvent,
conversion expansion and Hilbert-Schmidt quadratic projection independently.
Also regenerates the full n=23 integer certificate in a private temporary tree.
"""
from pathlib import Path
import importlib.util
import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from fractions import Fraction as Q
from math import factorial

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / 'codex/work_C33'

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

bond = load(ROOT / 'contrib/explicit_candidate/bond23.py', 'literal_bond')
lift = load(ROOT / 'contrib/counterexample_campaign/doping/bond24.py', 'literal_lift')

def quadratic(H, m):
    dim = 2**m
    cs = []
    for j in range(m):
        C = np.zeros((dim, dim))
        for bits in itertools.product((0, 1), repeat=m):
            if not bits[j]: continue
            final = list(bits); final[j] = 0
            col = sum(x << (m-1-k) for k,x in enumerate(bits))
            row = sum(x << (m-1-k) for k,x in enumerate(final))
            C[row,col] = (-1)**sum(bits[:j])
        cs.append(C)
    h = np.empty((m,m), complex)
    for i in range(m):
        for j in range(m):
            op = cs[j].T @ cs[i]
            if i == j: op -= .5*np.eye(dim)
            h[i,j] = 4*np.trace(H @ op)/dim
    scalar = np.trace(H)/dim - .5*np.trace(h)
    residual = H - scalar*np.eye(dim)
    for i in range(m):
        for j in range(m): residual -= h[i,j]*(cs[i].T @ cs[j])
    return h, float(np.linalg.norm(residual, 2))

def literal_leading(n, D, td):
    E0 = -2*np.cos(np.pi/(n+2))
    p = bond.params(Delta=4+E0,Dx=D,DM=D,td=td,gp=0,gm=0,
                    eps=0,musel=0,phx={'0':np.pi/(2*n),'+':0,'-':0},
                    phd={'0':np.pi,'+':np.pi/2,'-':np.pi/2})
    cell = bond.cell_from_bond(n,p)
    E0, vecs, residual = bond.code_vectors(cell,p)
    G = np.stack(list(vecs.values()),axis=1)
    ex = cell['exc']; ne = int(ex.sum())
    He = sp.csc_matrix(cell['H'][np.ix_(ex,ex)]-E0*np.eye(ne))
    pp = dict(p,gp=1.,gm=1.)
    hp = bond.cell_from_bond(n,pp)['H']
    C = sp.csc_matrix((hp-cell['H'])[np.ix_(ex,ex)])
    A = sp.kronsum(He,He,format='csc')
    C2 = sp.kronsum(C,C,format='csc')
    source = np.zeros((ne*ne,4**n),complex)
    for f,k,phase in [('0','0',1),('+','-',1j),('-','+',-1j)]:
        source += phase*np.kron((cell['FM'][f]@G)[ex],(cell['FX'][k]@G)[ex])
    solver = sla.splu(A)
    X = solver.solve(source)
    Hzero = -source.conj().T@X
    CX = C2@X
    Htwo = -CX.conj().T@solver.solve(CX)
    h,res = quadratic(Htwo,2*n)
    zero_scalar_error = float(np.max(abs(Hzero-np.trace(Hzero)/len(Hzero)*np.eye(len(Hzero)))))
    return h[:n,n:], h[:n,:n]+h[n:,n:], dict(
        n=n,D=D,td=td,cell_code_residual=float(residual),
        g_zero_scalar_error=zero_scalar_error,g_squared_nonquadratic_norm=res)

def small_checks(temp):
    sys.path.insert(0,str(SRC))
    integer = load(SRC/'full_loop_integer.py','integer_loop')
    records=[]
    for n,D,eta in [(1,4,12),(2,4,10)]:
        det0,det1=1,D+4
        for _ in range(2,n+2):det0,det1=det1,(D+4)*det1-det0
        beta=1/det1;td=beta*2.**(-eta)
        T,M,rec=literal_leading(n,D,td)
        old=Path.cwd();os.chdir(temp)
        try:F,MI,_,_=integer.compute(n=n,D=D,bits=128,eta_power=eta)
        finally:os.chdir(old)
        alpha=np.sqrt(2/(n+2))*np.sin(np.pi/(n+2))
        unit=2*alpha**4*beta**2
        Tf=-np.array(F,dtype=float)/(2.*2.**256)*unit
        Mf=1j*np.array(MI,dtype=float)/(2.*2.**256)*unit
        rec.update(T_relative_error=float(np.max(abs(T-Tf))/np.max(abs(Tf))),
                   M_absolute_error=float(np.max(abs(M-Mf))),
                   M_relative_error=(float(np.max(abs(M-Mf))/np.max(abs(Mf))) if np.max(abs(Mf))>1e-25 else None),
                   M_imaginary_relative_error=(float(np.max(abs(M.imag-Mf.imag))/np.max(abs(Mf))) if np.max(abs(Mf))>1e-25 else None),
                   T_literal=[[str(x) for x in row] for row in T],
                   M_literal=[[str(x) for x in row] for row in M])
        assert rec['T_relative_error']<1e-9,rec
        assert rec['M_absolute_error']<1e-16,rec
        assert rec['g_squared_nonquadratic_norm']<1e-16,rec
        records.append(rec)
    p=bond.params(Delta=2.5,Dx=7,DM=8,td=.013,gp=.031,gm=.027,
                  eps=.017,musel=1.1,K=10,lam=.023)
    h=bond.bond(p);hv=lift.vacancy_lift(h,hopping=1.3,internal_scale=.0027,chemical_potential=.017)
    ix=np.array([a*24+b for a in range(1,24) for b in range(1,24)])
    number=np.r_[0,np.ones(23)]
    hx=hv+.017/2*np.diag(np.add.outer(number,number).ravel())
    checks=bond.bond_checks(p)
    checks.update(lift_hermiticity=float(np.max(abs(hv-hv.conj().T))),
                  lift_occupied_error=float(np.max(abs(hx[np.ix_(ix,ix)]-.0027*h))))
    assert checks['lift_occupied_error']<1e-16 and checks['lift_hermiticity']==0
    for name,values in [('number',number),('charge',np.r_[0,[bond.charge(x) for x in bond.LET]]),
                        ('marker',np.r_[0,[int(x in bond.MARKERS) for x in bond.LET]])]:
        checks['lift_'+name+'_commutator']=lift.diagonal_commutator_error(hv,values)
        assert checks['lift_'+name+'_commutator']==0
    directed=0
    moves=bond.two_site_kinetic(p)+bond.port(p)
    for out,inp in zip(*np.nonzero(moves)):
        li,ri=bond.LET[inp//23],bond.LET[inp%23]
        lo,ro=bond.LET[out//23],bond.LET[out%23]
        ai,bi=bond.ARROW[li[0]],bond.ARROW[ri[0]]
        ao,bo=bond.ARROW[lo[0]],bond.ARROW[ro[0]]
        assert ai[1]==bi[0] and ao[1]==bo[0]
        assert ai[0]==ao[0] and bi[1]==bo[1]
        directed+=1
    checks['directed_pair_moves_preserving_all_bad_bonds']=directed
    V=bond.port(dict(p,lam=1.))
    checks['port_norm']=float(np.max(abs(np.linalg.eigvalsh(V))))
    assert abs(checks['port_norm']-np.sqrt(3))<1e-13
    return dict(literal_second_order=records,bond=checks)

def selector():
    def arctan(x,n):
        s=sum(((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(n)),Q(0))
        nextterm=(-1)**n*x**(2*n+1)/Q(2*n+1)
        return min(s,s+nextterm),max(s,s+nextterm)
    a,b=arctan(Q(1,5),100),arctan(Q(1,239),30)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    def cosine(x):
        s=sum(((-1)**j*x**(2*j)/factorial(2*j) for j in range(60)),Q(0))
        return s,s+x**120/factorial(120)
    def aa(l):
        low=cosine(pi[1]/(l+1))[0];high=cosine(pi[0]/(l+1))[1]
        return 2-2*high,2-2*low
    # Each interval operation below is outward, without float decisions.
    a23,a24,a25=aa(23),aa(24),aa(25)
    d23=(24*a23[0]-23*a24[1],24*a23[1]-23*a24[0])
    d24=(25*a24[0]-24*a25[1],25*a24[1]-24*a25[0])
    margin=(d23[0]-d24[1])/48
    Bupper=Q(2,2**20)
    assert margin>5*Bupper
    kappa=a24[0]-Q(1,2**136)-Q(1,2**365)
    assert kappa>0 and Q(10)>4+7*Bupper
    # gap = 2cos(pi/25)-2cos(2pi/25), independent of logical bits.
    gap=2*cosine(pi[1]/25)[0]-2*cosine(2*pi[0]/25)[1]
    assert gap>Q(47,1000)
    witness=Q((margin.numerator*2**256)//margin.denominator,2**256)
    assert witness>5*Bupper
    return dict(gamma_lower_rational=str(witness),gamma_lower=float(margin),
                kappa_lower=float(kappa),vacuum_gap_lower=float(gap),
                port_bound_upper=float(Bupper),
                wrong_length_gap_lower=float(margin-5*Bupper),
                illegal_gap_lower=float(Q(6)-7*Bupper))

def full_certificate(temp):
    work=temp/'c33';work.mkdir()
    for path in SRC.glob('*.py'):shutil.copy2(path,work/path.name)
    for name in ['derivative_n23.npz','certificate_n23.json']:shutil.copy2(SRC/name,work/name)
    logs=[]
    for name in ['full_loop_integer.py','certify_roots.py']:
        run=subprocess.run([sys.executable,name],cwd=work,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),
                           text=True,capture_output=True,check=True,timeout=240)
        logs.append(dict(script=name,stdout=run.stdout,stderr=run.stderr))
    original=np.load(SRC/'full_loop_n23.npz',allow_pickle=True)
    replay=np.load(work/'full_loop_n23.npz',allow_pickle=True)
    equal={k:bool(np.array_equal(original[k],replay[k])) for k in original.files}
    assert all(equal.values())
    cert=json.loads((work/'root_certificate_n23.json').read_text())
    exact=json.loads((work/'root_certificate_exact_n23.json').read_text())
    assert cert['root_count']==46 and cert['field_normalized_power']==-214
    assert cert['finite_g_power']==-365
    return dict(array_agreement=equal,certificate=cert,exact=exact,logs=logs)

def main():
    start=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='cd-microscopic-independent-') as td:
        temp=Path(td)
        result=small_checks(temp)
        result['selector']=selector()
        (HERE/'small_checks.json').write_text(json.dumps(result,indent=2)+'\n')
        result['c33']=full_certificate(temp)
    result['seconds']=time.monotonic()-start
    (HERE/'audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(seconds=result['seconds'],small=result['literal_second_order'],
                         c33=result['c33']['certificate']),indent=2))

if __name__=='__main__':main()
