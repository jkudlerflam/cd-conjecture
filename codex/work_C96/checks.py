"""Independent small checks of C96's norm and thermal statements."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json, math
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import quad
from tagged_kernel import bath,momenta,contact_ff
HERE=Path(__file__).resolve().parent
assert Path.cwd().resolve()==HERE


def logcosh(x):
    return np.logaddexp(x,-x)-math.log(2)


def toy(j,temp,t=1.):
    d=math.sqrt(t*t+j*j/4)
    delta=2*(logcosh(d/temp)-logcosh(t/temp))
    return -temp*(np.logaddexp(0.,delta)-math.log(2))


def covariances(h,operators,temp):
    en,un=eigh(h)
    en=en-en[0]
    mats=np.array([un.conj().T@op@un for op in operators])
    z=np.exp(-en/temp).sum()
    means=np.einsum('a,iaa->i',np.exp(-en/temp)/z,mats).real
    def corr(tau):
        weight=np.exp(-(1/temp-tau)*en[:,None]-tau*en[None,:])/z
        out=np.einsum('ab,iab,jba->ij',weight,mats,mats)
        return out-means[:,None]*means[None,:]
    return means,corr


def thermal_checks():
    hb,xs,psi,en,aa,a=bath(4,4)
    ops=[np.diag(aa[:,j]) for j in range(4)]
    eye=np.eye(2);x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]]);z=np.diag([1.,-1.])
    local=[.7*x+.2*z,.4*y-.3*z]
    matops=[np.kron(o,eye) for o in local]+[np.kron(eye,o) for o in local]
    hm=sum(matops)
    rows=[]
    for j,temp in [(.03,.2),(.03,.03),(.03,.003),(.03,.0003)]:
        mb,cb=covariances(hb.toarray(),ops,temp)
        mm,cm=covariances(j*mb[0]*hm,matops,temp)
        beta=1/temp
        endpoint=cm(0).real
        for tau in [0.,beta/100,beta/10,beta/2,beta]:
            c=cm(tau).real
            assert np.linalg.eigvalsh(c).min()>-1e-10
            assert np.linalg.eigvalsh(endpoint-c).min()>-1e-10
            assert np.linalg.eigvalsh(cb(tau).real).min()>-1e-10
        # Symmetry allows integration on half of the thermal circle.
        ch=np.empty((4,4))
        for i in range(4):
            for k in range(4):
                ch[i,k]=quad(lambda tau:cb(tau)[i,k].real,0,beta/2,epsabs=2e-12)[0]
        bconn=quad(lambda tau:np.sum(cb(tau).real*cm(tau).real),0,beta/2,epsabs=2e-12)[0]
        static=float(np.sum(ch*endpoint))
        upper=float(np.linalg.eigvalsh(ch)[-1]*np.trace(endpoint))
        assert -1e-10<=bconn<=static+1e-10<=upper+2e-10
        rows.append(dict(J=j,T=temp,B_connected=bconn,static_bound=static,operator_bound=upper))
    return rows


def thermal_determinant_checks():
    n=3;v=3;l=n+v
    states=list(momenta(n,v))
    energies=np.array([-2*np.cos(p).sum() for _,_,p in states])
    angles=np.array([2*math.pi*m/n for m,_,_ in states])
    ff=np.array([[contact_ff(p,k,n,l) for _,_,k in states]
                 for _,_,p in states])
    hb,_,_,_,aa,_=bath(n,v)
    ops=[np.diag(aa[:,i]) for i in range(n)]
    checks=[]
    for temp in [.5,.1]:
        weights=np.exp(-(energies-energies.min())/temp)
        weights/=weights.sum()
        means,cb=covariances(hb.toarray(),ops,temp)
        divided=np.empty_like(ff.real)
        for i in range(len(states)):
            for j in range(len(states)):
                de=energies[j]-energies[i]
                divided[i,j]=weights[i]/temp if abs(de)<1e-12 else (weights[i]-weights[j])/de
        row=np.array([(.5*np.sum(divided*abs(ff)**2*np.exp(1j*(angles[None,:]-angles[:,None])*r))
                       -means[0]**2/(2*temp)).real for r in range(n)])
        direct=np.array([quad(lambda tau:cb(tau)[0,r].real,0,1/(2*temp),epsabs=1e-12)[0]
                         for r in range(n)])
        error=float(np.max(abs(row-direct)))
        assert error<2e-11
        checks.append(dict(N=n,V=v,T=temp,thermal_kernel_error=error))
    return checks


def norm_checks():
    rows=[]
    # Exact periodization of 1/(1+r^2), using the elementary coth sum.
    for n in [32,64,128,256,512,1024]:
        r=np.arange(n)
        ang=2*math.pi*r/n
        # cosh(x)-cos(y)=2 sinh(x/2)^2+2 sin(y/2)^2.
        kernel=(math.pi/n*math.sinh(2*math.pi/n)
                /(2*math.sinh(math.pi/n)**2+2*np.sin(ang/2)**2))
        dist=np.minimum(r,n-r)
        lap=2*kernel-np.roll(kernel,1)-np.roll(kernel,-1)
        first=float((1+dist)@abs(kernel))
        weak=float(abs(kernel).sum()+9*(1+dist)@abs(lap))
        assert weak<36+math.pi/math.tanh(math.pi)+1e-9
        rows.append(dict(N=n,first_moment=first,weaker_norm_b3=weak))
    assert rows[-1]['first_moment']>rows[0]['first_moment']+6
    # A fixed wrong-length pair amid selected period-three scalar cells.
    b=3; n=36; pattern=np.array([.7,-.5,.1]); means_good=np.tile(pattern,n//b)
    lengths=[3,3,2,4,3,3,3,3,3,3,3,3]
    assert sum(lengths)==n
    means=[]
    for ell in lengths:
        means.extend(pattern if ell==3 else np.linspace(-.6,.8,ell))
    means=np.array(means);r=np.arange(n);dist=np.minimum(r,n-r)
    k=1/(1+dist**2);lap=2*k-np.roll(k,1)-np.roll(k,-1);c=.23
    chi=k-c/n
    matrix=np.array([[chi[(jj-ii)%n] for jj in range(n)] for ii in range(n)])
    assert np.linalg.eigvalsh(matrix).min()>0
    delta=abs(means@matrix@means-means_good@matrix@means_good)
    bound=32*b*6*(abs(k).sum()+b*b*(1+dist)@abs(lap)+c)
    assert delta<bound
    return rows,dict(scalar_difference=float(delta),proved_upper_bound=float(bound))


def main():
    thermal=thermal_checks();norms,comparison=norm_checks()
    thermal_spectral=thermal_determinant_checks()
    diag=[]
    for ratio in [10,1,.1,.01,.001]:
        j=.02;temp=ratio*j*j
        r=toy(j,temp);r2=-j*j*math.tanh(1/temp)/8
        assert -j*j/4-1e-12<=r<=1e-12
        diag.append(dict(T_over_J2=ratio,R_over_J2=r/j**2,R2_over_J2=r2/j**2))
    # Check the fourth-order coefficient at a fixed temperature where subtraction is stable.
    temp=.2
    j=.02
    r2=-j*j*math.tanh(1/temp)/8
    quartic=(math.tanh(1/temp)-1/temp)/128
    difference=(toy(j,temp)-r2)/j**4
    assert abs(difference-quartic)<2e-6
    out=dict(thermal=thermal,thermal_spectral=thermal_spectral,norms=norms,comparison=comparison,diagnostic=diag,
             fourth_order_coefficient=quartic,fourth_order_check=difference)
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    tagged=json.loads((HERE/'tagged_checks.json').read_text())
    lines=['# C96 checks','','[Numerical] These are finite checks of proved formulas and bounds. No uniform tagged-kernel estimate or full d=13 trace is inferred.','',
           '| N | v | Bath dimension | Kernel discrepancy | First moment after c_dens subtraction |',
           '| :-- | :-- | :-- | :-- | :-- |']
    for row in tagged:
        lines.append(f"| {row['N']} | {row['V']} | {row['dimension']} | {row['row_error']:.4g} | {row['K']:.12g} |")
    lines+=['','[Numerical] Complete determinant basis, eigenvalue, and form-factor checks are included through N=4. The largest full-basis error is '+f"{max(row.get('eigen_residual',0) for row in tagged):.4g}."+' At every checked size the contact spectral weight agrees with the projector identity.','',
            '| J | T | Resummed connected coefficient | Static covariance upper bound | Operator upper bound |',
            '| :-- | :-- | :-- | :-- | :-- |']
    for row in thermal:
        lines.append(f"| {row['J']} | {row['T']} | {row['B_connected']:.12g} | {row['static_bound']:.12g} | {row['operator_bound']:.12g} |")
    lines+=['','[Numerical] The material sources include a complex Pauli matrix and do not commute. Matrix positivity and the imaginary-time covariance inequality pass at five times for each row. The independent thermal determinant formula (33a) agrees with direct bath correlations at N=v=3 and T=0.5,0.1, with maximum error '+f"{max(r['thermal_kernel_error'] for r in thermal_spectral):.4g}."+'','',
            '| N | First moment of periodized 1/(1+r squared) | Weaker norm, b=3 |',
            '| :-- | :-- | :-- |']
    for row in norms:
        lines.append(f"| {row['N']} | {row['first_moment']:.10g} | {row['weaker_norm_b3']:.10g} |")
    lines+=['','[Scope] This table concerns the explicit test kernel in PROOF equation (29), not the tagged contact kernel. Its first moment diverges logarithmically while the weaker norm stays bounded. A deterministic scalar-cell comparison also satisfies (19).','',
            '| T / J squared, t=1 | Exact R / J squared | Quadratic R2 / J squared |',
            '| :-- | :-- | :-- |']
    for row in diag:
        lines.append(f"| {row['T_over_J2']} | {row['R_over_J2']:.12g} | {row['R2_over_J2']:.12g} |")
    lines+=['','[Numerical] The fourth-order coefficient at T=0.2 is '+f"{quartic:.12g}; the finite-difference value at J=0.02 is {difference:.12g}."+' The exact residual obeys its uniform bound on every tested row.','',
            '[Check correction] The first run stopped at the provisional assertion weaker_norm<36, with value 39.0379841724 at N=32. That test threshold omitted the absolute row sum. The analytic bound is 36+pi*coth(pi)<40, as proved by telescoping the second difference. The corrected check passes. No contact-kernel theorem depended on the provisional threshold.','',
            '[Reproduction] From work_C96 run python3 -B tagged_kernel.py and python3 -B checks.py. Maximum bath dimension is 6435. No parent script is executed, no large interacting calculation is run, and all outputs stay in this folder.','']
    (HERE/'CHECKS.md').write_text('\n'.join(lines))
    print('All C96 checks passed. See CHECKS.md.')

if __name__=='__main__':main()
