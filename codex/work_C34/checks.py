"""Small finite checks for C34. NumPy/SciPy only; outputs stay here."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import json
import numpy as np
from scipy.linalg import null_space

P = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]],complex), np.diag([1,-1]).astype(complex)]

def dm(d):
    x,y,z=d
    return np.array([[0,z,-y],[-z,0,x],[y,-x,0.]])

def cert(J,B):
    s=np.linalg.norm(J,2)
    M=np.vstack((J+s*np.eye(3),J.T+s*np.eye(3)))
    E=null_space(M, rcond=1e-11)
    return float(2*np.linalg.norm(E.T@B)), E

def chain(J,B,L):
    N=1<<L
    H=np.zeros((N,N),complex)
    bits=np.arange(N)
    # Site 0 is the least significant bit. Y|0>=i|1>.
    for j in range(L):
        k=(j+1)%L
        for a in range(3):
            if B[a]:
                z=1-2*((bits>>j)&1)
                target=bits if a==2 else bits^(1<<j)
                amp=z if a==2 else (1j*z if a==1 else np.ones(N))
                H[target,bits]+=B[a]*amp
            for b in range(3):
                if J[a,b]==0: continue
                zj=1-2*((bits>>j)&1); zk=1-2*((bits>>k)&1)
                target=bits.copy()
                if a!=2: target^=1<<j
                if b!=2: target^=1<<k
                aa=zj if a==2 else (1j*zj if a==1 else np.ones(N))
                ab=zk if b==2 else (1j*zk if b==1 else np.ones(N))
                H[target,bits]+=J[a,b]*aa*ab
    assert np.max(np.abs(H-H.conj().T))<1e-12
    return H

def signs(S,D,B):
    # Only used below with diagonal S having distinct eigenvalues.
    unit=[]; anti=[]
    for a in [-1,1]:
      for b in [-1,1]:
       for c in [-1,1]:
        R=np.diag([a,b,c]); det=a*b*c
        if np.linalg.norm(R@B-B)>1e-10: continue
        if np.linalg.norm(det*R@D-D)>1e-10: continue
        (unit if det==1 else anti).append([a,b,c])
    return unit,anti

examples={
 'TR_only': (np.array([1.,.6,-.2]),np.array([.31,.23,.17]),np.zeros(3)),
 'reflection_DM': (np.array([.7,.2,-1.4]),np.array([0.,0.,.3]),np.array([.4,.6,0.])),
 'reflection_symmetric_gap': (np.array([-1.4,.3,.6]),np.zeros(3),np.array([.4,.7,0.])),
 'none_misaligned': (np.array([.7,.2,-1.4]),np.array([.11,.07,.3]),np.array([.4,.6,.5])),
 'none_aligned_gap': (np.array([.7,.2,-1.4]),np.array([0.,0.,.3]),np.array([.4,.6,.5])),
 'none_symmetric_gap': (np.array([-1.4,.3,.6]),np.zeros(3),np.array([.4,.7,.5])),
 'none_aligned_fail': (np.array([.7,.2,1.4]),np.array([0.,0.,.3]),np.array([.4,.6,.5])),
}
result={'examples':{}}
for name,(ev,d,b) in examples.items():
    J=np.diag(ev)+dm(d)
    lower,E=cert(J,b)
    unit,anti=signs(np.diag(ev),d,b)
    row={'unitary_actions':unit,'antiunitary_actions':anti,'certificate':lower,'axis_dim':E.shape[1]}
    if lower>0:
        row['gaps']={}
        for L in [4,6,8]:
            vals=np.linalg.eigvalsh(chain(J,b,L))
            gap=float(vals[1]-vals[0]); assert gap>=lower-1e-10
            row['gaps'][L]=gap
    result['examples'][name]=row

# Check degenerate maximal singular spaces and the averaging proof.
# J is minus an orthogonal rotation about z, admitting rotating local axes.
L=6; q=2*np.pi/L
R=np.array([[np.cos(q),-np.sin(q),0],[np.sin(q),np.cos(q),0],[0,0,1.]])
J=-R; n0=np.array([.6,0,.8]); B=np.array([.1,.2,.7])
axes=np.array([np.linalg.matrix_power(R,-j)@n0 for j in range(L)])
err=max(np.linalg.norm(J@axes[(j+1)%L]+axes[j]) for j in range(L))
err=max(err,max(np.linalg.norm(J.T@axes[j]+axes[(j+1)%L]) for j in range(L)))
mean=axes.mean(axis=0); lower,E=cert(J,B)
assert err<1e-12 and np.min(axes@B)>0
assert np.linalg.norm(J@mean+mean)<1e-12
assert lower>=2*np.min(axes@B)-1e-12
result['rotating_axes']={'L':L,'bond_error':float(err),'mean_axis':mean.tolist(),'rotating_bound':float(2*np.min(axes@B)),'uniform_bound':lower}

# Norm formula for axial transverse block.
rng=np.random.default_rng(34028); maxerr=0.
for _ in range(200):
    x,y,d=rng.normal(size=3)
    T=np.array([[x,d],[-d,y]])
    formula=np.sqrt(((x+y)/2)**2+d*d)+abs((x-y)/2)
    maxerr=max(maxerr,abs(formula-np.linalg.norm(T,2)))
assert maxerr<1e-12
result['block_norm_error']=maxerr

# Exact product eigenstate conditions, including a TR-only uniaxial case.
S=np.diag([.4,.4,-1.1]); J=S+dm([.3,0,.2]); B=np.zeros(3)
H=chain(J,B,6); vac=np.zeros(64); vac[0]=1
leak=np.linalg.norm(H@vac-H[0,0]*vac)
assert leak<1e-12
result['TR_product_eigenstate']={'L':6,'leakage':float(leak),'eigenvalue':float(H[0,0].real),'ground_energy':float(np.linalg.eigvalsh(H)[0])}

# Uniform product-state energy and analytic vacuum-preemption threshold.
I=np.eye(2); n=np.diag([0.,1.]); ket=np.eye(4)
def bond(phi,V,mu,gL,gR):
    h=V*np.kron(n,n)-mu*(np.kron(n,I)+np.kron(I,n))/2
    h=h.astype(complex)
    for out,inn,amp in [(2,1,-np.exp(1j*phi)),(3,1,gL),(3,2,gR)]:
        h[out,inn]+=amp; h[inn,out]+=np.conjugate(amp)
    return h
maxenergy=0.; maxmin=0.
for _ in range(100):
    phi=rng.uniform(-np.pi,np.pi); V=rng.uniform(-3,3); mu=rng.normal()
    gL=rng.normal()+1j*rng.normal(); gR=rng.normal()+1j*rng.normal()
    G=abs(gL+gR); A=V+2*np.cos(phi)
    root=np.sqrt(A*A+4*G*G); density=(1-A/root)/2
    threshold=-2*np.cos(phi)+(A-root)/2
    # Matrix convention gives <11|h|01>=gL, so phase=arg(gsum)+pi.
    phase=np.angle(gL+gR)+np.pi
    state=np.array([np.sqrt(1-density),np.exp(1j*phase)*np.sqrt(density)])
    pair=np.kron(state,state)
    exact=float(np.vdot(pair,bond(phi,V,mu,gL,gR)@pair).real)
    formula=(-2*np.cos(phi)-mu)*density+A*density*density-2*G*density**1.5*np.sqrt(1-density)
    maxenergy=max(maxenergy,abs(exact-formula))
    maxmin=max(maxmin,abs(exact/density-(threshold-mu)))
assert maxenergy<1e-11 and maxmin<1e-11
result['product_witness']={'random_bonds':100,'two_site_energy_error':maxenergy,'threshold_error':maxmin}

# New domain-wall bound, beyond the common-axis gap criterion.
wall_tests=[]
for a,c,dz,d,b in [(.4,-1.4,.2,.3,0.),(.4,-1.4,.2,.3,.23),(.4,-1.4,.2,.3,-.23),(.4,-1.4,0.,.3,.23),(-.3,-1.6,-.4,.25,0.),(.4,-(np.hypot(.4,.2)+.6),.2,.3,.23)]:
    J=np.diag([a,a,c])+dm([d,0,dz]); B=np.array([0.,0.,b])
    kappa=abs(c)-np.hypot(a,dz)-2*d
    assert kappa>=-1e-12 and (b!=0 or kappa>0)
    source_bound,_=cert(J,B); assert source_bound<1e-12
    for L in [4,6,8]:
        H=chain(J,B,L); bits=np.arange(1<<L)
        walls=np.zeros(1<<L)
        for j in range(L): walls+=((bits>>j)&1)!=((bits>>((j+1)%L))&1)
        zsum=L-2*np.array([int(i).bit_count() for i in bits])
        Eref=c*L-abs(b)*L
        lower_diag=2*kappa*walls + b*zsum+abs(b)*L
        remainder=H-Eref*np.eye(1<<L)-np.diag(lower_diag)
        mineig=float(np.linalg.eigvalsh(remainder)[0]); assert mineig>-1e-10
        energies=np.linalg.eigvalsh(H)
        ground_dim=2 if b==0 else 1
        gap=float(energies[ground_dim]-energies[0])
        bound=4*kappa if b==0 else 2*abs(b)
        assert abs(energies[0]-Eref)<1e-10 and gap>=bound-1e-10
        wall_tests.append({'parameters':[a,c,dz,d,b],'L':L,'operator_remainder_min':mineig,'gap':gap,'bound':bound})
result['wall_bound']=wall_tests

with open('checks_results.json','w') as f: json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
