"""C110 checks. Small matrices only, no computation above 256 states."""
from pathlib import Path
import numpy as np

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])

def hbond(s):
    return -0.5*(np.kron(Z,I)+np.kron(I,Z))+s*(np.kron(X,Z)-np.kron(Z,X))

def ring(L,s,w):
    H=np.zeros((2**L,2**L))
    for state in range(2**L):
        z=np.array([1-2*((state>>j)&1) for j in range(L)])
        for j in range(L):
            k=(j+1)%L
            H[state,state]-=0.5*w[j]*(z[j]+z[k])
            H[state^(1<<j),state]+=s*w[j]*z[k]
            H[state^(1<<k),state]-=s*w[j]*z[j]
    assert np.max(np.abs(H-H.T))<1e-13
    return H

out=['# C110 checks','','[Numerical] All matrices have dimension at most 256. These checks test the explicit derivations, not the general bound.','']
s=0.2
h=hbond(s).reshape(2,2,2,2)
pt_left=np.trace(h,axis1=0,axis2=2)
pt_right=np.trace(h,axis1=1,axis2=3)
assert np.max(np.abs(pt_left-pt_right))<1e-14
assert abs(np.trace(hbond(s)))<1e-14
out+=['[Numerical] The interacting gapped bond has trace zero and equal partial traces, so it lies in the Hilbert-Schmidt canonical gauge.','','| L | uniform energy | gap | rigorous gap lower bound |','| :-- | --: | --: | --: |']
for L in [4,6,8]:
    eig=np.linalg.eigvalsh(ring(L,s,np.ones(L)))
    assert abs(eig[0]+L)<1e-11
    assert eig[1]-eig[0]>=2-4*s-1e-11
    out.append(f'| {L} | {eig[0]:.10f} | {eig[1]-eig[0]:.10f} | {2-4*s:.10f} |')
out+=['','| L | L times exact weighted deficit | L times trial deficit | rigorous upper bound times L |','| :-- | --: | --: | --: |']
amp=0.4
for L in [4,6,8]:
    w=1+amp*np.cos(2*np.pi*np.arange(L)/L)
    b=(w+np.roll(w,1))/2
    delta=w-np.roll(w,1)
    x=-s*delta/b
    z=np.sqrt(1-x*x)
    trial=-np.sum(b*z)+s*np.sum(w*(x*np.roll(z,-1)-z*np.roll(x,-1)))
    deficit=-np.sum(w)-np.linalg.eigvalsh(ring(L,s,w))[0]
    trial_deficit=-np.sum(w)-trial
    upper=s*s/(2-4*s)*np.sum(delta*delta/b)
    assert deficit>=trial_deficit-1e-11
    assert deficit<=upper+1e-11
    assert deficit>0
    out.append(f'| {L} | {L*deficit:.10f} | {L*trial_deficit:.10f} | {L*upper:.10f} |')
I_f=4*np.pi**2*(1-np.sqrt(1-amp**2))
out+=['',f'[Proved formula, numerical evaluation] For f(x)=1+0.4 cos(2 pi x), the integral of f prime squared divided by f is {I_f:.12f}. The rigorous asymptotic lower and upper coefficients are {s*s*I_f/2:.12f} and {s*s*I_f/(2-4*s):.12f}.','']
# Product-state formulas only beyond eight sites.
for L in [64,128,256,512]:
    w=1+amp*np.cos(2*np.pi*np.arange(L)/L)
    b=(w+np.roll(w,1))/2
    delta=w-np.roll(w,1)
    x=-s*delta/b
    z=np.sqrt(1-x*x)
    # Stable evaluation of 1-z avoids cancellation.
    trial_deficit=-np.sum(b*x*x/(1+z))-s*np.sum(w*(x*np.roll(z,-1)-z*np.roll(x,-1)))
    assert abs(L*trial_deficit-s*s*I_f/2)<0.001
    gauge_deficit=np.sum(s*s*delta*delta/(np.sqrt(b*b+s*s*delta*delta)+b))
    assert abs(L*gauge_deficit-s*s*I_f/2)<0.001
out+=['[Numerical] Scalar sums at L=64,128,256,512 verify the product trial coefficient and the exactly solvable transverse-coboundary coefficient. No large spin matrices were formed.','']
# Spin XX with flux. Check particle-number sector formula against small matrices.
out+=['| L | flux | XX exact energy | finite-sum energy |','| :-- | --: | --: | --: |']
for L in [4,6,8]:
    for phi in [0.,np.pi/3,np.pi]:
        H=np.zeros((2**L,2**L),dtype=complex)
        for state in range(2**L):
            for j in range(L):
                k=(j+1)%L
                a,b=(state>>j)&1,(state>>k)&1
                if a!=b:
                    phase=1 if j<L-1 else np.exp(1j*phi*(1 if a else -1))
                    H[state^(1<<j)^(1<<k),state]-=phase
        exact=np.linalg.eigvalsh(H)[0]
        predicted=-2*np.cos(phi/L)/np.sin(np.pi/L)
        assert abs(exact-predicted)<1e-11
        out.append(f'| {L} | {phi:.8f} | {exact:.10f} | {predicted:.10f} |')
out+=['','[Numerical] Nine XX flux matrices agree with the even-length formula. At flux pi, Q_L approaches -4 for hopping t=1.']
for L in [64,128,256]:
    Q=6/np.pi*L*(-2*L/np.pi+2/np.tan(np.pi/L))
    assert abs(Q+4)<0.001
# Twist gauge-orbit formula, using only a two-state local spectrum.
for r in [0.,0.25,1.]:
    a,b=-0.5+r,-0.5-r
    maxnorm=abs(a)+abs(b)
    observed=max(np.linalg.norm(a*Z+b*(np.cos(theta)*Z+np.sin(theta)*X),2) for theta in np.linspace(0,np.pi,101))
    assert abs(observed-maxnorm)<1e-12
out+=['','[Numerical] The one-site twist optimization agrees with max(1,2|r|), including the gauge-dependent linear divergence at r=1.','','[Work record] All checks passed. No heavy numerics, network shell commands, or cluster jobs were used.','']
Path('CHECKS.md').write_text('\n'.join(out))
print('All C110 checks passed.')
