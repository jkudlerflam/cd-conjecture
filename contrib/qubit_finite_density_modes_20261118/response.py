"""Full asymmetric covariance/response interface; no p153c fit is invented.

Phi=(phi_A,phi_B,theta_A,theta_B), covariance of differences = G log r.
X=(-dphi_A,-dphi_B,dtheta_A,dtheta_B)/pi; h couples as -h dot X.
"""
import numpy as np
from scipy import linalg

ETA=np.block([[np.zeros((2,2)),np.eye(2)],[np.eye(2),np.zeros((2,2))]])
SIGNS=np.diag([-1.,-1.,1.,1.])


def positive_power(matrix,power):
    values,vectors=linalg.eigh(matrix)
    if min(values)<=0:raise ValueError('quadratic response is not positive definite')
    return (vectors*values**power)@vectors.T


def gaussian(H):
    H=np.asarray(H,float)
    if not np.allclose(H,H.T,atol=1e-12):raise ValueError('nonsymmetric H')
    half=positive_power(H,.5); inverse=positive_power(H,-.5)
    eigen,vectors=linalg.eigh(half@ETA@half)
    G=inverse@(vectors*abs(eigen))@vectors.T@inverse
    return dict(H=H,G=G,signed_velocities=eigen,
                canonical_error=float(np.max(abs(G@ETA@G-ETA))))


def from_response(chi,matching=None,contact=None):
    """chi must be the bulk response with UV contact and operator matching.

    If O=matching X, then chi_O=matching chi_X matching^T+contact.
    A bare microscopic susceptibility cannot be passed as chi_X silently.
    """
    chi=np.asarray(chi,float)
    if contact is not None:chi=chi-np.asarray(contact,float)
    if matching is not None:
        Z=np.asarray(matching,float)
        if np.linalg.matrix_rank(Z)!=4:raise ValueError('probe matching has rank below four')
        inverse=np.linalg.inv(Z);chi=inverse@chi@inverse.T
    positive_power(chi,1)
    return gaussian(SIGNS@np.linalg.inv(chi)@SIGNS/np.pi)


def dimension(G,n,m):
    ell=np.r_[-2*np.asarray(m,float),np.asarray(n,float)]
    return float(ell@G@ell/4)


def mismatch(k,rho,n,m):
    return float(np.angle(np.exp(1j*(np.dot(n,k)+2*np.pi*np.dot(m,rho)))))


def dangerous_vertices(G,k,rho,xmax=2.,G_error_norm=0.,dk=(0.,0.),drho=(0.,0.)):
    """Exhaust every parity-even vertex potentially below xmax.

    Uses a proven Euclidean radius from a lower eigenvalue bound. Input
    uncertainty is a spectral norm bound on G, not inferred from old data.
    """
    lower=float(min(np.linalg.eigvalsh(G))-G_error_norm)
    if lower<=0:raise ValueError('covariance error leaves no positive lower eigenvalue bound')
    radius2=4*xmax/lower; lim=int(np.floor(np.sqrt(radius2)))
    out=[]
    import itertools
    for na,nb,ma,mb in itertools.product(range(-lim,lim+1),repeat=4):
        if (na+nb)%2:continue
        v=(na,nb,ma,mb)
        first=next((x for x in v if x),0)
        if first<=0:continue   # remove zero and Hermitian-conjugate duplicates
        norm2=na*na+nb*nb+4*(ma*ma+mb*mb)
        if norm2>radius2:continue
        n=(na,nb);m=(ma,mb);x=dimension(G,n,m);dx=G_error_norm*norm2/4
        if x-dx>xmax:continue
        q=mismatch(k,rho,n,m)
        dq=float(np.dot(abs(np.array(n)),dk)+2*np.pi*np.dot(abs(np.array(m)),drho))
        ell=np.r_[-2*np.array(m),np.array(n)]
        out.append(dict(n=n,m=m,x=x,dx=dx,q=q,dq=dq,
                        mismatch_lower_bound=max(0.,abs(q)-dq),
                        self_null=bool(ell@ETA@ell==0)))
    return sorted(out,key=lambda row:(row['x'],abs(row['q'])))


def response_error_bound(chi,error_norm):
    """Operator norm bound for H inferred from a matched susceptibility."""
    smallest=float(min(np.linalg.eigvalsh(chi)))
    if error_norm>=smallest:raise ValueError('response uncertainty reaches singularity')
    return float(error_norm/(np.pi*smallest*(smallest-error_norm)))


def gaussian_error_bounds(H,H_error_norm):
    """Conservative spectral-norm bounds, conditional on the supplied H error.

    Matrix square-root inequalities propagate the error through the actual
    chiral formula. Floating point evaluation is not interval arithmetic.
    """
    m,M=np.linalg.eigvalsh(H)[[0,-1]]; e=float(H_error_norm)
    if not 0<=e<m:raise ValueError('H error must preserve strict positivity')
    droot=e/(np.sqrt(m)+np.sqrt(m-e))
    dinverse=droot/np.sqrt(m*(m-e))
    dB=droot*(np.sqrt(M)+np.sqrt(M+e))
    dabs=(2*M+e)*dB/(2*m-e)
    dG=dinverse*((M+e)/np.sqrt(m-e)+M/np.sqrt(m))+dabs/np.sqrt(m*(m-e))
    return dict(H_minimum_lower_bound=float(m-e),
                covariance_error_norm=float(dG),signed_velocity_error=float(dB))
