"""Connected pi-pair and complete integer-lattice Gaussian diagnostics.

The default run calibrates analytic controls only. It does not open candidate
MPS files or fit any V4 data. A pair exponent alone cannot certify two modes.
"""
from quartet import np, la, threadpool_limits
from pathlib import Path
import json
import time


def covariance(H):
    H=np.asarray(H,float)
    e,U=la.eigh(H)
    if min(e)<=0: raise ValueError('The Gaussian Hamiltonian must be positive.')
    eta=np.block([[np.zeros((2,2)),np.eye(2)],[np.eye(2),np.zeros((2,2))]])
    root=(U*np.sqrt(e))@U.T; inverse=(U/np.sqrt(e))@U.T
    e,U=la.eigh(root@eta@root)
    return inverse@(U*abs(e))@U.T@inverse


def lattice_stability(Q, entry_errors=(0.0,0.0,0.0)):
    """All integer (r,s), with x=2*(a*r*r+b*s*s-2*h*r*s).

    Q is the measured/derived covariance on (theta_c,phi_s), not an
    arbitrary microscopic model. entry_errors bounds |da|, |db|, |dh|.
    The eigenvalue bound makes both searches exhaustive, not a fixed cutoff.
    Amplitudes, commensurate extra operators, and Gaussian validity are
    separate inputs; 'irrelevant' below is only this conditional lattice test.
    """
    Q=np.asarray(Q,float)
    a,b,h=Q[0,0],Q[1,1],Q[0,1]
    if Q.shape!=(2,2) or np.max(abs(Q-Q.T))>1e-12:
        raise ValueError('Q must be a symmetric 2 by 2 covariance.')
    lam=float(la.eigvalsh(Q)[0])
    if lam<=0: raise ValueError('Q must be positive.')
    bound=int(np.ceil(np.sqrt(min(a,b)/lam)))+1
    best=(float('inf'),None)
    for r in range(-bound,bound+1):
        for s in range(-bound,bound+1):
            if r<0 or (r==0 and s<=0): continue
            val=a*r*r+b*s*s-2*h*r*s
            if val<best[0]-1e-12: best=(float(val),[r,s])
    da,db,dh=map(float,entry_errors)
    if min(da,db,dh)<0: raise ValueError('Error bounds must be nonnegative.')
    lam_lower=lam-float(la.norm(np.array([[da,dh],[dh,db]]),2))
    status='undetermined_with_errors'
    witnesses=[]
    if lam_lower>0:
        # Outside this box the lower quadratic form is >1 by its eigenvalue.
        radius=int(np.ceil(1/np.sqrt(lam_lower)))+1
        all_irrelevant=True
        for r in range(-radius,radius+1):
            for s in range(-radius,radius+1):
                if r<0 or (r==0 and s<=0): continue
                val=a*r*r+b*s*s-2*h*r*s
                err=da*r*r+db*s*s+2*dh*abs(r*s)
                if val-err<=1+1e-12: all_irrelevant=False
                if val+err<1-1e-12: witnesses.append([r,s])
        if witnesses: status='has_relevant_allowed_vertex'
        elif all_irrelevant: status='all_generic_vertices_strictly_irrelevant'
        elif abs(best[0]-1)<1e-12 and da+db+dh==0: status='marginal_vertex_present'
    return dict(Q=Q.tolist(),minimum_dimension=2*best[0],
                shortest_vertex=best[1],search_bound=bound,
                positive_eigenvalue=lam,error_bounds=[da,db,dh],
                status=status,strictly_relevant_witnesses=witnesses,
                quartet_dimension=2*a,pi_pair_exponent=a,
                scope='Conditional Gaussian test; no microscopic or CFT inference.')


def connected_pi_pair(r, raw, mean_left=0j, mean_right=0j, radius=4):
    """Subtract the measured disconnected part before pi demodulation.

    raw[r]=<P_0^dagger P_r>, P a specified local physical spin-pair operator.
    mean_left=<P_0>, mean_right[r]=<P_r>. All complex phases are retained.
    A binomial convolution suppresses other momenta; it cannot prove overlap
    with the desired dressed pair or separate nearby harmonics by itself.
    """
    import math
    r=np.asarray(r,int); raw=np.asarray(raw,complex)
    if len(r)!=len(raw) or not np.all(np.diff(r)==1):
        raise ValueError('Provide a consecutive, equally spaced distance grid.')
    connected=raw-np.conj(mean_left)*np.asarray(mean_right)
    demodulated=np.exp(-1j*np.pi*r)*connected
    weights=np.array([math.comb(2*radius,j) for j in range(2*radius+1)],float)/2**(2*radius)
    smooth=np.convolve(demodulated,weights,mode='valid')
    return r[radius:len(r)-radius],smooth


def pair_exponent(r, connected_projected, window):
    """A window estimator, with no claim of L/chi convergence or CFT status."""
    lo,hi=window; r=np.asarray(r); C=np.asarray(connected_projected)
    use=(r>=lo)&(r<=hi)
    if use.sum()<24 or hi/lo<2:
        raise ValueError('A short window cannot support this diagnostic.')
    x=np.log(r[use]); y=np.log(abs(C[use]))
    fit=np.polyfit(x,y,1); residual=y-np.polyval(fit,x)
    return dict(window=[lo,hi],eta_pair=float(-fit[0]),
                inferred_quartet_dimension=float(-2*fit[0]),
                log_residual_rms=float(np.sqrt(np.mean(residual**2))),
                scope='Window estimate only; verify operator overlap and L/chi convergence separately.')


def xx_correlator(distance):
    """Exact infinite half-filled ferromagnetic XX <b_0^dagger b_m>."""
    if distance==0: return 0.5
    idx=np.arange(distance)
    n=idx[:,None]-idx[None,:]+1
    G=np.zeros_like(n,dtype=float)
    odd=n%2!=0
    G[odd]=2*(-1.)**((n[odd]-1)//2)/(np.pi*n[odd])
    return float(la.det(G)/2)


def calibrate():
    H=np.eye(4); H[0,3]=H[3,0]=H[1,2]=H[2,1]=0.6
    G=covariance(H)
    assert np.max(abs(G-np.eye(4)))<1e-12
    tests=[('free_opposite_drift',G[np.ix_([2,1],[2,1])],'marginal_vertex_present',[0,1]),
           ('stable_gaussian_control',np.diag([1.5,1.2]),'all_generic_vertices_strictly_irrelevant',[0,1]),
           ('mixed_trap',np.array([[1.1,.9],[.9,1.1]]),'has_relevant_allowed_vertex',[1,1]),
           ('higher_integer_trap',np.array([[2.,6.],[6.,18.2]]),'has_relevant_allowed_vertex',[3,1])]
    stability={}
    for name,Q,status,vertex in tests:
        result=lattice_stability(Q)
        assert result['status']==status and result['shortest_vertex']==vertex
        stability[name]=result
    uncertain=lattice_stability(np.diag([1.5,1.2]),(.05,.05,.02))
    assert uncertain['status']=='all_generic_vertices_strictly_irrelevant'
    g=np.array([xx_correlator(m) for m in range(137)])
    assert np.all(np.isfinite(g)) and np.all(g>0)
    assert abs(g[1]-1/np.pi)<1e-14 and abs(g[2]-2/np.pi**2)<1e-14
    r=np.arange(1,265)
    # Two independent XX legs, with b_j -> exp(i*pi*j/2)b_j.
    # The physical adjacent inter-leg pair has momentum pi and eta=1.
    exact=(-1.)**r*g[r//2]*g[(r+1)//2]
    mean=.19+.11j
    contaminant=.05*np.exp(.31j*r)/r**.3
    raw=exact+abs(mean)**2+contaminant
    pr,pc=connected_pi_pair(r,raw,mean,mean)
    fits=[pair_exponent(pr,pc,w) for w in ((32,128),(48,160),(64,192))]
    assert max(abs(f['eta_pair']-1) for f in fits)<.01
    shifted_raw=raw+abs(.37-.21j)**2-abs(mean)**2
    _,shifted=connected_pi_pair(r,shifted_raw,.37-.21j,.37-.21j)
    assert np.max(abs(pc-shifted))<1e-14
    return dict(gaussian_controls=stability,
                robust_error_bound_control=uncertain,
                XX_two_leg_pi_pair=dict(exact_eta_pair=1,exact_quartet_dimension=2,
                    fits=fits,constant_background_invariance_error=float(np.max(abs(pc-shifted))),
                    extra_harmonic_wavevector=.31,filter_radius=4,
                    note='Range-two analytic calibration, not a proposed nearest-neighbour qubit candidate.'),
                V4_data_loaded=False)


if __name__=='__main__':
    start=time.perf_counter()
    with threadpool_limits(limits=1): out=calibrate()
    out['seconds']=time.perf_counter()-start
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
