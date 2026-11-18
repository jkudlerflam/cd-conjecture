"""Reproduce the physical projection, complete vertex, and coupled instability."""
import itertools,json,os,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from physical import parts,bond,profiles,onsite,verify_projections,CONTROLS,CHARGE,MARKER
from mechanism import scalar,pair_coefficients,normal_order,fermi_modes
from literal import literal_coefficients,finite_conversion
from flow import match,beta,complete_beta,EDGES
from exact import certificate


def encode(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x))


def main():
    assert os.environ.get('OPENBLAS_NUM_THREADS')=='1'
    started=time.perf_counter();out=dict(physical_d=9,logical_flavours=3,physical_sites_per_cell=5,workers=1)
    out['projection_errors']=verify_projections()
    out['local_checks']=[];q=(CHARGE[:,None]+CHARGE[None,:]).ravel();ml=np.repeat(MARKER,9);mr=np.tile(MARKER,9)
    for name,H in parts().items():
        errors=[np.max(abs(H-H.T)),np.max(abs(H*(q[:,None]-q[None,:]))),np.max(abs(H*(ml[:,None]-ml[None,:]))),np.max(abs(H*(mr[:,None]-mr[None,:])))]
        assert max(errors)==0
        out['local_checks'].append(dict(part=name,hermiticity_charge_marker_error=max(errors),nonzero_entries=int(np.count_nonzero(H))))
    c=pair_coefficients(3);literal,dimensions=literal_coefficients(3);T=scalar(3)[2]
    errors=[np.max(abs(c[k]-literal[k])) for k in c]
    q2=normal_order(c[2],6);q4=normal_order(c[4],6)
    errors.extend([q2['gaussian_residual'],q4['higher_residual'],np.max(abs(q2['h'][:3,3:]-T)),np.max(abs(q2['h'][:3,:3]+q2['h'][3:,3:]))])
    assert max(errors)<1e-12
    out.update(dimensions=dimensions,max_coefficient_matching_error=max(errors),quartic_W=q4['W'],quartic_pairs=q4['pairs'],K4_quadratic=q4['h'],K4_norm=np.linalg.norm(c[4],2))
    exact,Te,Me,poly=certificate();out['exact']=exact
    assert np.max(abs(Te.astype(float)-T))<1e-15 and np.max(abs(Me.astype(float)-onsite()))<1e-15
    primary=match(T,onsite(),q4);out['primary']=primary
    # Compare the cancellation-sensitive determinant interpolation to exact coefficients.
    exact_roots=np.polynomial.polynomial.polyroots([float(x) for x in poly])
    root_error=np.max(abs(exact_roots-primary['cos_roots']));assert root_error<1e-7
    out['root_vs_exact_polynomial_error']=root_error
    assert max(primary['D'])<0 and np.prod(primary['c'])<0
    kappa=-max(primary['D']);xmin=min(abs(primary['c']))
    out['instability_cone']=dict(minimum_negative_D_margin=kappa,minimum_transfer_magnitude=xmin,
        transfer_sign_product=float(np.prod(primary['c'])),formal_one_loop_blowup_time_upper_bound=np.log1p(kappa/xmin)/kappa,
        interpretation='Analytic invariant-cone instability of the leading flow; no inference about the ultimate microscopic central charge.')
    # Check the reduced six-variable flow against all forward and transfer equations.
    v=primary['velocities'];C=np.zeros((3,3))
    for x,(i,j) in zip(primary['c'],EDGES):C[i,j]=C[j,i]=x
    a=4*v[:,None]*v[None,:]/(v[:,None]+v[None,:])**2
    db,df,dc=complete_beta(primary['b'],primary['f'],C,a)
    dD=db[:,None]+db[None,:]-2*df
    combined=np.r_[[dD[i,j] for i,j in EDGES],[dc[i,j] for i,j in EDGES]]
    reduced=beta(0,np.r_[primary['D'],primary['c']],primary['a'])
    out['full_vs_reduced_flow_error']=np.max(abs(combined-reduced));assert np.max(abs(combined-reduced))<1e-9
    # Stop while physical couplings can still be chosen uniformly weak.
    initial=np.r_[primary['D'],primary['c']];limit=10*max(abs(initial))
    def stop(t,y):return limit-max(abs(y))
    stop.terminal=True
    sol=solve_ivp(lambda t,y:beta(t,y,primary['a']),(0,.1),initial,rtol=1e-10,atol=1e-11,events=stop,dense_output=True)
    assert sol.success and sol.t_events[0].size==1
    ts=np.linspace(0,sol.t[-1],21)
    out['flow_trace']=dict(time_variable='tau=g² ell',tau=ts,values=sol.sol(ts).T,stop_reason='Tenfold growth, before using the one-loop equations at strong physical coupling.')
    # A box certificate is exact for roots; corner RG values are finite numerical checks.
    width=float(1e-9);corners=[]
    for signs in itertools.product((-1,1),repeat=5):
        r=match(T,onsite(CONTROLS+width*np.array(signs)),q4)
        corners.append(dict(signs=signs,min_velocity=min(r['velocities']),negative_D_margin=-max(r['D']),transfer_sign_product=np.prod(r['c'])))
    assert min(r['negative_D_margin'] for r in corners)>0 and max(r['transfer_sign_product'] for r in corners)<0
    out['control_box_corner_checks']=corners
    # Finite conversion errors: a direct floating check at visible g plus a rigorous bound.
    delta=2-2*np.cos(np.pi/6);out['finite_conversion']=[]
    for g in (.05,.025):
        err=np.linalg.norm(finite_conversion(3,g)-c[0]-g*g*c[2]-g**4*c[4],2)
        bound=(g/delta)**6/(delta*(1-(g/delta)**2));assert err<bound
        out['finite_conversion'].append(dict(g=g,error=err,bound=bound))
    diagnostic_g=1e-6
    tail_in_quadratic_units=diagnostic_g**4/(delta**7*(1-(diagnostic_g/delta)**2))
    # Compression and normal ordering give ||delta W||<=6||delta K||;
    # each two-cell Bloch wedge has norm <=2, hence |delta Gamma|<=24||delta K||.
    vertex_error_coefficient=24*tail_in_quadratic_units/diagnostic_g**2
    D_error_at_fixed_bands=2*vertex_error_coefficient/(np.pi*min(v))
    assert D_error_at_fixed_bands<kappa/10
    corrected=match(T+diagnostic_g**2*q4['h'][:3,3:],onsite()+diagnostic_g**2*(q4['h'][:3,:3]+q4['h'][3:,3:]),q4)
    out['K4_quadratic_band_shift']=dict(g=diagnostic_g,
        maximum_cosine_root_change=np.max(abs(corrected['cos_roots']-primary['cos_roots'])),
        maximum_relative_velocity_change=np.max(abs(corrected['velocities']/v-1)),
        maximum_D_coefficient_change=np.max(abs(corrected['D']-primary['D'])),
        scope='Includes the generated vacuum-normal-ordered K4 bilinear only, not a self-consistent Fermi-sea shift or higher port corrections.')
    assert out['K4_quadratic_band_shift']['maximum_D_coefficient_change']<kappa/100
    out['error_comparison']=dict(g=diagnostic_g,minimum_velocity_in_t0_units=min(v),
        full_K4_in_quadratic_units=diagnostic_g**2*np.linalg.norm(c[4],2),
        conversion_tail_in_quadratic_units=tail_in_quadratic_units,
        fixed_band_vertex_error_coefficient=vertex_error_coefficient,
        fixed_band_D_error_coefficient_bound=D_error_at_fixed_bands,
        leading_negative_D_margin=kappa,
        scope='Finite two-cell conversion estimate at fixed band data; does not certify dressed roots, port errors, finite K, or global sector selection.')
    full=bond(delta,1.,100.,.05,.01,1e-7)
    assert np.max(abs(full-full.T))==np.max(abs(full*(q[:,None]-q[None,:])))==0
    assert np.max(abs(full*(ml[:,None]-ml[None,:])))==np.max(abs(full*(mr[:,None]-mr[None,:])))==0
    out['seconds']=time.perf_counter()-started
    path=Path(__file__).parent;(path/'audit.json').write_text(json.dumps(out,indent=2,default=encode)+'\n')
    (path/'exact.json').write_text(json.dumps(exact,indent=2)+'\n')
    print(json.dumps(dict(seconds=out['seconds'],matching_error=max(errors),roots=primary['cos_roots'],velocities=v,
        D=primary['D'],c=primary['c'],other_pair_sum_margin=primary['other_pair_sum_chord_margin'],
        instability=out['instability_cone'],error_comparison=out['error_comparison']),indent=2,default=encode))


if __name__=='__main__':main()
