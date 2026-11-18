"""One-worker literal-cluster, connected-coefficient and norm checks."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
from itertools import product
import json
import time
from math import log10
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import eigsh
from cluster import Cluster, fock_extract, clock_field, model
from bounds import normal_form_majorant, explicit_fourth_bound, car_weight

HERE = Path(__file__).resolve().parent


def params(n, g=0.):
    return dict(q=.21, w=.4, D=2., Delta=3., g=g, phiA=np.pi/(2*n),
                phiB=0., psiA=np.pi, psiB=np.pi/2, td=.08)


def matrix_record(A):
    return dict(real=np.asarray(A).real.tolist(), imaginary=np.asarray(A).imag.tolist())


def conversion():
    out = []
    for n in (1, 2):
        c = Cluster(n, 2, params(n))
        series = c.conversion_series(6)
        second = fock_extract(series[2], 2*n)
        fourth = fock_extract(series[4], 2*n)
        qnorm = la.norm(fourth['quartic_operator'], 2)
        assert qnorm > 1e-12
        assert la.norm(second['quartic']) < 1e-11*la.norm(series[2])
        assert la.norm(fourth['higher']) < 1e-10*qnorm
        for k in (1, 3, 5):
            assert la.norm(series[k]) < 1e-20
        zeta = (c.v[0]*c.v[-1])**2
        assert la.norm(series[4], 2) <= zeta/32
        reports = []
        for g in (.16, .08):
            full = Cluster(n, 2, params(n, g)).coefficients()['K2']
            approximation = series[0]+g*g*series[2]+g**4*series[4]
            residual = la.norm(full-approximation, 2)
            bound = zeta*g**6/(128*(1-g*g/4))
            assert residual <= bound
            reports.append(dict(g=g, residual=float(residual), sixth_order_tail_bound=bound))
        out.append(dict(n=n, literal_states=len(c.states), zeta=zeta,
                        fourth_coefficient_norm=float(la.norm(series[4], 2)),
                        quartic_operator_norm=float(qnorm), pair_indices=fourth['pair_indices'],
                        quartic_tensor=matrix_record(fourth['quartic']),
                        higher_Fock_residual=float(la.norm(fourth['higher'])), finite_g=reports))
    return out


def field_distortion():
    out = []
    for n in (2, 3):
        q = .21
        data = clock_field(n, q)
        assert np.min(data['C']) > 0
        current_error = la.norm(data['C']-data['current_C'])
        assert current_error < 1e-12
        c = Cluster(n, 1, params(n))
        exact = c.coefficients()['field2']
        expected = np.diag([-np.asarray(b)@data['C']@np.asarray(b) for b in c.bits])
        assert la.norm(exact-expected) < 1e-12
        # Endpoint-source derivative, checked against independently re-diagonalized clocks.
        derivative_errors = []
        for i in range(n):
            h = 1e-5
            base = -q*(np.eye(n+1, k=1)+np.eye(n+1, k=-1)); base[0, 0] = -1
            d = np.diag(np.arange(n+1) > i)
            _, vp = la.eigh(base+h*d); _, vm = la.eigh(base-h*d)
            finite = ((vp[0, 0]*vp[-1, 0])**2-(vm[0, 0]*vm[-1, 0])**2)/(2*h)
            derivative_errors.append(abs(finite-data['visibility_derivative'][i]))
        assert max(derivative_errors) < 1e-10
        pair_coefficients=-2*data['C'].copy(); np.fill_diagonal(pair_coefficients,0.)
        out.append(dict(n=n, C=data['C'].tolist(), density_pair_coefficients_offdiagonal=pair_coefficients.tolist(),
                        one_body_coefficients=(-np.diag(data['C'])).tolist(),
                        current_formula_error=float(current_error),
                        full_literal_field_squared_error=float(la.norm(exact-expected)),
                        visibility_derivative=data['visibility_derivative'].tolist(),
                        visibility_derivative_errors=derivative_errors))
    return out


def connected_port():
    n = 1; p = params(n, .13)
    c = Cluster(n, 3, p)
    all_terms = c.coefficients()
    left = c.coefficients({0}); right = c.coefficients({1})
    connected = all_terms['port4']-left['port4']-right['port4']
    fock = fock_extract(connected, 3)
    two = Cluster(n, 2, p)
    k_two = two.coefficients()
    embed_left = np.kron(k_two['port4'], np.eye(2))
    embed_right = np.kron(np.eye(2), k_two['port4'])
    spectator_error = max(la.norm(left['port4']-embed_left), la.norm(right['port4']-embed_right))
    assert spectator_error < 1e-12
    recursion_errors = []
    for field in (0., .2):
        formula = c.coefficients(field_scale=field)
        k2, k4 = c.bloch_fourth(field_scale=field)
        recursion_errors.append(max(la.norm(k2-formula['K2']), la.norm(k4-formula['K4'])))
    assert max(recursion_errors) < 1e-12
    zeta = (c.v[0]*c.v[-1])**2
    bound = explicit_fourth_bound(n, zeta, .2)
    assert la.norm(connected, 2) < bound['three_cell_connected_port_bound']
    assert la.norm(k_two['port4'], 2) < bound['two_cell_port_bound']
    return dict(n=n, cells=3, states=len(c.states), g=.13,
                connected_norm=float(la.norm(connected, 2)),
                connected_norm_after_scalar_removal=float(la.norm(connected-fock['scalar']*np.eye(8),2)),
                connected_scalar=matrix_record(np.array([fock['scalar']])),
                connected_quadratic=matrix_record(fock['quadratic']),
                connected_quartic=matrix_record(fock['quartic']),
                higher_Fock_norm=float(la.norm(fock['higher'])),
                spectator_error=float(spectator_error), independent_Bloch_recursion_errors=recursion_errors,
                connected_full_matrix=matrix_record(connected), bounds=bound)


def exact_low_cluster():
    c = Cluster(1, 2, params(1, .13)); chi = .2
    coefficients = c.coefficients(field_scale=chi)
    rows = []
    for lam in (.12, .06):
        H = c.H0+lam*c.V+lam*lam*chi*c.F
        energies, U = eigsh(H, k=len(c.bits), which='SA', tol=1e-13,
                            v0=np.linspace(1., 2., H.shape[0]))
        A = c.P.conj().T@U
        e, S = la.eigh(A@A.conj().T)
        polar = (S@np.diag(1/np.sqrt(e))@S.conj().T)@A
        effective = polar@np.diag(energies)@polar.conj().T
        approx = lam*lam*coefficients['K2']+lam**4*coefficients['K4']
        residual = la.norm(effective-approx, 2)
        rows.append(dict(lam=lam, residual=float(residual), scaled_sixth_order=float(residual/lam**6),
                         eigensolver_residual=float(la.norm(H@U-U*energies))))
    assert rows[1]['residual'] < rows[0]['residual']/40
    return dict(n=1, cells=2, g=.13, diagnostic_readout_scale=chi,
                status='diagnostic coefficient validation, not a certified band member', rows=rows)


def bound_constants():
    n=2; c=clock_field(n,.21); zeta=(c['v'][0]*c['v'][-1])**2
    chi=.13**2*zeta*2.**(-2*n-2)*.1/.21**2
    return dict(n=n,D=2.,g=.13,chi=chi,
                local=normal_form_majorant(n,2.,chi,kappa=0.),
                CAR_weight=car_weight(n),
                CAR=normal_form_majorant(n,2.,chi,kappa=car_weight(n)),
                fourth=explicit_fourth_bound(n,zeta,chi))


def frozen_sectors():
    """Full physical rings, including headless and zero-data competitors."""
    from collections import defaultdict
    out=[]
    for n in (1,2):
        L=n+1; p=params(n,.13); q=p['q']
        ev,_,_=model.clock(n,q)
        e_next=model.clock(n+1,q)[0]+1+q*q
        En=ev+1+q*q
        gamma_low=(1-q*q)**2*q**(2*n+4)
        lam=gamma_low/(10*np.sqrt(2))
        p.update(nu_head=1+q*q-(n+2)*En,lam=lam,hL=.02)
        W,o,_=model.bond(p)
        byin=defaultdict(list)
        for r,c in zip(*np.nonzero(W)):byin[int(c)].append((int(r),W[r,c]))
        ring=list(product(range(9),repeat=L)); indices={x:i for i,x in enumerate(ring)}
        sectors=defaultdict(list)
        labels=[]
        for idx,word in enumerate(ring):
            boundaries=tuple(int(word[j] not in (0,1) and word[(j+1)%L] not in (2,3)) for j in range(L))
            bad=tuple(int(word[j] in (0,1) and word[(j+1)%L] in (2,3)) for j in range(L))
            key=(boundaries,bad);labels.append(key);sectors[key].append(idx)
        values={}
        for key,inds in sectors.items():
            local={idx:j for j,idx in enumerate(inds)}
            H=np.zeros((len(inds),)*2,complex)
            for col,idx in enumerate(inds):
                word=ring[idx];H[col,col]+=sum(o[x] for x in word)
                for j in range(L):
                    k=(j+1)%L
                    for pair,value in byin[9*word[j]+word[k]]:
                        new=list(word);new[j]=pair//9;new[k]=pair%9
                        dest=indices[tuple(new)]
                        assert labels[dest]==key
                        H[local[dest],col]+=value
            assert la.norm(H-H.conj().T)<1e-14
            values[key]=float(la.eigvalsh(H,subset_by_index=[0,0])[0])
        selected=[value for (B,bad),value in values.items() if sum(B)==1 and not any(bad)]
        S=min(selected); minimum_excess=float('inf');count=0
        for (B,bad),value in values.items():
            ends=[i for i,b in enumerate(B) if b]
            if not ends:
                assert value>=-1e-14 and S<0
                continue
            wrong=0
            for j,end in enumerate(ends):
                previous=ends[j-1]
                length=(end-previous)%L or L
                sites=[(previous+1+t)%L for t in range(length)]
                dead=any(bad[s] for s in sites)
                wrong+=int(length!=n+1 or dead)
            if wrong:
                lower=wrong*(e_next-5*np.sqrt(2)*lam)
                assert value-S>=lower-1e-12
                minimum_excess=min(minimum_excess,value-S-lower);count+=1
        out.append(dict(n=n,L=L,physical_words=9**L,sectors=len(sectors),
                        max_sector_size=max(map(len,sectors.values())),lambda_value=lam,
                        readout=.02,selector_gamma=e_next,selected_energy=S,
                        competing_sectors_checked=count,minimum_margin_beyond_bound=minimum_excess))
    return out


def saved_band_scales():
    """Reuse saved roots; no new root search or parameter scan."""
    import importlib.util
    base=HERE.parent/'markerless_all_n_20261118'
    spec=importlib.util.spec_from_file_location('saved_band_response',base/'response.py')
    response=importlib.util.module_from_spec(spec);spec.loader.exec_module(response)
    record=json.loads((base/'checks.json').read_text())['high_precision_roots']
    member=record['checks'][0]
    n=record['n'];D=float(record['D']);q=float(record['q']);w=float(record['w'])
    epsilon=float(record['epsilon']);eta=float(member['eta'])
    ev,v,gamma=model.clock(n,q);zeta=(v[0]*v[-1])**2
    T,M=response.coefficients(n,D,w,3-ev,eta/D**n)
    T*=D**(2*n+2);M*=D**(2*n+2);F=epsilon*np.diag(gamma/q**2)
    speeds=[];residuals=[]
    for value in member['roots']:
        k=float(value);H=T*np.exp(1j*k)+T.T*np.exp(-1j*k)+M+F
        e,U=la.eigh(H);j=np.argmin(abs(e));u=U[:,j]
        speeds.append(abs(np.vdot(u,(1j*T*np.exp(1j*k)-1j*T.T*np.exp(-1j*k))@u).real))
        residuals.append(abs(e[j]))
    vmin=min(speeds);s0=zeta*D**(-2*n-2);g=2.**-100;lam=2.**-320
    chi=g*g*s0*epsilon/q**2;kappa=car_weight(n)
    nf=normal_form_majorant(n,D,chi,kappa=kappa)
    fourth=explicit_fourth_bound(n,zeta,chi,kappa=kappa)['total_local_strength_bound']
    # Twice the full bounds allow scalar removal, not counting scalars as interactions.
    log_I2=log10(2)+2*kappa/np.log(10)+(2*n+2)*log10(D)+2*log10(g)-log10(16*(1-g*g/4))
    log_I4=log10(2*fourth)+2*log10(lam)-2*log10(g)-log10(s0)
    log_R5=log10(2*nf['remainder_constant'])+3*log10(lam)-2*log10(g)-log10(s0)
    theta=min(float(member['delta2']),float(member['delta4']))/100
    assert max(log_I2,log_I4,log_R5)-log10(vmin)<log10(theta)-2
    assert lam<nf['analytic_radius']/2
    assert lam<nf['block_ground_radius']
    gamma_low=(1-q*q)**2*q**(2*n+4)
    assert lam<gamma_low/(10*np.sqrt(2))
    logW=2*log10(lam)+2*log10(g)+log10(s0)
    return dict(status='numerical scale audit of a previously saved band point; not an interval or RG certificate',
                n=n,D=D,q=q,w=w,epsilon=epsilon,eta=eta,
                normalized_speeds=speeds,min_normalized_speed=vmin,
                largest_saved_root_energy_residual=max(residuals),zeta=zeta,payload_factor_without_lambda_g=s0,
                delta2_cell_units=float(member['delta2']),delta4_cell_units=float(member['delta4']),
                physical_site_momentum_margin=float(member['delta4'])/(n+1),
                g_power_of_two=-100,lambda_power_of_two=-320,
                diagnostic_threshold=theta,
                log10_correction_over_W=dict(conversion=log_I2,fourth=log_I4,unprojected_remainder=log_R5),
                log10_W=logW,log10_smallest_physical_site_velocity=log10(n+1)+logW+log10(vmin),
                log10_physical_L1_field=2*log10(lam)+log10(chi),
                log10_largest_to_L1_coefficient_ratio=log10(D+4+q*q-(n+2)*(ev+1+q*q))-(2*log10(lam)+log10(chi)))


def main():
    start=time.monotonic()
    out=dict(workers=1, scope='finite connected coefficients and uniform local bounds; no thermodynamic CFT claim')
    for name,run in [('conversion_quartic',conversion),('readout_distortion',field_distortion),
                     ('connected_port',connected_port),('finite_lambda_check',exact_low_cluster),
                     ('majorant_constants',bound_constants),('frozen_sector_selection',frozen_sectors),
                     ('saved_band_scale_audit',saved_band_scales)]:
        out[name]=run();print(name,'passed',flush=True)
        (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    out['runtime_seconds']=time.monotonic()-start
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('All checks passed in',round(out['runtime_seconds'],2),'seconds.')


if __name__=='__main__':main()
