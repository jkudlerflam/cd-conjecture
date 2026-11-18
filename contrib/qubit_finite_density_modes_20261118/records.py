"""Read-only p153c record audit and finite-window carrier measurements.

No exponent or thermodynamic phase fit; no MPS or production code is run.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import json
import hashlib
import time
import numpy as np
from scipy.optimize import root_scalar
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def wrap(x):return np.angle(np.exp(1j*x))


def carrier(rs, values, seed, lo, hi, radius=4):
    """Unweighted complex phase slope after a binomial carrier filter.

    This returns a finite-window carrier, not a dressed momentum bound.
    """
    rs=np.asarray(rs); values=np.asarray(values)
    h=values*np.exp(-1j*seed*rs)
    for _ in range(radius):h=np.convolve(h,[.25,.5,.25],mode='same')
    mask=(rs>=lo)&(rs<=hi)&(rs>=rs[0]+radius)&(rs<=rs[-1]-radius)
    x=rs[mask]; phase=np.unwrap(np.angle(h[mask]))
    if len(x)<8:raise ValueError('fewer than eight phase points')
    slope,intercept=np.polyfit(x,phase,1)
    residual=phase-slope*x-intercept
    rms=float(np.sqrt(np.mean(residual**2)))
    # This is a fit residual only; correlated/boundary/model errors are larger.
    formal=float(np.sqrt(np.sum(residual**2)/(len(x)-2)/np.sum((x-x.mean())**2)))
    return dict(k=float(wrap(seed+slope)),points=len(x),lo=int(x[0]),hi=int(x[-1]),
                phase_rms=rms,formal_slope_error=formal,
                minimum_filtered_magnitude=float(np.min(abs(h[mask]))))


def bare_band():
    phase=1.53
    e=lambda k:2*np.cos(k-phase)-2*np.cos(2*k)
    de=lambda k:-2*np.sin(k-phase)+4*np.sin(2*k)
    dde=lambda k:-2*np.cos(k-phase)+8*np.cos(2*k)
    ks=[root_scalar(de,x0=k,x1=k+.001).root for k in (-2.9,-.25)]
    return dict(phi=phase,mu_eff=(40+2*(-19.97589125)-3/40)*80,
                k=ks,energy=[e(k) for k in ks],curvature=[dde(k) for k in ks],
                pair_shift_from_minus_pi=float(wrap(sum(ks)+np.pi)),
                quartet_mismatch=float(wrap(2*sum(ks))),physical_t2=1/80)


def standard_paths():
    base=ROOT/'codex/work_C190/gpu/typhon_out'
    out=[]
    for folder,glob in (('outG','p153c_L*_P*_rand1.json'),('outG','p153c_L64_P0_prod14.json'),
                        ('outX1024','p153c_L256_P*.json'),('outX1536','p153c_L256_P*.json')):
        out.extend(sorted((base/folder).glob(glob)))
    out.extend(sorted((ROOT/'codex/work_C195/typhon_out/out512').glob('p153c_*.json')))
    out.extend(sorted((ROOT/'codex/work_C195/typhon_out_x2048').glob('p153c_*.json')))
    return out


def ordinary_record(path):
    data=json.loads(path.read_text()); cfg=data['cfg']; m=data.get('meas',{})
    if not m or data['status']!='done':return None
    assert cfg['model']=='p153' and cfg['hz']==-19.97589125
    L=cfg['L']; hist=data['hist']; last=hist[-1]
    previous=[x for x in hist if x.get('alpha',0)==0 and x.get('chi',0)<last['chi'] and 'Sbar6' in x]
    prevS=previous[-1]['Sbar6'] if previous else cfg.get('continues',{}).get('Sbar6')
    result=dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                L=L,P=cfg['P'],chi=last['chi'],hz=cfg['hz'],
                energy=last.get('E_phys',last['E']+data.get('E_off',0)),
                n_mid=float(np.mean(m['dens'][L//4:3*L//4])),
                last_entropy_change=last['Sbar6']-hist[-2]['Sbar6'] if len(hist)>1 else None,
                chi_entropy_change=data['Sbar6']-prevS if prevS is not None else None,
                discarded=last['trunc'],local_residual=last['res_max'],
                gwin=m['gwin'],has_long_pair=False,has_density_response=False,has_current_response=False)
    G=np.array(m['Gb_re'])+1j*np.array(m['Gb_im'])
    rs=np.arange(1,len(G)); vals=np.array([np.mean(np.diag(G,r)) for r in rs])
    windows=[(8,24),(8,32),(12,36)] if len(G)>=48 else [(6,20)]
    result['carrier_windows']=[]
    for lo,hi in windows:
        row=[carrier(rs,vals,seed,lo,hi) for seed in (-2.9,-.24)]
        result['carrier_windows'].append(dict(window=[lo,hi],valleys=row,
             pair_shift_from_minus_pi=float(wrap(sum(r['k'] for r in row)+np.pi)),
             quartet_mismatch=float(wrap(2*sum(r['k'] for r in row)))))
    return result


def load_long(path):
    data=json.loads(path.read_text());m=data['meas']; L=m['L']
    G={}
    for key,(real,imag) in m['Gb_rows'].items():
        i=int(key)
        for r,(re,im) in enumerate(zip(real,imag),1):G[i,i+r]=re+1j*im
    return data,m,G


def long_record(path):
    data,m,G=load_long(path);L=m['L'];dens=np.array(m['dens'])
    assert data['cfg']['model']=='p153' and data['cfg']['hz']==-19.97589125
    result=dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                L=L,P=data['cfg']['P'],chi=data['Dmax'],hz=data['cfg']['hz'],
                n_mid=float(np.mean(dens[L//4:3*L//4])),pair_rows=[],centered_boson=[])
    if L>=128:
        windows=[(12,32),(16,40)] if L>=256 else [(8,20)]
        for key,row in m['Dpair'].items():
            d,i=map(int,key.split('_'));rs=np.arange(len(row))+d+1
            raw=np.array([complex(*v) for v in row]); mean_left=complex(*m['mm'][key])
            values=raw-mean_left*np.array([complex(*m['pp'][str(d)][str(i+r)]) for r in rs])
            for lo,hi in windows:
                estimate=carrier(rs,values,-np.pi,lo,hi)
                estimate.update(d=d,anchor=i,shift_from_minus_pi=float(wrap(estimate['k']+np.pi)),
                    edge_clearance=int(min(i,L-1-(i+hi+d))),
                    clearance_over_rmax=float(min(i,L-1-(i+hi+d))/hi))
                result['pair_rows'].append(estimate)
        rs=[];values=[]
        for r in range(1,L//2):
            entries=[v for (i,j),v in G.items() if j-i==r and abs(i+j-(L-1))<=8]
            if len(entries)<3:break
            rs.append(r);values.append(np.mean(entries))
        bwin=[(16,40),(24,64)] if L>=256 else [(8,24),(12,32)]
        for lo,hi in bwin:
            row=[carrier(rs,values,seed,lo,hi) for seed in (-2.9,-.24)]
            result['centered_boson'].append(dict(window=[lo,hi],valleys=row,
                pair_shift_from_minus_pi=float(wrap(sum(r['k'] for r in row)+np.pi))))
    # A bulk representative of the derived second-order density dressing.
    local=[]
    for key,connected in m['Cnn'].items():
        i=int(key)
        D2=connected[1]+dens[i]*dens[i+2]
        D3=connected[2]+dens[i]*dens[i+3]
        B=1-3.5*dens[i]+1.5*D2+11/18*D3+.5*G[i,i+2].real
        local.append(dict(anchor=i,n=float(dens[i]),D2=float(D2),D3=float(D3),
                          real_T2=float(G[i,i+2].real),B_density=float(B),
                          delta_n_second_order=float(2/40**2*B)))
    result['density_dressing_bulk_representatives']=local
    result['mean_dressing_representative']=float(np.mean([r['delta_n_second_order'] for r in local]))
    result['leading_dressing_upper_bound_at_nmid']=float(2/40**2*(1-8/9*result['n_mid']))
    return result


def export_interface(out):
    strongest=max(out['standard_records'],key=lambda r:(r['L'],r['chi']))
    w=np.exp(1.53j)/80
    contract=dict(
        date='2026-11-18',status='undecided; no matched two-boson medium parameters yet',
        physical_point=dict(microscopic_d=2,range=1,Jx=(1+w.real)/2,Jy=(-1+w.real)/2,
                            D=w.imag/2,Jz=10,hz=-19.97589125,
                            end_convention='same hz on every open-chain site',
                            Delta=1,V=40,abs_t2=1/80,mu=40+2*(-19.97589125)),
        convention=dict(fields=['phi_A','phi_B','theta_A','theta_B'],
            X='(-dphi_A,-dphi_B,dtheta_A,dtheta_B)/pi',
            difference_covariance='G log(|x|/a0)',vertex='exp i(n.theta-2m.phi)',
            dimension='ell.T G ell / 4; ell=(-2m_A,-2m_B,n_A,n_B)',
            mismatch='wrap(n.k + 2pi m.rho)',momentum_unit='radians per physical site',
            physical_energy_unit='Delta=1',effective_energy_unit='abs_t2=1/80',
            signed_speed='eigenvalues of eta H; omega/q has opposite sign for exp i(qx-omega t)'),
        matched_medium=dict(rho=None,rho_errors=None,k=None,k_errors=None,
            chi_X=None,chi_error_norm=None,probe_matching_Z=None,UV_contact=None,
            H=None,H_error_norm=None,G=None,G_error_norm=None,
            signed_velocities=None,velocity_errors=None,a0=None,
            quartet_mismatch=None,quartet_mismatch_error=None,
            both_valleys_occupied=None,positive_relative_stiffness=None),
        strongest_finite_window_diagnostic=strongest,
        finite_window_warning='Not a matched momentum or a thermodynamic uncertainty bound.',
        bare_band_comparison_only=out['bare_band'],
        missing_input='Rank-four normalized in-medium static response and compactification, '
                      'with analytic UV contribution, errors and independent dynamic/occupation checks.',
        amplitude_owner='contrib/qubit_detuned_quartet_20261118; no amplitude inferred here')
    (HERE/'interface.json').write_text(json.dumps(contract,indent=2)+'\n')


def main():
    start=time.monotonic()
    with threadpool_limits(1):
        standard=[ordinary_record(p) for p in standard_paths()]
        standard=[r for r in standard if r]
        paths=sorted((ROOT/'codex/work_C193/ckpt_meas/typhon_out/out').glob('p153c_*_meas.json'))
        long=[long_record(p) for p in paths]
        out=dict(bare_band=bare_band(),standard_records=standard,long_records=long,
                 runtime_seconds=time.monotonic()-start)
    (HERE/'record_results.json').write_text(json.dumps(out,indent=2)+'\n')
    export_interface(out)
    print('Bare:',out['bare_band'])
    for r in standard:
        print('window',r['L'],r['P'],r['chi'],'n',r['n_mid'],'dSchi',r['chi_entropy_change'],
              'Q+pi',[round(x['pair_shift_from_minus_pi'],7) for x in r['carrier_windows']])
    for r in long:
        estimates=[x['shift_from_minus_pi'] for x in r['pair_rows']]
        print('long',r['L'],r['P'],'pair Q+pi span', [min(estimates),max(estimates)] if estimates else None,
              'dressing',r['mean_dressing_representative'])
    print('seconds',out['runtime_seconds'])


if __name__=='__main__':main()
