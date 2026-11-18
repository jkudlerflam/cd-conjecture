#!/usr/bin/env python3
"""Saved-record accuracy, boundary, shell and momentum diagnostics."""
import time,sys,resource
import numpy as np
from threadpoolctl import threadpool_limits
from common import HERE,ROOT,PRIOR,OLD,read,write,sources,candidate_records


def chord(L):
    l=np.arange(1,L);return np.log(L/np.pi*np.sin(np.pi*l/L))


def summary(row):
    L=row['L'];l=np.arange(1,L);S=np.array(row['S']);x=chord(L)
    den=np.array(row['operators']['density']['profile'] if 'operators' in row else row['density'])
    row=dict(row);row['density']=den.tolist();row['averages']={}
    for name,m in [(f'W{w}',abs(l-L//2)<=w) for w in (0,1,6,12)]+[
            ('quarter',(l>=3*L/8)&(l<=5*L/8)),('half',(l>=L/4)&(l<=3*L/4))]:
        row['averages'][name]=dict(S=float(S[m].mean()),x=float(x[m].mean()))
    row['n_half']=float(den[L//4:3*L//4].mean());row['n_quarter']=float(den[3*L//8:5*L//8].mean())
    row['N_expectation']=float(den.sum());return row


def pair(a,b):
    return dict(tags=[a['tag'],b['tag']],lengths=[a['L'],b['L']],parities=[a.get('P'),b.get('P')],
        density_half=[a['n_half'],b['n_half']],delta_density=b['n_half']-a['n_half'],
        estimates={k:dict(logL=float(6*(b['averages'][k]['S']-a['averages'][k]['S'])/np.log(b['L']/a['L'])),
            chord=float(6*(b['averages'][k]['S']-a['averages'][k]['S'])/(b['averages'][k]['x']-a['averages'][k]['x'])))
            for k in a['averages']})


def peak(y):
    y=np.array(y);M=len(y);q=np.linspace(2*np.pi/M,np.pi,4097)
    z=np.exp(-1j*np.outer(q,np.arange(M)))@y;k=np.argmax(abs(z))
    return dict(q_over_pi=float(q[k]/np.pi),amplitude=float(2*abs(z[k])/M),resolution_over_pi=2/M)


def regress(S,profiles,L,edge):
    l=np.arange(1,L);m=(l>=edge*L)&(l<=(1-edge)*L);cuts=l[m]
    x=chord(L)[m];X=np.column_stack([np.ones(len(x)),x]);y=np.array(S)[m]
    beta=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@beta
    out=dict(edge_fraction=edge,c_profile=float(6*beta[1]),residual_RMS=float(np.sqrt(np.mean(res**2))),
        entropy_residual_peak=peak(res),cuts=cuts.tolist(),entropy_residual=res.tolist(),channels={})
    for name,p in profiles.items():
        v=np.array(p)[m];vb=np.linalg.lstsq(X,v,rcond=None)[0];vr=v-X@vb
        rho=float(res@vr/np.sqrt((res@res)*(vr@vr))) if vr@vr>1e-24 else 0.
        full=np.column_stack([X,v]);coef=np.linalg.lstsq(full,y,rcond=None)[0];rr=y-full@coef
        cv=[]
        for side in (False,True):
            train=cuts<L/2 if side else cuts>L/2;test=cuts>L/2 if side else cuts<L/2
            b0=np.linalg.lstsq(X[train],y[train],rcond=None)[0];b1=np.linalg.lstsq(full[train],y[train],rcond=None)[0]
            e0=float(np.sqrt(np.mean((y[test]-X[test]@b0)**2)));e1=float(np.sqrt(np.mean((y[test]-full[test]@b1)**2)))
            cv.append(dict(train_side='left' if side else 'right',baseline_RMS=e0,augmented_RMS=e1))
        out['channels'][name]=dict(partial_correlation=rho,residual_R2=rho*rho,nuisance_coefficient=float(coef[2]),
            c_with_profile=float(6*coef[1]),augmented_residual_RMS=float(np.sqrt(np.mean(rr**2))),
            profile_residual_peak=peak(vr),half_holdout=cv,
            scope='Descriptive same-state regression; reflected halves are not independent samples or a validated entropy subtraction.')
    return out


def controls():
    rows=[]
    for d in read(PRIOR/'controls.json')['rows']:
        if 'S' not in d:continue
        L=d['L'];legs=int(d['kind'].split('_')[-2]);m=L//legs;N=int(np.floor(.296*(m+1)))
        u=np.sqrt(2/(m+1))*np.sin(np.pi*np.arange(1,m+1)[:,None]*np.arange(1,N+1)/(m+1));C=u@u.T
        den=np.empty(L)
        for leg in range(legs):den[leg::legs]=np.diag(C)
        crossing=[]
        for cut in range(1,L):
            e=0.
            for leg in range(legs):
                nleft=len(range(leg,cut,legs))
                if 0<nleft<m:e+=-2*C[nleft-1,nleft]
            crossing.append(e)
        row=summary(dict(tag=d['kind']+f'_L{L}',kind=d['kind'],L=L,S=d['S'],density=den.tolist()))
        ps=dict(density_cut=((den[:-1]+den[1:])/2).tolist(),crossing_kinetic=crossing)
        row['residual_tests']=[regress(d['S'],ps,L,e) for e in (1/8,1/4)]
        rows.append(row)
    return rows


def momentum(d):
    m=d['measurement'];M=m['gwin'][1];q=np.linspace(-np.pi,np.pi,8192,endpoint=False)
    phase=np.exp(-1j*np.outer(q,np.arange(M)));out={}
    for label in ('b','f'):
        G=np.array(m['G'+label+'_re'])+1j*np.array(m['G'+label+'_im'])
        nk=np.real(np.einsum('qi,ij,qj->q',phase.conj(),G,phase,optimize=True))/M
        assert min(nk)>-1e-8
        peaks=[];weights=[]
        for mask in (np.cos(q)>=0,np.cos(q)<0):
            ids=np.flatnonzero(mask);k=ids[np.argmax(nk[mask])]
            peaks.append(dict(k=float(q[k]),k_over_pi=float(q[k]/np.pi),height=float(nk[k])))
            weights.append(float(nk[mask].sum()/nk.sum()))
        Q=abs(np.angle(np.exp(1j*(peaks[0]['k']-peaks[1]['k']))))
        K=np.angle(np.exp(1j*(peaks[0]['k']+peaks[1]['k'])))
        # Symmetry-related two peaks are not two measured gapless modes.
        out[label]=dict(peaks=peaks,partition_weights=weights,peak_difference_over_pi=float(Q/np.pi),
            peak_sum_over_pi=float(K/np.pi),q_grid=q[::32].tolist(),n_q=nk[::32].tolist(),
            mean_occupation_integral=float(nk.mean()),density=float(np.trace(G).real/M))
        assert abs(nk.mean()-np.trace(G).real/M)<1e-10
    return dict(tag=d['tag'],L=d['L'],P=d['P'],chi=d['chi'],window=m['gwin'],resolution_over_pi=2/M,
        definitions='Gij=<a_i^dagger a_j>, n(k)=sum exp[ik(i-j)] Gij/M; physical pre-common-rotation frame; no peak count implies c.',channels=out)


def main():
    start=time.monotonic();before=sources();rows=[summary(r) for r in candidate_records()]
    residuals=[]
    for r in rows:
        den=np.array(r['density']);profiles=dict(density_cut=(den[:-1]+den[1:])/2)
        if 'operators' in r:
            profiles.update(bond=r['operators']['bond']['profile'],kinetic=r['operators']['kinetic']['profile'])
        residuals.append(dict(tag=r['tag'],state_status=r['state_status'],tests=[regress(r['S'],profiles,r['L'],e) for e in (1/8,1/4)]))
    complete=[r for r in rows if r['state_status']=='completed'];table={(r['L'],r['P']):r for r in complete}
    pairs=[]
    for L,L2 in ((48,96),(96,144)):
        for P in (0,1):pairs.append(dict(sequence=f'parity_{P}',**pair(table[L,P],table[L2,P])))
        a=min((table[L,P] for P in (0,1)),key=lambda r:r['E']);b=min((table[L2,P] for P in (0,1)),key=lambda r:r['E'])
        pairs.append(dict(sequence='ground_completed',**pair(a,b)))
    interrupted=next(r for r in rows if r['state_status']=='interrupted')
    pairs += [dict(sequence='odd_interrupted',**pair(table[96,1],interrupted)),
              dict(sequence='ground_interrupted',**pair(table[96,0],interrupted))]
    den=np.log(1.5);switch=6*(table[96,1]['averages']['W6']['S']-table[96,0]['averages']['W6']['S'])/den
    base=pair(table[96,1],table[144,1])['estimates']['W6']['logL']
    accuracy=6*(interrupted['averages']['W6']['S']-table[144,1]['averages']['W6']['S'])/den
    final=pair(table[96,0],interrupted)['estimates']['W6']['logL'];assert abs(final-base-switch-accuracy)<1e-12
    decomposition=dict(fixed_odd_completed=base,change_from_96_parity_selection=switch,
        change_from_144_chi_increase=accuracy,ground_interrupted=final,
        residual_above_two_in_fixed_odd=base-2,
        scope='Exact finite-record arithmetic, not a causal separation or license to replace the ground state by an excited parity.')
    sector=[]
    for L in (48,96,144):
        a,b=table[L,0],table[L,1];sector.append(dict(L=L,E_odd_minus_even=b['E']-a['E'],
            L_times_abs_parity_gap=L*abs(b['E']-a['E']),density_odd_minus_even=b['n_half']-a['n_half'],
            Sbar6_odd_minus_even=b['averages']['W6']['S']-a['averages']['W6']['S']))
    stages=[]
    paths=list((PRIOR/'outputs').glob('*.json'))+[OLD/'outputs'/f'legacy_V4_a1_L96_P{P}.json' for P in (0,1)]
    for p in paths:
        d=read(p);hs=[r for r in d['hist'] if 'Sbar6' in r and r['alpha']==0]
        if not hs:continue
        stages.append(dict(tag=p.stem,L=d['cfg']['L'],P=d['cfg']['P'],status=d['status'],
            sweeps=[{k:r[k] for k in ('stage','chi','E_phys','Sbar1','Sbar6','Sbar12','n_mid','t','trunc','res_max')} for r in hs]))
    con=controls();cpairs=[]
    for kind in sorted({r['kind'] for r in con}):
        group=sorted([r for r in con if r['kind']==kind],key=lambda r:r['L'])
        cpairs += [dict(kind=kind,**pair(a,b)) for a,b in zip(group[:-1],group[1:])]
    momenta=[momentum(read(p)) for p in sorted((HERE/'onebody').glob('*.json'))]
    kinematics=[]
    for r in rows:
        if 'operators' not in r:continue
        q=r['operators']['density']['half']['profile_peak']['q_over_pi'];n=r['n_half']
        fold=lambda x:float(abs(np.angle(np.exp(1j*np.pi*x)))/np.pi)
        predictions=dict(one_charge1_first=fold(2*n),one_charge1_second=fold(4*n),
            pair_liquid_first=fold(n),pair_liquid_third=fold(3*n),balanced_two_intra=fold(n),
            bare_two_valley_coherence=1-2*np.arcsin(.25)/np.pi)
        md=next((a for a in momenta if a['tag']==r['tag']),None)
        if md:predictions['measured_boson_coherence']=md['channels']['b']['peak_difference_over_pi']
        kinematics.append(dict(tag=r['tag'],L=r['L'],P=r['P'],raw_density=n,density_profile_q_over_pi=q,
            density_connected_q_over_pi=r['operators']['density']['half']['connected_peak']['q_over_pi'],
            half_window_resolution_over_pi=4/r['L'],predictions=predictions,
            caveat='Density formulas assume conserved dressed charge; raw microscopic N is not conserved. Higher-harmonic coincidences are not mode counts.'))
    compact=[{k:v for k,v in r.items() if k not in ('S','operators','density')} for r in rows]
    write(HERE/'audit.json',dict(states=compact,pairs=pairs,decomposition=decomposition,sectors=sector,
        residual_tests=residuals,accuracy_stages=stages,controls=con,control_pairs=cpairs,
        momenta=momenta,kinematics=kinematics,asymptotic_c=None,
        runtime_seconds=time.monotonic()-start,peak_RSS_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9))
    after=sources();assert before==after
    write(HERE/'source_manifest.json',dict(source_hashes=before,archive_unchanged=True))
    print('decomposition',decomposition)
    for r in residuals:
        t=r['tests'][1];print('residual',r['tag'],t['c_profile'],{k:round(v['partial_correlation'],6) for k,v in t['channels'].items()})
    for m in momenta:print('momenta',m['tag'],m['channels']['b']['peaks'],m['channels']['b']['peak_difference_over_pi'])
    print('seconds',time.monotonic()-start)


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
