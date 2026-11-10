"""Read-only provenance and convergence audit; partial runs remain explicit."""
import argparse,hashlib,json,pickle
from collections import Counter
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--final',action='store_true');a=ap.parse_args()
errors=[];rows=[];counts=Counter()
for f in sorted((HERE/'results').glob('*.json')):
    r=json.loads(f.read_text());counts[r['status']]+=1
    if a.final and r['status']=='running':errors.append(f'{f.name}: still running')
    for source,sha in r['source_hashes'].items():
        snapshots=list((HERE/'source_snapshots').glob(sha+'_*'))
        if not snapshots or any(hashlib.sha256(s.read_bytes()).hexdigest()!=sha for s in snapshots):errors.append(f'{f.name}: missing/bad source {sha}')
    row=dict(source=f.name,status=r['status'],converged=r.get('converged'))
    if r['status']=='done':
        ck=Path(r['checkpoint'])
        if not ck.exists() or hashlib.sha256(ck.read_bytes()).hexdigest()!=r['checkpoint_sha256']:
            errors.append(f'{f.name}: missing or mismatched checkpoint');continue
        with ck.open('rb') as h:saved=pickle.load(h)
        for key in ['L','geometry','parity','couplings']:
            if saved[key]!=r[key]:errors.append(f'{f.name}: wrong checkpoint {key}')
        if saved.get('gauge',False)!=r.get('gauge',False):errors.append(f'{f.name}: wrong checkpoint gauge')
        if saved['checkpoint_kind']!='final':errors.append(f'{f.name}: not a final checkpoint')
        psi=saved['psi']
        err=max(float(np.max(np.abs(psi.entanglement_entropy(n=n)-r['entropy'][str(n)]))) for n in [1,2,3])
        if err>1e-9:errors.append(f'{f.name}: checkpoint entropy differs by {err}')
        if r['converged']!=all(r['convergence_checks'].values()):errors.append(f'{f.name}: convergence flags inconsistent')
        inf=r['geometry']=='infinite';st=r['sweep_stats'];cc=r['convergence_checks']
        expected=dict(recent_energy_stable=len(st['E'][-4:])==4 and np.ptp(st['E'][-4:])<(1e-8 if inf else 1e-7),
            recent_entropy_stable=len(st['S'][-4:])==4 and np.ptp(st['S'][-4:])<1e-5,
            local_residual_small=r['local_solver'][-1]['max_absolute_residual']<1e-8,
            norm_before_small=r['norm_before']<1e-5,
            final_sweep_norm_small=st['norm_err'][-1]<1e-5,
            energy_agreement=abs(r['E_returned']-r['E_direct_shifted'])<(1e-7 if inf else 1e-8))
        for key,value in cc.items():
            if expected[key]!=value:errors.append(f'{f.name}: wrongly evaluated gate {key}')
        shift=r['couplings']['Jz']+r['couplings']['hz'] if inf else (r['L']-1)*r['couplings']['Jz']+r['L']*r['couplings']['hz']
        if abs(r['E_physical']-r['E_direct_shifted']-shift)>1e-11:errors.append(f'{f.name}: wrong physical energy shift')
        row['entropy_reproduction_error']=err
        row['checkpoint_verified']=True
        row['all_current_gates_pass']=bool(all(expected.values()))
    rows.append(row)
warm=[]
for f in sorted(HERE.glob('warm_snapshot*.json')):
    data=json.loads(f.read_text());data=data if isinstance(data,list) else [data]
    for r in data:
        target=ROOT/r['snapshot'];ok=target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()==r['sha256']
        if not ok:errors.append(f'{f.name}: warm snapshot missing or changed')
        warm.append(dict(snapshot=r['snapshot'],verified=ok))
spectra=[]
for f in sorted(HERE.glob('dispersion_L*.json')):
    d=json.loads(f.read_text());maximum=max(x['residual'] for x in d['sectors'])
    if a.final and not d['complete']:errors.append(f'{f.name}: incomplete momentum scan')
    if maximum>=1e-8:errors.append(f'{f.name}: failed Ritz residual')
    if len({(x['parity'],x['k']) for x in d['sectors']})!=len(d['sectors']):errors.append(f'{f.name}: duplicate sector')
    spectra.append(dict(file=f.name,complete=d['complete'],max_residual=maximum))
vumps=[];vcounts=Counter()
for f in sorted((HERE/'vumps_results').glob('*.json')):
    r=json.loads(f.read_text());vcounts[r['status']]+=1
    if a.final and r['status']=='running':errors.append(f'{f.name}: VUMPS still running')
    for name,sha in r['source_hashes'].items():
        snaps=list((HERE/'source_snapshots').glob(sha+'_*'))
        if not snaps or any(hashlib.sha256(s.read_bytes()).hexdigest()!=sha for s in snaps):errors.append(f'{f.name}: bad VUMPS source archive')
    source=ROOT/r['source']
    if hashlib.sha256(source.read_bytes()).hexdigest()!=r['source_sha256']:errors.append(f'{f.name}: changed source record')
    source_record=json.loads(source.read_text())
    if source_record['checkpoint_sha256']!=r['source_checkpoint_sha256']:errors.append(f'{f.name}: wrong source checkpoint hash')
    if source_record['converged']!=r['source_converged']:errors.append(f'{f.name}: wrong source convergence flag')
    for key in ['model','L','chi','gauge','couplings']:
        if source_record.get(key,False)!=r[key]:errors.append(f'{f.name}: changed source {key}')
    row=dict(source=f.name,status=r['status'],converged=r.get('converged'))
    if r['status']=='done':
        ck=Path(r['checkpoint'])
        if hashlib.sha256(ck.read_bytes()).hexdigest()!=r['checkpoint_sha256']:errors.append(f'{f.name}: changed VUMPS checkpoint')
        with ck.open('rb') as h:saved=pickle.load(h)
        for key in ['L','couplings','gauge']:
            if saved[key]!=r[key]:errors.append(f'{f.name}: wrong VUMPS {key}')
        if saved['checkpoint_kind']!='final' or saved['algorithm']!='VUMPS':errors.append(f'{f.name}: wrong VUMPS checkpoint kind')
        err=max(float(np.max(np.abs(saved['psi'].entanglement_entropy(n=n)-r['entropy'][str(n)]))) for n in [1,2,3])
        if err>1e-9:errors.append(f'{f.name}: VUMPS entropy mismatch')
        st=r['sweep_stats'];ee=[complex(*x) if isinstance(x,list) else complex(x) for x in st['E'][-4:]]
        expected=dict(recent_energy_stable=len(ee)==4 and np.ptp(np.real(ee))<1e-8,
            recent_entropy_stable=len(st['S'][-4:])==4 and np.ptp(st['S'][-4:])<1e-5,
            norm_small=r['norm']<1e-7,final_sweep_norm_small=st['norm_err'][-1]<1e-7,
            split_small=st['max_split_err'][-1]<1e-7,tangent_small=r['max_tangent']<1e-7,
            energy_agreement=max(abs(r['E_returned']-r['E_direct_shifted']),max(abs(np.array(r['environment_energies_shifted'])-r['E_direct_shifted'])))<1e-8)
        if expected!=r['convergence_checks'] or r['converged']!=all(expected.values()):errors.append(f'{f.name}: inconsistent VUMPS gates')
        shift=r['couplings']['Jz']+r['couplings']['hz']
        if abs(r['E_physical']-r['E_direct_shifted']-shift)>1e-11:errors.append(f'{f.name}: wrong VUMPS physical energy shift')
        row['entropy_reproduction_error']=err
    vumps.append(row)
if a.final:
    needed=[f'candidate_infinite_L2_P0_s17_chi{c}_g'+('_one' if c<=128 else '')+'_vumps.json' for c in [64,96,128,192,256]]
    needed+=['ising_critical_infinite_L2_P0_s23_chi24_a1_vumps.json']
    for name in needed:
        if not (HERE/'vumps_results'/name).exists():errors.append('Missing planned VUMPS check: '+name)
result=dict(counts=dict(counts),errors=errors,records=rows,warm_snapshots=warm,spectra=spectra,vumps_counts=dict(vcounts),vumps_records=vumps)
(HERE/'audit_report.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(counts=dict(counts),errors=errors)))
raise SystemExit(bool(errors))
