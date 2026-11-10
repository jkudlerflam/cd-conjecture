#!/usr/bin/env python3
"""Incremental read-only diagnostics for completed continuation states."""
import argparse
import hashlib
import json
import math
import time
import analyze_saved as A
import numpy as np
from window_order import window_square,check_tiny
Z=A.Z

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def annotate_snapshots(output):
    snapshots=[]
    for filename in ['state_preservation.json','final_state_preservation.json']:
        manifest=A.HERE.parent/'bridge_scan'/filename
        if manifest.exists():
            snapshots.extend(json.loads(manifest.read_text())['preserved'])
    largest=A.HERE.parent/'bridge_scan/largest_size_preservation.json'
    if largest.exists():
        m=json.loads(largest.read_text())
        snapshots.append(dict(record=m['raw_record'],record_sha256=m['raw_record_sha256'],
            checkpoint=m['preserved_checkpoint'],checkpoint_sha256=m['checkpoint_sha256']))
    for row in output['rows']+output.get('superseded_rows',[]):
        for snap in snapshots:
            if row['record_sha256']==snap['record_sha256'] and row['checkpoint_sha256']==snap['checkpoint_sha256']:
                rec=A.Path(snap['record']);ck=A.Path(snap['checkpoint'])
                if sha(rec)==row['record_sha256'] and sha(ck)==row['checkpoint_sha256']:
                    row['source_record_snapshot']=str(rec.relative_to(A.ROOT))
                    row['checkpoint_snapshot']=str(ck)
                    row['checkpoint_replayable']=True
        if row['name']=='b_r5_2_d1_10_q0_hm7_8' and row['L']==96 and row['parity']==0 and row['chi']==192:
            row['checkpoint_replayable']=False
            row['checkpoint_note']='Original chi192 NPZ overwritten before preservation; stable pre-upgrade diagnostics and raw JSON remain.'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--family',choices=['attractive','bridge'],default='attractive')
    ap.add_argument('--preservation-manifest',help='Measure one immutable raw-record/checkpoint pair from its preservation manifest.')
    args=ap.parse_args()
    source=A.HERE.parent/'attractive_scan/results';dest=A.HERE/(args.family+'_diagnostics.json')
    output=json.loads(dest.read_text()) if dest.exists() else dict(tiny_validation=check_tiny(),rows=[])
    if args.family=='attractive':
        paths=list(source.glob('pref_*_L*.json'))+list(source.glob('attr_ref_L*.json'))
        paths+=list((A.OLD/'dmrg').glob('reference_peotta_exploratory_L*.json'))
    else:
        paths=list((A.HERE.parent/'bridge_scan/dmrg').glob('b_r5_2_d1_10_q0_hm7_8_L*.json'))
    override={}
    if args.preservation_manifest:
        m=json.loads(A.Path(args.preservation_manifest).read_text())
        path=A.Path(m['raw_record']);ck=A.Path(m['preserved_checkpoint'])
        assert sha(path)==m['raw_record_sha256'] and sha(ck)==m['checkpoint_sha256']
        paths=[path];override[str(path)]=ck
    for path in sorted(paths):
        raw=json.loads(path.read_text())
        if raw['status']!='done':continue
        rhash=sha(path);old=next((r for r in output['rows'] if r['record_sha256']==rhash),None)
        if old is not None:continue
        ck=override.get(str(path),A.Path(raw['checkpoint']));chash=sha(ck)
        with np.load(ck,allow_pickle=False) as z:meta=json.loads(str(z['meta']))
        if meta['couplings']!=raw['couplings'] or len(meta['history'])!=len(raw['history']):
            print('SKIP mismatched checkpoint',path.name,flush=True);continue
        t0=time.time();p=Z.MPS.load(ck);L=p.L;W=L//2
        cuts=A.schmidt_diagnostics(p)
        assert max(abs(np.array([c['S'] for c in cuts])-np.array(raw['Sprof'])))<1e-9
        summaries={field:A.D.entropy_stats([c[field] for c in cuts],L)
                   for field in ['S','S_binary','S_conditional']}
        _,g,_=A.D.build_model(L,raw['couplings'],raw.get('distance',1))
        ops={'X':[gj*Z.SP+gj.conjugate()*Z.SM for gj in g],
             'Y':[-1j*gj*Z.SP+1j*gj.conjugate()*Z.SM for gj in g],
             'Z':[Z.ZZ_]*L}
        order={};qgrid=np.linspace(0,np.pi,2001)
        for axis in ops:
            dist,y=np.array(raw['measurements']['long_connected_correlations'][axis]).T
            q=float(qgrid[np.argmax(np.abs(np.exp(-1j*np.outer(qgrid,dist))@y))])
            for label,k in [('zero',0.),('pi',np.pi),('adaptive',q)]:
                weights=np.exp(1j*k*np.arange(L));weights[:L//4]=0;weights[3*L//4:]=0
                sq,mean=window_square(p,ops[axis],weights)
                order[f'{axis}_{label}']=dict(q=k,square=sq,order_square=sq/W**2,
                    connected_order_square=(sq-abs(mean)**2)/W**2,structure_factor=sq/W,
                    mean_real=mean.real,mean_imag=mean.imag)
        group_name='attr_ref' if raw['name']=='reference_peotta_exploratory' else raw['name']
        for prior in output['rows']:
            if prior['name']==group_name:
                assert prior['couplings'].keys()==raw['couplings'].keys()
                assert max(abs(prior['couplings'][k]-raw['couplings'][k]) for k in raw['couplings'])<1e-12
        row=dict(name=group_name,original_name=raw['name'],L=L,parity=raw['parity'],seed=raw['seed'],chi=raw['chi'],
            E=raw['E_physical'],couplings=raw['couplings'],source_record=str(path.relative_to(A.ROOT)),
            record_sha256=rhash,checkpoint=str(ck),record_checkpoint=raw['checkpoint'],checkpoint_sha256=chash,cuts=cuts,
            summaries=summaries,envelope=A.envelope(raw),density=A.density_diagnostics(raw),
            window=[L//4,3*L//4],order=order,seconds=time.time()-t0)
        if sha(path)!=rhash or sha(ck)!=chash:
            print('SKIP changed input during measurement',path.name,flush=True);continue
        def same_state(prior):
            return all(prior[k]==row[k] for k in ['name','L','parity','seed'])
        for prior in output['rows']:
            if same_state(prior):
                oldpath=A.ROOT/prior['source_record']
                for snapshot in (oldpath.parent.parent/'validation_snapshots').glob(oldpath.stem+'*.json'):
                    if sha(snapshot)==prior['record_sha256']:
                        prior['source_record_snapshot']=str(snapshot.relative_to(A.ROOT))
                output.setdefault('superseded_rows',[]).append(prior)
        output['rows']=[r for r in output['rows'] if not same_state(r)]+[row]
        output['size_slopes']=[]
        for name in sorted(set(r['name'] for r in output['rows'])):
            rs=[r for r in output['rows'] if r['name']==name]
            for selection in ['ground_variational','even','odd']:
                pick=[]
                for length in sorted(set(r['L'] for r in rs)):
                    candidates=[r for r in rs if r['L']==length and
                        (selection=='ground_variational' or r['parity']==int(selection=='odd'))]
                    if candidates:pick.append(min(candidates,key=lambda r:r['E']))
                for a,b in zip(pick,pick[1:]):
                    entry=dict(name=name,selection=selection,sizes=[a['L'],b['L']],
                        parities=[a['parity'],b['parity']],slopes={})
                    for field in ['S','S_binary','S_conditional']:
                        entry['slopes'][field]={key:6*(b['summaries'][field][key]-a['summaries'][field][key])/math.log(b['L']/a['L'])
                            for key in ['S_central3','S_fractional']}
                    output['size_slopes'].append(entry)
        annotate_snapshots(output)
        A.D.atomic_json(dest,output)
        print(path.name,'X adaptive',order['X_adaptive'],'S',summaries['S']['S_fractional'],
              'binary',summaries['S_binary']['S_fractional'],flush=True)
    annotate_snapshots(output)
    A.D.atomic_json(dest,output)

if __name__=='__main__':main()
