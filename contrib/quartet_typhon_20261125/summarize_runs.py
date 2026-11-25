"""Recompute matched central-quarter entropies from retained full profiles."""
import argparse
import json
from pathlib import Path
import numpy as np
from common import HERE,write_json

# Rounded same-chi references, all L128 even parity, open_chain_ladders.txt.
REFERENCES={1.783:{64:1.50606,128:1.67048,256:1.75439,384:1.78305,512:1.79175},
            1.772:{64:.75290,128:.75453,256:.75457,384:.75458}}


def summarize():
    rows=[]
    for f in sorted((HERE/'results').glob('*/results/*_summary.json')):
        m=json.loads(f.read_text());p=f.with_name(f.name.replace('_summary.json','_profiles.npz'))
        if not p.exists():continue
        with np.load(p) as v:S=v['S']
        L=m['L'];mid=np.abs(np.arange(1,L)-L/2)<=L/8
        row={k:m.get(k) for k in ('L','h','parity','chi','E','Smid','converged','shelved','sweeps','truncation','dmrg_seconds','total_seconds','max_rss_mb','measured')}
        row.update(run=f.parents[1].name,Sbar_central_quarter=float(S[mid].mean()),
                   summary_Sbar=m.get('Sbar'),summary_Sbar_window=m.get('Sbar_window','central half; initial wrapper label corrected here'),
                   abs_dE_last=None if m['dE_last'] is None else abs(m['dE_last']),
                   abs_dS_last=None if m['dS_last'] is None else abs(m['dS_last']))
        if L==128 and m['parity']==0:
            ref=REFERENCES.get(m['h'],{}).get(m['chi'])
            row.update(reference_Smid_same_chi=ref,Smid_minus_same_chi_reference=None if ref is None else m['Smid']-ref)
            if m['h']==1.783:
                row['E_minus_chi1024_reference']=m['E']-(-133.564501760)
                row['energy_comparison_note']='Higher-chi reference, not same-chi convergence certificate'
        rows.append(row)
    write_json(HERE/'analysis/run_summary.json',{'rows':rows,'Sbar_definition':'mean over |bond-L/2|<=L/8, matching ana_obc.py',
        'note':'First pilot/control snapshot used the label Sbar for a central-half average. Raw files preserved; this table recomputes the correct central-quarter statistic from S profiles.'})
    for x in rows:print({k:x.get(k) for k in ('run','h','chi','E','Smid','Sbar_central_quarter','converged','measured')})


if __name__=='__main__':summarize()
