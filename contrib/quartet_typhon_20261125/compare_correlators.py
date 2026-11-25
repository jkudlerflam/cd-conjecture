"""Stage-to-stage raw and small-quadrature changes, with absolute scales."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from common import write_json
from audit_saved import historical,raw_matrix,filter_matrix,eigenvalues


def metric(a,b):
    return {'max_abs_change':float(np.max(abs(b-a))),
            'relative_l2_change':float(np.linalg.norm(b-a)/max(np.linalg.norm(b),1e-30)),
            'new_max_abs':float(np.max(abs(b))),'new_min_abs':float(np.min(abs(b)))}


def main(old,new,output):
    a=json.loads(gzip.decompress(Path(old).read_bytes()));b=json.loads(gzip.decompress(Path(new).read_bytes()))
    assert all(a[k]==b[k] for k in ('L','h','parity','refs'))
    A=historical();kappas=[A.fit_kappa(x['P']['2'],x['refs'][0]) for x in (a,b)]
    raw=[];filtered=[];fixed_momentum=[]
    for kind in ('P','Q'):
        for key in a[kind]:
            for ref in a['refs']:
                aa=a[kind][key][str(ref)];bb=b[kind][key][str(ref)]
                r=np.arange(len(aa['md']))+aa['j0']-ref
                for lo,hi in ((12,32),(32,56)):
                    s=(r>=lo)&(r<=hi)
                    if not s.any():continue
                    for lab in ('md','dm','mm','dd','mz'):
                        va=np.array(aa[lab]);vb=np.array(bb[lab]);va=va[:,0]+1j*va[:,1];vb=vb[:,0]+1j*vb[:,1]
                        raw.append({'kind':kind,'operator':key,'ref':ref,'window':[lo,hi],'ordering':lab,**metric(va[s],vb[s])})
                for width in (24,32,48):
                    vals=[];rc=None
                    for data,kappa in zip((a,b),kappas):
                        rr,C=raw_matrix(data[kind][key],ref);rc,K,_=filter_matrix(rr,kappa,width,2)
                        if not len(K):break
                        F=np.einsum('ar,rij->aij',K,C);large,small=eigenvalues(F)
                        vals.append(large if kind=='P' else small)
                    if len(vals)<2:continue
                    s=(rc>=24)&(rc<=44)
                    if not s.any():continue
                    filtered.append({'kind':kind,'operator':key,'ref':ref,'width':width,'window':[24,44],
                                     **metric(vals[0][s],vals[1][s])})
                    for stage,kappa in zip(('old','new'),kappas):
                        fixed=[]
                        for data in (a,b):
                            rr,C=raw_matrix(data[kind][key],ref);rc,K,_=filter_matrix(rr,kappa,width,2)
                            F=np.einsum('ar,rij->aij',K,C);large,small=eigenvalues(F)
                            fixed.append(large if kind=='P' else small)
                        fixed_momentum.append({'kind':kind,'operator':key,'ref':ref,'width':width,
                            'window':[24,44],'momentum_from_stage':stage,'kappa_pi':float(kappa/np.pi),
                            **metric(fixed[0][s],fixed[1][s])})
    write_json(output,{'old':str(old),'new':str(new),'old_chi':a['chi'],'new_chi':b['chi'],
        'kappa_pi':[float(k/np.pi) for k in kappas],'raw_changes':raw,'filtered_changes':filtered,
        'filtered_fixed_momentum_changes':fixed_momentum,
        'note':'Finite-L convergence diagnostics. filtered_changes uses each stage kappa; the fixed-momentum variants hold it common. Relative L2 changes use the new-state norm in the denominator. No independent-sample errors or emergent-charge identification.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('old');p.add_argument('new');p.add_argument('--output',required=True);a=p.parse_args()
    main(a.old,a.new,a.output)
