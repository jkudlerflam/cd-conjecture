"""Enlarged-matrix diagnostics without assigning an emergent charge.

At finite separation C is not a positive Gram matrix. Singular values and
fixed probe subtractions below are sensitivity diagnostics, not a spectral
projection onto primary fields. All raw normal/anomalous entries are kept.
"""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from common import write_json
from measure_matrix import hermitian_basis
from audit_saved import historical,filter_matrix,eigenvalues,signed_fit,forward_fit


def main(matrix_file,quartet_file,output):
    m=json.loads(gzip.decompress(Path(matrix_file).read_bytes()))
    q=json.loads(gzip.decompress(Path(quartet_file).read_bytes()))
    A=historical();kappa=A.fit_kappa(q['P']['2'],q['refs'][0])
    names,T=hermitian_basis(centered=m.get('centered',False));rows=[]
    for ref,raw in m['rows'].items():
        r=np.array(raw['r'],float);C=np.array(raw['C']);C=C[...,0]+1j*C[...,1]
        H=np.einsum('ai,rij,bj->rab',T,C,T)
        for width in (24,32,40,48):
            rc,K,cond=filter_matrix(r,kappa,width,2)
            if not len(K):continue
            F=np.einsum('ar,rij->aij',K,H)
            sym=(F.real+np.swapaxes(F.real,-1,-2))/2
            region=(rc>=24)&(rc<=44)
            scales=np.sqrt(np.maximum(np.mean(abs(np.diagonal(sym,axis1=1,axis2=2)),axis=0),1e-30))
            normalized=sym/scales[None,:,None]/scales[None,None,:]
            sv=np.linalg.svd(normalized,compute_uv=False)
            item={'ref':int(ref),'width':width,'r':rc.tolist(),'condition':cond,
                  'imaginary_max':float(np.max(abs(F.imag))),
                  'antisymmetric_max':float(np.max(abs(F.real-np.swapaxes(F.real,-1,-2)))),
                  'constant_scales':scales.tolist(),'normalized_singular_values':sv.tolist(),
                  'C_hermitian_real':sym.tolist(),'subtractions':[]}
            for probe_set in ((4,),(7,),(4,7),(4,5,6,7)):
                # One constant linear combination across the fitted window,
                # never a different optimized eigenvector at each distance.
                for lo,hi in ((24,40),(28,44)):
                    use=(rc>=lo)&(rc<=hi)
                    if use.sum()<3:continue
                    ps=list(probe_set)
                    N=sym[:,ps,:][:,:,ps]
                    weights=[];residuals=[]
                    for a in range(4):
                        y=sym[:,a,ps]
                        X=N[use].reshape(-1,len(ps));v=y[use].reshape(-1)
                        coef=np.linalg.lstsq(X,v,rcond=1e-10)[0]
                        weights.append(coef)
                        residuals.append(float(np.linalg.norm(X@coef-v)/max(np.linalg.norm(v),1e-30)))
                    V=np.eye(8)
                    for a in range(4):V[a,ps]-=weights[a]
                    transformed=np.einsum('ai,rij,bj->rab',V,sym,V)
                    signals={}
                    for shape,inds in (('2_5',[0,1]),('3_6',[2,3])):
                        _,before=eigenvalues(sym[:,inds,:][:,:,inds])
                        _,after=eigenvalues(transformed[:,inds,:][:,:,inds])
                        signals[shape]={'before':before.tolist(),'after':after.tolist(),
                            'change_relative_l2':float(np.linalg.norm(after-before)/max(np.linalg.norm(before),1e-30)),
                            'before_signed_fit':signed_fit(rc,before,24,44),
                            'after_signed_fit':signed_fit(rc,after,24,44)}
                    item['subtractions'].append({'probes':[names[p] for p in ps],'fit_window':[lo,hi],
                             'weights':np.array(weights).tolist(),'relative_cross_residuals':residuals,'signals':signals})
            rows.append(item)
    result={'matrix_file':str(matrix_file),'quartet_file':str(quartet_file),'L':m['L'],'chi':m['chi'],
            'kappa_pi':float(kappa/np.pi),'basis':names,'rows':rows,
            'interpretation':'Finite-r sensitivity diagnostic only. No charge assignment or primary scaling estimate.'}
    write_json(output,result)
    print('matrix diagnostics',m['L'],m['chi'],'rows',len(rows))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('matrix');p.add_argument('quartet');p.add_argument('--output',required=True);a=p.parse_args()
    main(a.matrix,a.quartet,a.output)
