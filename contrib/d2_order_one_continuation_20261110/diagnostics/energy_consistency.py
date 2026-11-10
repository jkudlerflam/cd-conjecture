#!/usr/bin/env python3
"""Post-hoc open-chain ground-energy consistency, not a dynamical exponent."""
import json
import math
import analyze_saved as A
import numpy as np

LENGTHS=[24,48,96,192]

def analyze(energies):
    E=dict(zip(LENGTHS,energies))
    d={L:E[4*L]-3*E[2*L]+2*E[L] for L in [24,48]}
    ratio=d[48]/d[24] if d[24]!=0 else None
    exponent=-math.log2(ratio) if ratio is not None and ratio>0 else None
    X=np.array([[L,1,1/L] for L in LENGTHS[:3]],dtype=float)
    ebulk,ebound,a=np.linalg.solve(X,energies[:3])
    predicted=float(192*ebulk+ebound+a/192)
    return dict(energies={str(L):e for L,e in E.items()},D={str(L):x for L,x in d.items()},
        ratio=ratio,p_eff=exponent,fit_fixed_inverse_size=dict(e_bulk=float(ebulk),
        e_boundary=float(ebound),inverse_size_coefficient=float(a),predicted_E192=predicted,
        actual_E192=E[192],prediction_residual=E[192]-predicted),
        interpretation='Post-hoc energy correction exponent; not a dynamical exponent or velocity.')

def main():
    xx=[float(sum(-2*np.cos(np.pi*np.arange(1,L//2+1)/(L+1)))) for L in LENGTHS]
    out=dict(exact_xx=analyze(xx),cases=[])
    for file,name in [('bridge_diagnostics.json','b_r5_2_d1_10_q0_hm7_8'),('attractive_diagnostics.json','attr_ref')]:
        data=json.loads((A.HERE/file).read_text());rows=[r for r in data['rows'] if r['name']==name]
        if not all(any(r['L']==L and r['parity']==P for r in rows) for L in LENGTHS for P in [0,1]):
            print('WAIT both parities at all lengths:',name,flush=True);continue
        ground=[min([r for r in rows if r['L']==L],key=lambda r:r['E']) for L in LENGTHS]
        result=analyze([r['E'] for r in ground]);result['name']=name
        result['inputs']=[{k:r[k] for k in ['L','parity','chi','E','source_record','record_sha256','checkpoint_sha256']} for r in ground]
        out['cases'].append(result);print(name,result['D'],result['p_eff'],result['fit_fixed_inverse_size'],flush=True)
    A.D.atomic_json(A.HERE/'energy_consistency.json',out)
    print('XX',out['exact_xx']['p_eff'],flush=True)

if __name__=='__main__':main()
