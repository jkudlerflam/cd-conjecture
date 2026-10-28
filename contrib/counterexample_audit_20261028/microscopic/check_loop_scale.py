"""Explain the finite-loop momentum shifts using saved certified matrices.
This is a numerical scale diagnostic, not a replacement for root certification.
"""
from pathlib import Path
from fractions import Fraction
import json
import math
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from scipy.linalg import eig

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SRC=ROOT/'codex/work_C33'

def main():
    raw=np.load(SRC/'full_loop_n23.npz',allow_pickle=True)
    white=np.load(SRC/'full_loop_whitened_n23.npz',allow_pickle=True)
    derivative=np.load(SRC/'derivative_n23.npz',allow_pickle=True)
    roots=json.loads((SRC/'full_loop_bands_n23.json').read_text())['roots']
    exact=json.loads((SRC/'root_certificate_exact_n23.json').read_text())
    bits=int(raw['bits']);den=2*(1<<(2*bits));n=len(raw['F'])
    Tf=-np.array([[int(v)/den for v in row] for row in raw['F']])
    Mf=1j*np.array([[int(v)/den for v in row] for row in raw['M']])
    T=white['T'];M=white['mu'];field=white['field']
    W=field/2.**(-214)
    I=np.eye(n);Z=np.zeros_like(I)
    def locate(onsite):
        zs=eig(np.block([[-onsite,-T.T],[I,Z]]),np.block([[T,Z],[Z,I]]),right=False)
        return np.sort(np.angle(zs))
    reference=locate(field)
    finite=locate(field+M)
    assert np.max(abs(finite-roots))<1e-12
    centers=(finite[:n]+finite[n:][::-1])/2
    reference_centers=(reference[:n]+reference[n:][::-1])/2
    rows=[]
    for idx,k in enumerate(reference):
        H=field+T*np.exp(1j*k)+T.T*np.exp(-1j*k)
        vals,U=np.linalg.eigh(H);u=U[:,np.argmin(abs(vals))]
        derivativeH=1j*T*np.exp(1j*k)-1j*T.T*np.exp(-1j*k)
        derivative_whitened=float(np.vdot(u,derivativeH@u).real)
        metric=float(np.vdot(u,W@u).real)
        delta_whitened=float(np.vdot(u,M@u).real)
        linear_shift=-delta_whitened/derivative_whitened
        rows.append(dict(index=idx,reference_root=float(k),finite_loop_root=float(finite[idx]),
                         exact_numeric_shift=float(finite[idx]-k),linear_shift=linear_shift,
                         metric=metric,normalized_velocity=abs(derivative_whitened)/metric,
                         onsite_expectation_physical_normalized=delta_whitened/metric,
                         onsite_expectation_whitened=delta_whitened,
                         derivative_whitened=derivative_whitened))
    intervals=[]
    for pair in exact['intervals']:
        intervals.append([2*math.atan(float(Fraction(int(a),int(b)))) for a,b in pair])
    center_bounds=[]
    for i in range(n):
        j=2*n-1-i
        center_bounds.append([(intervals[i][0]+intervals[j][0])/2,
                              (intervals[i][1]+intervals[j][1])/2])
    order=sorted(range(n),key=lambda i:centers[i])
    minseparation=min(center_bounds[j][0]-center_bounds[i][1] for i,j in zip(order,order[1:]))
    # Inspect the tree filter metric itself in its ordinary coordinates.
    tree_scale=1<<(2*int(derivative['bits']))
    Ftree=np.array([[int(v)/tree_scale for v in row] for row in derivative['F_numerator']])
    data=dict(status='scale_diagnostic_using_saved_certified_entries',
              beta=1/int(raw['detC']),loop_ratio=2.**(-136),
              log10_td=-math.log10(int(raw['detC']))-136*math.log10(2),
              normalized_T_norm=float(np.linalg.norm(Tf,2)),
              normalized_M_norm=float(np.linalg.norm(Mf,2)),
              normalized_M_to_T_norm=float(np.linalg.norm(Mf,2)/np.linalg.norm(Tf,2)),
              normalized_M_per_loop_ratio=float(np.linalg.norm(Mf,2)/2.**(-136)),
              whitened_T_norm=float(np.linalg.norm(T,2)),
              whitened_M_norm=float(np.linalg.norm(M,2)),
              whitened_M_per_loop_ratio=float(np.linalg.norm(M,2)/2.**(-136)),
              min_velocity_normalized=min(r['normalized_velocity'] for r in rows),
              max_velocity_normalized=max(r['normalized_velocity'] for r in rows),
              max_reference_center_abs=float(np.max(abs(reference_centers))),
              root_centers=centers.tolist(),certified_bracket_center_bounds=center_bounds,
              bracket_min_center_separation=minseparation,
              min_center=float(min(centers)),max_center=float(max(centers)),
              max_linear_shift_error=max(abs(r['exact_numeric_shift']-r['linear_shift']) for r in rows),
              root_details=rows)
    (HERE/'loop_scale_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['root_details','root_centers','certified_bracket_center_bounds']},indent=2))
    print(json.dumps(dict(witness_centers=[dict(index=i,center=centers[i],bounds=center_bounds[i]) for i in [7,8]])))

if __name__=='__main__':main()
