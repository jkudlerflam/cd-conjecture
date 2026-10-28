"""Aggregate the two independent criticality diagnostics without claiming c>2."""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
Ls=np.array([8,12,16])
data={L:json.loads((HERE/f'result_L{L}.json').read_text()) for L in Ls}
responses={L:json.loads((HERE/f'response_L{L}.json').read_text()) for L in [12,16]}
out={'status':'numerical_c2_control_not_a_counterexample','local_dimension':3,'sizes':Ls.tolist(),'runs':[]}
for ix in range(3):
    ds=[data[L]['runs'][ix] for L in Ls]
    J=ds[0]['J']; assert all(d['J']==J for d in ds)
    ss=np.array([d['entropy'][str(L//2)] for L,d in zip(Ls,ds)])
    cfit=float(3*np.polyfit(np.log(Ls),ss,1)[0])
    cpair=float(3*(ss[-1]-ss[-2])/np.log(Ls[-1]/Ls[-2]))
    centers=[float(3*(d['entropy'][str(L//2)]-d['entropy'][str(L//2-2)])/(-np.log(np.cos(2*np.pi/L)))) for L,d in zip(Ls,ds)]
    es=np.array([d['energies'][0] for d in ds])/Ls
    poly=np.polyfit(1/Ls.astype(float)**2,es,2)
    casimir=float(-6/np.pi*poly[1])
    casimir2=float(-6/np.pi*np.polyfit(1/Ls[-2:].astype(float)**2,es[-2:],1)[0])
    response=responses[16]['runs'][ix]
    v_sum=response['velocity_sum']
    max_flux_difference=max(abs(v[0]['vK']-v[1]['vK']) for v in response['flux'].values())
    entry={'J':J,'half_chain_entropy':ss.tolist(),'c_entropy_three_sizes':cfit,
        'c_entropy_L12_L16':cpair,'center_shape_c_by_L':centers,
        'casimir_velocity_sum_Lminus1_Lminus3_fit':casimir,
        'casimir_velocity_sum_two_largest_sizes':casimir2,'fitted_bulk_energy':float(poly[2]),
        'response_velocity_sum_L16':v_sum,'response_charge_velocity_L16':response['charge']['v'],
        'response_spin_velocity_L16':response['spin']['v'],'casimir_to_response_ratio':casimir/v_sum,
        'velocity_sum_first_order':float(2+J*(2-12/np.pi**2)),
        'max_flux_vK_change_phi_halving':float(max_flux_difference),
        'energy_remainder_over_J2_per_L':[float(d['energy_remainder_over_J2']/L) for L,d in zip(Ls,ds)]}
    out['runs'].append(entry)
    assert all(max(d['residuals'])<1e-7 for d in ds)
    print(json.dumps(entry,indent=2))
out['interpretation']='Two independent numerical diagnostics support one charge and one payload boson at fixed J. Sizes are small; no thermodynamic theorem or violating central charge is established.'
(HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
