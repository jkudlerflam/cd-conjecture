"""Read-only provenance and observable audit of the actual p153c records."""
from pathlib import Path
import hashlib
import json
import math

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def read(path):
    data=(ROOT/path).read_bytes()
    return dict(path=path,sha256=hashlib.sha256(data).hexdigest()),data


def main():
    rows=[]
    for path in ('codex/work_C195/typhon_out/out512/p153c_L512_P0_rand1.json',
                 'codex/work_C195/typhon_out_x2048/p153c_L512_P0_rand1.json'):
        row,data=read(path); r=json.loads(data); m=r['meas']
        row.update(cfg=r['cfg'],measurement_keys=list(m),one_body_window=m['gwin'],
                   CP_anchor=m['CP_i0'],CP_distance_range=[m['CP'][0][0],m['CP'][-1][0]],
                   CP_operator='P_j=XX+YY=2(S+_j S-_(j+1)+S-_j S+_(j+1)); connected neutral bond correlator.',
                   has_four_time_source_response=False,
                   has_amputating_dynamic_single_response=False)
        rows.append(row)
    sources=[]
    for path in ('codex/work_C195/models.py','codex/work_C195/z2dmrg.py',
                 'contrib/qubit_stable_two_modes_20261118/quartet.py',
                 'contrib/qubit_finite_density_modes_20261118/interface.json'):
        row,data=read(path);sources.append(row)
    w=complex(math.cos(1.53),math.sin(1.53))/80
    coeffs=dict(XX=(1+w.real)/2,YY=(-1+w.real)/2,XY=w.imag/2,
                YX=-w.imag/2,ZZ=10,Z=-19.97589125)
    result=dict(date='2026-11-18',records=rows,sources=sources,
        fixed_microscopic_coefficients=coeffs,
        literal_monomials=6,independent_real_controls_before_scale=5,
        independent_real_controls_after_scale=4,
        nonzero_coefficient_spread=max(map(abs,coeffs.values()))/min(map(abs,coeffs.values())),
        absent_information=['dynamic connected six-insertion source response',
          'full interacting external propagators for amputation',
          'local matching and UV contact subtraction',
          'single-field and fused quartet vertex normalization',
          'dressed covariance and velocities with errors',
          'a positive thermodynamic mismatch lower bound'],
        conclusion='Equal-time one-body and neutral-bond records do not determine the connected static quartet response; no missing correlator is inferred by Wick factorization in this interacting medium.')
    (HERE/'source_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print('audited',len(rows),'p153c records and',len(sources),'source files')
    print('microscopic coefficient spread',result['nonzero_coefficient_spread'])


if __name__=='__main__': main()
