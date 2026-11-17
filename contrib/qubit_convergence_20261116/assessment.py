#!/usr/bin/env python3
"""Assemble the current assessment without substituting old data for new stages."""
import json
import math
from pathlib import Path
from audit import HERE,ROOT,sha


def main():
    a=json.loads((HERE/'audit.json').read_text())
    v=json.loads((HERE/'validation.json').read_text())
    control=json.loads((HERE/'controls.json').read_text())
    resources=json.loads((HERE/'resources.json').read_text())
    corr=json.loads((HERE/'correlations.json').read_text())
    current=[c for c in a['comparisons'] if c['key']=='sym40_n316' and c['Lsmall']==256][0]
    r=next(r for r in a['runs'] if r['key']=='sym40_n316' and r['L']==512 and r['chi']==2048)
    # Two-point linear discarded-weight extrapolation. Preserve raw final S.
    old=next(q for q in a['runs'] if q['file']==r['C4']['predecessor_file'])
    ea=old['checks']['trunc']; eb=r['checks']['trunc']
    extrap={w: current['c_two'][w]+6/math.log(2)*(r['Sbar'][w]-old['Sbar'][w])*eb/(ea-eb) for w in ('1','6','12')}
    out=dict(verdict='UNDECIDED',new_chi2560_production='NOT RUN',
        registered_test='REGISTRATION_20261116.md',registration_sha256=sha(HERE/'REGISTRATION_20261116.md'),
        gate_status=dict(saved_arithmetic='PASS',local_hamiltonian_and_resume=v['validation_status'],
            saved_same_solver_controls_through_256='PASS against exact finite-size entropy',
            new_matched_solver_controls_256_512='NOT EVALUATED',
            new_both_parities_2560='NOT EVALUATED',old_L512_C4='FAIL',new_L512_C4='NOT EVALUATED',
            new_Ks_size_drift='NOT EVALUATED'),
        current_ground_branch=current,current_L512_T1=r['checks'],
        linear_discarded_weight_extrapolation_c_two=extrap,
        extrapolation_note='Heuristic through 1536/2048 at L512 only, fixed L256@1536 reference; neither rigorous lower nor upper bound.',
        blockers=resources['host_observed'],
        physics='Finite-size evidence remains consistent with c about 2. This work establishes neither asymptotic c>1 nor a concrete counterexample.',
        correlation_status=corr['new_256_512_Ks_test'])
    (HERE/'assessment.json').write_text(json.dumps(out,indent=2)+'\n')
    sources=['codex/work_C195/'+f for f in ('models.py','z2dmrg.py','z2dmrg_xp.py','run195.py')]
    sources+=['codex/work_C193/ckpt_meas/'+f for f in ('meas193.py','ana193.py','frozen.json')]
    (HERE/'source_manifest.json').write_text(json.dumps({p:sha(ROOT/p) for p in sources},indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
