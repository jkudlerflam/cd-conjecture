"""Reproduce final saved-state convergence and centred-matrix diagnostics.

No primary dimension is assigned. Repeated measurements of the unchanged
control are checked separately from changes between optimized chi stages.
"""
import csv
import gzip
import json
from pathlib import Path
import numpy as np
from common import HERE, ROOT, write_json, source_hashes, sha256
from compare_correlators import main as compare
from analyze_matrix import main as matrix
from summarize_runs import summarize


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_text())


def main():
    results = HERE / 'results'
    final = results / 'quartet-final-measure-20261125a' / 'results'
    analysis = HERE / 'analysis'
    states = read(final / 'measured_states.json')
    assert len(states) == 2 and all(x['measurement_complete'] for x in states)
    for x in states:
        assert abs(x['parity_expectation'] - 1) < 1e-9
        assert x['saved_state_norm_error'] < 1e-9
    summarize()
    selected = {}
    for tag, run in [('target', 'quartet-pilot-20261125a'), ('c1control', 'quartet-c1control-20261125a')]:
        stem = f'L128_p0_chi256_{tag}'
        old = results / run / 'results' / f'L128_p0_chi128_{tag}_quartet.json.gz'
        new = final / (stem + '_quartet.json.gz')
        output = analysis / f'{tag}_final_chi128_256_changes.json'
        compare(old, new, output)
        matrix(final / (stem + '_matrix.json.gz'), new,
               analysis / f'{tag}_centered_chi256_diagnostics.json')
        d = read(output)
        selected[tag] = {'kappa_pi': d['kappa_pi'],
            'changes_ref47_width32_r24_44': [x for x in d['filtered_changes']
                                           if x['ref'] == 47 and x['width'] == 32],
            'fixed_momentum_changes_ref47_width32_r24_44': [x for x in d['filtered_fixed_momentum_changes']
                                           if x['ref'] == 47 and x['width'] == 32]}
    # Same immutable control checkpoint, independently loaded and measured.
    before = read(results / 'quartet-c1control-20261125a' / 'results' /
                  'L128_p0_chi256_c1control_quartet.json.gz')
    after = read(final / 'L128_p0_chi256_c1control_quartet.json.gz')
    differences = []
    for kind in ('P', 'Q'):
        for key in before[kind]:
            for ref in before['refs']:
                for label in ('md', 'dm', 'mm', 'dd', 'mz'):
                    a = np.asarray(before[kind][key][str(ref)][label])
                    b = np.asarray(after[kind][key][str(ref)][label])
                    differences.append(float(np.max(abs(a - b))))
    repeat_error = max(differences)
    assert repeat_error < 1e-13, repeat_error
    jobs = []
    for run in sorted(results.iterdir()):
        accounting = run / 'steps_accounting.psv'
        if not accounting.exists():
            continue
        rows = list(csv.DictReader(accounting.read_text().splitlines(), delimiter='|'))
        batch = next(x for x in rows if x['JobIDRaw'].endswith('.batch'))
        assert batch['State'] == 'COMPLETED' and batch['ExitCode'] == '0:0'
        jobs.append({'run': run.name, **batch})
    write_json(analysis / 'final_summary.json', {
        'source_sha256': source_hashes(),
        'input_sha256': {str(p.relative_to(ROOT)): sha256(p)
                         for p in sorted(final.glob('*')) if p.is_file()},
        'jobs': jobs,
        'states': [{k: x[k] for k in ('tag', 'state', 'state_sha256',
                                    'parity_expectation', 'saved_state_norm_error')} for x in states],
        'same_control_checkpoint_repeat_max_abs': repeat_error,
        'stage_comparisons': selected,
        'interpretation': 'Conditional full-scale production gate not passed: target chi256 is unconverged, and neutral identification/filter systematics remain unresolved. No charge-four dimension or locking conclusion.'})
    print('Same-checkpoint control repeat maximum absolute difference:', repeat_error)
    for tag, d in selected.items():
        print(tag, 'kappa/pi', d['kappa_pi'], 'relative L2 changes:',
              {x['kind'] + ':' + x['operator']: round(x['relative_l2_change'], 6)
               for x in d['changes_ref47_width32_r24_44']})


if __name__ == '__main__':
    main()
