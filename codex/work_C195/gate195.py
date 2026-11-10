"""C195 gate comparison (Claude Opus 5.5 agent, 2026-11-10).

Usage: python3 gate195.py <reference.json> <test.json> [<test2.json> ...]

Sweep-by-sweep comparison of run outputs of the same configuration: energy, entropy at every cut (mixer-off
sweeps), <Z_j> profile, residual and matvec counts, then the final exact entropies and density.  The gate of the
spec is |dE| <= 1e-10 and |dS| <= 1e-8 at every sweep.
"""
import sys
import json
import numpy as np

ref = json.load(open(sys.argv[1]))
ok_all = True
for fn in sys.argv[2:]:
    t = json.load(open(fn))
    print(f'== {sys.argv[1]} against {fn}: {len(ref["hist"])} and {len(t["hist"])} sweeps, status {ref.get("status")}, '
          f'{t.get("status")}')
    wE = wS = wZ = 0.0
    for a, b in zip(ref['hist'], t['hist']):
        dE = abs(a['E'] - b['E'])
        dZ = float(np.max(np.abs(np.array(a['Zprof']) - np.array(b['Zprof']))))
        dS = float(np.max(np.abs(np.array(a['Sprof']) - np.array(b['Sprof'])))) if 'Sprof' in a else float('nan')
        wE = max(wE, dE); wZ = max(wZ, dZ)
        if dS == dS:
            wS = max(wS, dS)
        print(f"  chi {a['chi']:4d} alpha {a['alpha']:g}: E {a['E_phys']:.10f} dE {dE:.1e} dS {dS:.1e} dZ {dZ:.1e} "
              f"res {a['res_max']:.1e}/{b['res_max']:.1e} mv {a['mv_mean']:.2f}/{b['mv_mean']:.2f} "
              f"t {a['t']:.0f}/{b['t']:.0f}s dev_used {b.get('dev_used_GB', float('nan')):.2f} GB")
    ok = len(ref['hist']) == len(t['hist']) and wE <= 1e-10 and wS <= 1e-8
    print(f'  worst dE {wE:.2e}, worst dS (any cut, mixer-off sweeps) {wS:.2e}, worst dZ {wZ:.2e}')
    if ref.get('status') == 'done' and t.get('status') == 'done':
        dSf = float(np.max(np.abs(np.array(ref['S']) - np.array(t['S']))))
        dn = float(np.max(np.abs(np.array(ref['meas']['dens']) - np.array(t['meas']['dens']))))
        print(f'  final exact entropies max diff {dSf:.2e}, density max diff {dn:.2e}, Sbar6 {ref["Sbar6"]:.8f} '
              f'{t["Sbar6"]:.8f}')
        ok = ok and dSf <= 1e-8
    print('  GATE', 'PASS' if ok else 'FAIL')
    ok_all = ok_all and ok
print('ALL PASS' if ok_all else 'NOT ALL PASS')
