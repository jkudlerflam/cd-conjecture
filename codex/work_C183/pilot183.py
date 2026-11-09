#!/usr/bin/env python3
"""C183 local pilot (Claude Opus 5.5 agent, 2026-11-09): density ladder at L = 24 (chi <= 48) for one family, printing
the dressed density n_d = n(central half) - n_vac against x.  usage: pilot183.py r wa phi x1,x2,...  outdir"""
import sys, json, subprocess, numpy as np
import tune183 as TU
r, wa, phi = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
xs = [float(v) for v in sys.argv[4].split(',')]
od = sys.argv[5]
T = TU.tune(r, wa, phi)
json.dump(T, open(f'{od}_tune.json', 'w'), indent=1)
for x in xs:
    cp = TU.couplings(T, x)
    tag = f"lad_r{r:g}_w{wa:g}_p{phi:.3f}_x{x:g}"
    subprocess.run(['python3', 'fc.py', tag, '24', '--cp', json.dumps(cp), '--chis', '16', '32', '48', '48', '--outdir', od, '--nocorr'], capture_output=True)
    d = json.load(open(f'{od}/{tag}_24.json'))
    Z = np.array(d['Z']); L = len(Z)
    n = float(np.mean((1 - Z[L // 4: 3 * L // 4]) / 2))
    print(f"x={x:+.3f} hz={cp['hz']:.6f} n_central={n:.4f} n_d={n - T['n_vac']:.4f} Sbar={d['Sbar']:.4f} t={d['time_s']:.0f}s", flush=True)
