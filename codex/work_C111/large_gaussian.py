"""Large XX checks. Refuses execution outside a SLURM compute allocation."""
import json
import os
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh_tridiagonal


def main():
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURMD_NODENAME'):
        raise RuntimeError('Run only as a SLURM job on a Typhon compute node.')
    if 'login' in os.uname().nodename.lower():
        raise RuntimeError('Refusing a login node.')
    rows = []
    for n in [512, 1024, 2048, 4096, 8192, 16384]:
        x = np.arange(1, n, dtype=float)
        w = x*(n-x)/n
        ev = eigvalsh_tridiagonal(np.zeros(n), w/2, select='i',
                                  select_range=(0, n//2-1))
        defect_xx = -w.sum()/np.pi-ev.sum()
        q = np.sum(np.diff(w)**2/w[:-1])
        for delta in [.05, .01, .001]:
            right, left = np.r_[w, 0.], np.r_[0., w]
            u, v = (right+left)/2, right-left
            ds = np.sum(.01*v*v/(np.hypot(delta*u, .1*v)+delta*u))
            rows.append(dict(n=n, delta=delta, xx_deficit=float(defect_xx),
                             spectator_trial_deficit=float(ds),
                             total_deficit_lower_bound=float(defect_xx+ds),
                             Q=float(q), predicted_spectator_log_slope=.01/delta))
    Path('large_gaussian_results.json').write_text(json.dumps(rows, indent=2)+'\n')


if __name__ == '__main__':
    main()
