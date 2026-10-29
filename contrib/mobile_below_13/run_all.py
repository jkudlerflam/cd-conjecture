"""Run every check in this folder (all local, numpy/scipy only).
Approximate wall times on a laptop: parent, rings, identity, obstruction,
vacancy_scales under a minute each; perturb about 2 minutes; lift about 15
minutes; dictionary about 15 minutes."""
import subprocess, sys
for s in ['check_parent.py', 'check_rings.py', 'check_frozen_identity.py', 'obstruction.py',
          'perturb.py', 'vacancy_scales.py', 'check_lift.py', 'check_bitblind.py', 'check_dictionary.py']:
    print('==', s, flush=True)
    subprocess.run([sys.executable, '-B', s], check=True)
