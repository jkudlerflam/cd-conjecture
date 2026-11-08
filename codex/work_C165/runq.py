#!/usr/bin/env python3
"""C165 driver: build (J, B) from classical.fam keyword arguments and run fesx.py (two-site VUMPS, FES)
usage: runq.py 'a=0.8,Dz=1.02,g=0.3,...' out.json chi1 chi2 ... [--n 2]"""
import sys, os, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classical import fam, nfloat
kw = eval('dict(' + sys.argv[1] + ')'); out = sys.argv[2]
args = sys.argv[3:]; n = '2'
if '--n' in args:
    i = args.index('--n'); n = args[i + 1]; del args[i:i + 2]
J, B = fam(**kw) if 'nfloat' not in kw else nfloat(**kw['nfloat'])
cmd = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fesx.py'), 'gen', '--J'] + \
      [f'{x:.12g}' for x in J.ravel()] + ['--B'] + [f'{x:.12g}' for x in B] + ['--n', n, '--chis'] + args + ['--out', out]
sys.exit(subprocess.call(cmd))
