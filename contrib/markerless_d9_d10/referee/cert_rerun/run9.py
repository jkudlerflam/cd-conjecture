"""Markerless 9-letter graph: C84's integer certificate pipeline with the marker vertex deleted
(n+1 express vertices instead of n+2).  Usage: run9.py n eta_power"""
import os, sys, json
os.environ['OPENBLAS_NUM_THREADS']='1'
from certify_tree import certificate
from rational_bands import derivative, bands_from
from full_loop_integer import compute, analyze
from certify_roots import main as roots
n = int(sys.argv[1]); eta_power = int(sys.argv[2])
c = certificate(n); print(json.dumps(c), flush=True)
assert c['status']=='proved_positive'
json.dump(c, open(f'certificate_n{n}.json','w'), indent=2)
f,a,d,b = derivative(n); bands_from(f,a,d,b,n)
compute(n=n, eta_power=eta_power)
analyze(n)
roots(n=n, width=1e-9)
