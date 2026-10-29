"""New n=16 band certificate using attributed local copies of C75 routines."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
from certify_tree import certificate
from rational_bands import derivative
from full_loop_integer import compute, analyze
from certify_roots import main as roots

if __name__ == '__main__':
    n = 16
    c = certificate(n)
    assert c['status'] == 'proved_positive'
    with open('certificate_n16.json', 'w') as f:
        json.dump(c, f, indent=2)
    derivative(n)
    compute(n=n, eta_power=92)
    analyze(n)
    roots(n=n, width=1e-9)
