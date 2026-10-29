"""C83 n=13 certificate using C80 integer/rational source copies."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from certify_tree import certificate
from rational_bands import derivative
from full_loop_integer import compute, analyze
from certify_roots import main as roots
if __name__=='__main__':
    n=13
    c=certificate(n)
    assert c['status']=='proved_positive'
    with open('certificate_n13.json','w') as f: json.dump(c,f,indent=2)
    derivative(n)
    compute(n=n,eta_power=74)
    analyze(n)
    roots(n=n,width=1e-9)
