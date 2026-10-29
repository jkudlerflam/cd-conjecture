"""C84, thirteen active flavours, fourteen data sites, no E_1 marker."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from certify_tree import certificate
from rational_bands import derivative,bands_from
from full_loop_integer import compute,analyze
from certify_roots import main as roots
from finite_budget import main as budget
if __name__=='__main__':
 n=13
 c=certificate(n)
 assert c['status']=='proved_positive'
 json.dump(c,open(f'certificate_n{n}.json','w'),indent=2)
 f,a,d,b=derivative(n)
 bands_from(f,a,d,b,n)
 compute(n=n,eta_power=66)
 analyze(n)
 roots(n=n,width=1e-9)
 budget()
