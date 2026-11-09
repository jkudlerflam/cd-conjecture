"""Ring ED map of the parity class slice phi = phi0, h_z = 0 (mu = V), fermion units t = 1.
Per point and ring size: lowest 3 levels in each Pi sector, half-ring entropy of the global ground state.
Usage: python3 ring_map.py phi L1,L2 out.jsonl"""
import sys, json, time
import numpy as np
from model import bond, spin_from_fermion, ring_data
phi = float(sys.argv[1]); Ls = [int(x) for x in sys.argv[2].split(',')]; fout = sys.argv[3]
Ds = [float(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else [0.0,0.1,0.2,0.3,0.4,0.5,0.6]
Vs = [float(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else [-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2]
with open(fout, 'a') as f:
    for Dl in Ds:
        for V in Vs:
            p = spin_from_fermion(1.0, phi, Dl, V, V)
            h = bond(*p)
            for L in Ls:
                r = ring_data(h, L, k=3)
                f.write(json.dumps(dict(phi=phi, Delta=Dl, V=V, L=L, p=[float(x) for x in p], **r)) + '\n'); f.flush()
