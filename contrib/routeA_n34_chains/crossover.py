#!/usr/bin/env python3
"""Exact block entanglement S(l) of the Hartree-Fock free part of the TRUE chains at large l, to
locate the length beyond which S grows as (n/3) ln l.  usage: crossover.py OUT.json lmax SPEC6 ...
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import models34 as Mo
out = []
lmax = int(sys.argv[2])
for spec6 in sys.argv[3:]:
    K, Ws, Th, Mh, info = Mo.base(spec6)
    ells = [2 ** j for j in range(3, int(np.log2(lmax)) + 1)]
    S = Mo.gaussian_block_entropy(Th, Mh, ells, nk=max(65536, 16 * lmax))
    sl = [3 * (S[i + 1] - S[i]) / np.log(2) for i in range(len(S) - 1)]
    print(spec6, 'l', ells, '\n  S', np.round(S, 4).tolist(), '\n  c slopes', np.round(sl, 3).tolist(), flush=True)
    out.append(dict(spec=spec6, ells=ells, S=S, slopes=sl))
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
