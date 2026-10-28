#!/usr/bin/env python3
"""Exact open-chain entropies of the quadratic (free and Hartree-Fock) chains, Calabrese-Cardy fits at
several L: the calibration of the open-chain estimator.  usage: exact_open_table.py SPEC ...
Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28."""
import os, sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_chain as R
from analyze import cc_fit
for spec in sys.argv[1:]:
    K, meta = R.build(spec)
    row = dict(spec=spec)
    for L in (64, 96, 128, 256):
        ex = R.exact_open(K, L)
        row[str(L)] = cc_fit(ex['S'], L)
    print(json.dumps(row), flush=True)
