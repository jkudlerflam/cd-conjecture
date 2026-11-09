#!/usr/bin/env python3
"""C170: ground-state energy per site of the C165 parent (Delta = 2, Pauli units) in a field hz, from ring ED
sectors (min over magnetisation), to compare with C165's VUMPS 'gapped Neel' energies at hz = 1.0, 2.2, 2.6."""
import numpy as np
from plateau import ed_sector
L = 20; D = 2.0
E = {M: ed_sector(L, D, M) for M in range(L // 2 - 7, L // 2 + 1)}
for hz in (0.6, 1.0, 1.4, 1.8, 2.2, 2.6):
    vals = {M: (E[M] + hz * (2 * M - L)) / L for M in E}   # M up spins, field +hz sum mu^z
    Mb = min(vals, key=vals.get)
    print(f'hz={hz}: parent ring L={L} e = {vals[Mb]:.5f} at m = {(2*Mb-L)/L:+.2f}   (Neel sector e = {vals[L//2]:.5f})')
