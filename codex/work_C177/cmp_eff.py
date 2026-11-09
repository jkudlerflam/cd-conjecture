#!/usr/bin/env python3
"""C177: full-chain open-chain entanglement against C173's effective-model DMRG at the same L (block cuts at 3k).
usage: cmp_eff.py file.json ..."""
import sys, json, numpy as np
EFF = {48: [1.6315, 1.6759, 1.6315], 72: [1.7818, 1.8053, 1.7818], 96: [1.8638, 1.8804, 1.8639]}  # C173 bd_E/bd_K, block cuts L/2-3, L/2, L/2+3
for f in sys.argv[1:]:
    d = json.load(open(f)); L = d['L']; S = np.array(d['S']); c = L // 2
    blk = S[[c - 4, c - 1, c + 2]] if L % 3 == 0 else None
    P = np.array(d.get('P', [0])); PP = np.array(d.get('PP', [np.nan]))
    print(f"{d['tag']} L={L} Sbar(3 mid sites)={d['Sbar']:.4f} S(L/2-3,L/2,L/2+3)={np.round(blk, 4).tolist() if blk is not None else None} "
          f"mean={np.mean(blk) if blk is not None else np.nan:.4f} eff={EFF.get(L)} c_prof={d['c_prof']:.3f} trunc={d['trunc']:.1e} "
          f"n_mid={d['n_mid']:.4f} |P|max={np.max(np.abs(P)):.1e} PP(L/4)={PP[L // 4 - 2] if PP.size > L // 4 else np.nan:+.3e} "
          f"PP(end)={PP[-3] if PP.size > 3 else np.nan:+.3e}")
