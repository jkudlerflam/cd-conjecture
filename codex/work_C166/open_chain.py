#!/usr/bin/env python3
"""C166 tasks 3 and 4 (Claude Opus 5.5 agent, 2026-11-08): open-chain DMRG (dmrg.py, unchanged, through
c166lib.dmrg_edges) for a named bond at one length, with optional extra edge fields, and the C166 observables:
S(l) profile, c_prof, c_sm, S_mid (dmrg_c.profile_fits), energy and truncation history, <sigma^a_j>,
bond energies <h_{j,j+1}> and their mean over the central L/4 bonds (emid), and the correlation matrix
<sigma^a_{L/4} sigma^b_j> for j = L/2, 3L/4 - 1, 3L/4.
usage: open_chain.py NAME L --chis 16 32 ... [--left hx hy hz] [--right hx hy hz] [--tag T] [--outdir oc]"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import c166lib as X
import dmrg_c as DC

ap = argparse.ArgumentParser()
ap.add_argument('name'); ap.add_argument('L', type=int)
ap.add_argument('--chis', type=int, nargs='+', default=[16, 32, 64, 96, 128, 128])
ap.add_argument('--left', type=float, nargs=3, default=None)
ap.add_argument('--right', type=float, nargs=3, default=None)
ap.add_argument('--tag', default='')
ap.add_argument('--outdir', default='oc')
ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
fn = os.path.join(a.outdir, f"{a.name}{('_' + a.tag) if a.tag else ''}_{a.L}.json")
t0 = time.time()
r, h, v = X.dmrg_edges(a.name, a.L, tuple(a.chis), left=a.left, right=a.right, seed=a.seed, verbose=True)
S = r['S']
cp, cs, mid = DC.profile_fits(S, a.L)
out = dict(name=a.name, L=a.L, chis=a.chis, left=a.left, right=a.right, tag=a.tag, v=v.tolist(), E=r['E'],
           hist=r['hist'], trunc=r['hist'][-1][2], c_prof=cp, c_sm=cs, Smid=mid, S=[float(x) for x in S],
           t_dmrg=time.time() - t0)
out.update(X.measure(r['M'], h, a.L))
out['time_s'] = time.time() - t0
json.dump(out, open(fn + '.tmp', 'w'), default=float); os.replace(fn + '.tmp', fn)
C = np.array(out['corr'])
print(f"done {fn} E={r['E']:.10f} c_prof={cp:.3f} Smid={mid:.4f} trunc={out['trunc']:.1e} emid={out['emid']:.10f} "
      f"S(1)={S[0]:.4f} |m|max={np.max(np.abs(out['mag'])):.4f} zz(L/4,3L/4)={C[2, 2, -1]:+.4f} "
      f"t={out['time_s']:.0f}s", flush=True)
