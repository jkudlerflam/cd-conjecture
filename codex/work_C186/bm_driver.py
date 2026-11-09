"""C186 Typhon driver: tune (phi, mu_onset) of the mixed cusp model on an N-capped ring (pairbind2.Tuner, L = 48,
N_max = 4), then run the block DMRG (blockmix.py) of the open chain at mu = mu_onset + dmu for Lb blocks.
Usage: bm_driver.py u v rmode dmu Lb chis tag"""
import numpy as np, json, sys, time
from pairbind2 import Tuner
from blockmix import block_model, density_profile, LOC
from dmrg_d import dmrg
u, v, rmode, dmu, Lb = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3], float(sys.argv[4]), int(sys.argv[5])
chis = tuple(int(x) for x in sys.argv[6].split(','))
tag = sys.argv[7]
t0 = time.time()
tu = Tuner(48, 4, u, v, rmode)
phi, mu0, ea, ia, eb, ib = tu.tune()
r = 0.0 if rmode == "zero" else np.sqrt(v)
hb, hs = block_model(u * np.exp(1j * phi), V2=-2.0 + v, mu=mu0 + dmu, U=40.0, r=r)
out = dmrg(hb, Lb, chis=chis, cutoff=1e-11, seed=0, onsite=[hs] * Lb)
nloc = np.array([sum(s) for s in LOC], float)
dens = density_profile(out['M'], nloc)
S = out['S']
rec = dict(u=u, v=v, rmode=rmode, r=r, phi=phi, mu_onset=mu0, dmu=dmu, mu=mu0 + dmu, Lb=Lb, L=3 * Lb, chis=chis,
           E=out['E'], n_mid=float(dens[Lb // 4: 3 * Lb // 4].mean() / 3),
           Smid=float(np.mean(S[len(S) // 2 - 1: len(S) // 2 + 2])), S=[float(x) for x in S], hist=out['hist'],
           secs=time.time() - t0)
with open(f"bm_{tag}.jsonl", "a") as f:
    f.write(json.dumps(rec) + "\n")
print(json.dumps({k: rec[k] for k in rec if k not in ('S', 'hist')}), flush=True)
