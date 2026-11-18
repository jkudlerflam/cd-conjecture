"""Non-perturbative check: exact diagonalisation of two padded cells of the 9-letter chain (n=2, 3600 states,
all head sectors) at finite port lam and conversion g.  The lowest 4^n levels are compared with the
free-fermion prediction E0 + lam^2 h0 + lam^2 g^2 * (one-body spectrum filled in all ways)."""
import numpy as np, itertools, json
from twocell9 import build, analyse, one_particle
n = 2
base = dict(Delta=2.2-2*np.cos(np.pi/(n+2)), D=0.4, nu=0.0, f=0.0, phiA=np.pi/(2*(n+1)), td=0.3)
ref = analyse(n, dict(base, lam=1.0, g=1.0))
h, vacE2 = one_particle(ref['h2'], n, ref['bitlist'])
eps = np.linalg.eigvalsh((h+h.conj().T)/2)                       # one-body levels, coefficient of lam^2 g^2
E0 = 2*(-2*np.cos(np.pi/(n+2))); h0 = ref['h0'][0,0].real
occ = np.array(list(itertools.product((0,1), repeat=2*n)))
many = np.sort(occ @ eps)                                        # free-fermion many-body levels (coefficient)
out = []
for lam, g in [(0.2,0.2),(0.1,0.2),(0.05,0.2),(0.1,0.1),(0.05,0.1),(0.05,0.05)]:
    states, H, leaks = build(n, dict(base, lam=lam, g=g))
    ev = np.linalg.eigvalsh(H.toarray())[:4**n]
    pred = E0 + lam**2*(h0 + g**2*(vacE2.real + many))
    shift = ev.mean() - pred.mean()
    rel = np.abs((ev-ev.mean()) - (pred-pred.mean())).max()/(lam**2*g**2*np.ptp(many))
    out.append(dict(lam=lam, g=g, leaks=leaks, common_shift=shift, splitting_error_rel=rel, gap_to_next=None))
    print(f"lam={lam} g={g}: leaks={leaks} common shift={shift:.3e} (lam^4={lam**4:.1e}, lam^2 g^4={lam**2*g**4:.1e})  "
          f"max level error / payload bandwidth = {rel:.3e}")
json.dump(out, open('ed_two_cell.json','w'), indent=1)
