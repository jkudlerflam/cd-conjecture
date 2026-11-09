"""C168, Lemma D check (c): scaling of the exchange (passing) energy at low density.
Same model as imp_swap.py. At eta = 0 the swap lowers the ground energy by dE_x(L) = E0(s) - E0(0);
Lemma D uses dE_x ~ s n^3 (transmission ~ k at low relative momentum), to be compared with the
mass-imbalance ordering energy ~ eta^2 n^2 of imp_mix.py. Also reports w_alt at eta = 0.5, s = 0.2.
Claude Opus 5.5, 2026-11-08."""
import numpy as np, scipy.sparse.linalg as spla
from imp_swap import build
for L in (12, 16, 20, 24, 28, 32):
    out = []
    for (eta, s) in ((0.0, 0.0), (0.0, 0.2), (0.5, 0.2)):
        H, alt = build(L, 1 + eta, 1 - eta, s)
        w, v = spla.eigsh(H, k=1, which='SA', tol=1e-12)
        out.append((w[0], np.sum(np.abs(v[:, 0]) ** 2 * alt)))
    dEx = out[1][0] - out[0][0]
    print(f'L {L:3d}  dE_x(s=0.2) {dEx:+.4e}  dE_x*L^3 {dEx*L**3:+.4f}  dE_x*L^2 {dEx*L**2:+.4f}  w_alt(eta=0,s=.2) {out[1][1]:.4f}  w_alt(eta=.5,s=.2) {out[2][1]:.4f}', flush=True)
