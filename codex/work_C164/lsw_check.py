"""C164: check clat.lsw_freq against the exact one-magnon band of rings where the polarized state is exact
(U(1) about z: XXZ + D_z + field along z, strong field), using est.py ring ED. Fixes the momentum convention."""
import numpy as np, clat, est
rng = np.random.default_rng(3)
L = 10
for trial in range(4):
    Jp, Jz, Dz = rng.normal(size=3) * 0.5
    fz = -4.0
    v = np.array([Jp, Jp, Jz, 0, 0, Dz, 0, 0, fz])
    J, f = clat.J_f(v)
    n, e = clat.uniform_min(J, f)
    ks = 2 * np.pi * np.arange(L) / L
    om, mineig = clat.lsw_freq(J, f, n, ks)
    h = est.bond_from_JB(J, f)
    T = est.table(L)
    E0 = None; ex = []
    for m in range(L):
        w, _ = est.lowest(T.sector_matrix(h, m), 3)
        ex.append(w)
    ex = np.array(ex); E0 = ex.min()
    first = np.array([sorted(ex[m] - E0)[1] if m == 0 else ex[m, 0] - E0 for m in range(L)])
    print(f"stable={mineig.min()>0} " + f'n={np.round(n,3)} max|om(k)-ED(k)|={np.abs(om-first).max():.2e}  max|om(-k)-ED(k)|={np.abs(om[(-np.arange(L))%L]-first).max():.2e}')
