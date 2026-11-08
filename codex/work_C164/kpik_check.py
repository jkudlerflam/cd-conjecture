"""C164: the linear-spin-wave symmetry k -> pi - k at a_c = 0 is an exact symmetry of the Pi-symmetric bond
a(XX - YY) + s(XY + YX) + D(XY - YX) + Jz ZZ + f Z. U = I R_{pi/2}(z) prod_{j odd} R_pi(z)_j maps momentum k to
pi - k in the odd-parity sector of Pi = prod sigma^z (one-particle states) and k to -k in the even sector.
Check on rings: parity-resolved sector spectra; a field along x (which breaks Pi) breaks it."""
import numpy as np, est
rng = np.random.default_rng(2)
L = 12
T = est.table(L)
def par_levels(h, m, p):
    s = T.sector(m)
    Hk = T.sector_matrix(h, m).toarray()
    pc = np.array([bin(int(x)).count('1') % 2 for x in s['R']])
    sel = pc == p
    return np.sort(np.linalg.eigvalsh(Hk[np.ix_(sel, sel)]))[:6] if sel.any() else np.array([])
for trial in range(3):
    a, s, D, Jz, fz = rng.normal(size=5)
    J = np.array([[a, s + D, 0], [s - D, -a, 0], [0, 0, Jz]])
    for fx in (0.0, 0.3):
        h = est.bond_from_JB(J, np.array([fx, 0, fz]))
        if fx == 0:
            # parity of a basis state: number of 1 bits (site i is bit L-1-i); Pi is conserved
            pp = (L % 2)
            dodd = max(np.max(np.abs(par_levels(h, m, 1) - par_levels(h, (L // 2 - m) % L, 1))) for m in range(L))
            dev_ = max(np.max(np.abs(par_levels(h, m, 0) - par_levels(h, (-m) % L, 0))) for m in range(L))
            print(f'trial {trial} fx=0: odd sector max|E(m)-E(L/2-m)| = {dodd:.2e}, even sector max|E(m)-E(-m)| = {dev_:.2e}')
        else:
            lev = [np.sort(np.linalg.eigvalsh(T.sector_matrix(h, m).toarray()))[:6] for m in range(L)]
            dv = max(np.max(np.abs(lev[m] - lev[(L // 2 - m) % L])) for m in range(L))
            print(f'trial {trial} fx=0.3 (Pi broken): max|E(m)-E(L/2-m)| all levels = {dv:.2e}')
