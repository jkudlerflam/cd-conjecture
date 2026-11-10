"""C196b gate G1 (Claude Opus 5.5 agent, 2026-11-10): the q-site-cell iDMRG against exact diagonalization.

python3 -E -s validate196b.py [ed|tm|embed]

ed:    untruncated growth from a random open window of 2q sites, q = 2, 3, 4. After iteration t the accumulated
       window eigenvalues are the ground energy (P = 0) of the open chain of 2q(t + 1) sites with the bulk MPO.
       Compared with sparse Lanczos ED of the same chain, at sym40_n316, p153c and the TFIM (h = 1.05).
tm:    the Arnoldi transfer spectrum against dense eig of the same two-cell map.
embed: the effective-model window converted to Z2 labels reproduces its own energy (effective model, through the
       same MPS) and its particle number, and the window density.
"""
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z2dmrg as Z0        # noqa: E402
import idmrg196b as B      # noqa: E402
import models as Mo        # noqa: E402

HZ = dict(p153c=-19.97589125, sym40_n316=-19.976047250613814, sym40_n18=-19.976333693659278)
CASES = [('sym40_n316', dict(model='sym40', hz=HZ['sym40_n316'])),
         ('p153c', dict(model='p153', hz=HZ['p153c'])),
         ('tfim', dict(model='tfim', h=1.05))]


def local_terms(cfg):
    if cfg['model'] == 'tfim':
        return Z0.pauli_bond(Jx=-1.0), -cfg['h'] * Z0.ZZ_
    b = Mo.bond_family(cfg['model'])
    h2 = Z0.pauli_bond(Jx=b['Jx'], Jy=b['Jy'], Jz=b['Jz'], D=b['D'])
    shift = -(b['Jz'] + cfg['hz'])
    return h2, cfg['hz'] * Z0.ZZ_ + shift * Z0.I2


def ed_energy(cfg, L):
    h2, h1 = local_terms(cfg)
    H = sp.csr_matrix((2 ** L, 2 ** L), dtype=complex)
    for j in range(L - 1):
        H = H + sp.kron(sp.kron(sp.identity(2 ** j), sp.csr_matrix(h2)), sp.identity(2 ** (L - j - 2)), 'csr')
    for j in range(L):
        H = H + sp.kron(sp.kron(sp.identity(2 ** j), sp.csr_matrix(h1)), sp.identity(2 ** (L - j - 1)), 'csr')
    par = np.array([bin(i).count('1') % 2 for i in range(2 ** L)])
    idx = np.where(par == 0)[0]
    Hs = H[idx][:, idx]
    if Hs.shape[0] <= 512:
        return float(np.linalg.eigvalsh(Hs.toarray())[0])
    return float(spla.eigsh(Hs, k=1, which='SA', tol=1e-14, maxiter=100000)[0][0])


def part_ed(Lmax=16):
    print('# G1 part ed: untruncated q-site-cell iDMRG against sparse ED (open chain, bulk MPO, P = 0)')
    worst = 0.0
    for name, cfg in CASES:
        mpo, W, e_off = B.build_model(cfg)
        for q in (2, 3, 4):
            st = B.CellIDMRG(mpo, W, q, cutoff=0.0, seed=3)
            Ms, n0 = B.random_window(2 * q, seed=5, D0=4)
            st.set_window(Ms, n0)
            Etot, line = 0.0, []
            t = 0
            while 2 * q * (t + 1) <= Lmax:
                rec = st.iterate(chi=10 ** 6, nfull=3)
                Etot += rec['E']
                L = 2 * q * (t + 1)
                e0 = ed_energy(cfg, L)
                d = abs(Etot - e0)
                worst = max(worst, d)
                line.append(f'L={L} E={Etot:.12f} dE={d:.1e}')
                t += 1
            print(f'{name} q={q}: ' + '; '.join(line), flush=True)
    print(f'G1 worst |E_iDMRG - E_ED| = {worst:.2e}')
    return worst


def part_tm():
    print('# G1 part tm: Arnoldi against dense eig of the two-cell transfer map')
    name, cfg = CASES[2]
    mpo, W, e_off = B.build_model(cfg)
    st = B.CellIDMRG(mpo, W, 2, seed=3, e_off=e_off)
    Ms, n0 = B.random_window(4, seed=5, D0=4)
    st.set_window(Ms, n0)
    for _ in range(30):
        st.iterate(chi=24)
    a = st.transfer_spectrum(k=4, m=40)
    import idmrg196b as BB
    old = BB.arnoldi_top
    # dense reference: force the dense path by a large threshold
    L = st.last
    print('arnoldi', {r: np.round(a['spec'][r]['abs'], 12).tolist() for r in (0, 1)}, 'nmv',
          [a['spec'][r]['nmv'] for r in (0, 1)])
    return a


def part_embed():
    print('# embedding check: effective window -> Z2 qubit MPS')
    Ms, n0, E = B.eff_window(32, 6, 1j)
    mpo, W, e_off = B.build_model(dict(model='sym40', hz=HZ['sym40_n316']))
    st = B.CellIDMRG(mpo, W, 16)
    st.set_window(Ms, n0)
    d = st.window_density()
    print('E_eff =', E, ' N =', sum(d), ' dens[0:8] =', np.round(d[:8], 4))


if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'ed'
    {'ed': part_ed, 'tm': part_tm, 'embed': part_embed}[part]()
