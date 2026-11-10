"""C190 validation of z2dmrg.py against exact diagonalization (Claude Opus 5.5 agent, 2026-11-09).

For each test Hamiltonian at small L, both parity sectors: ground energy, entropy at every cut, <Z_j>,
P_bond, one-body matrices G_b and G_f on a window, and C_P(r), compared with sparse ED of the ORIGINAL
(unrotated) Hamiltonian built independently from Pauli Kronecker products.
"""
import sys
import time
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import z2dmrg as Z
import models as Mo


def site_op(op, j, L):
    return sps.kron(sps.kron(sps.identity(2 ** j), sps.csr_matrix(op)), sps.identity(2 ** (L - j - 1)), format='csr')


def ed_ham(L, h2s, h1s, nnn=None):
    H = sps.csr_matrix((2 ** L, 2 ** L), dtype=complex)
    for j in range(L - 1):
        H = H + sps.kron(sps.kron(sps.identity(2 ** j), sps.csr_matrix(h2s[j])), sps.identity(2 ** (L - j - 2)))
    for j in range(L):
        H = H + site_op(h1s[j], j, L)
    if nnn is not None:
        for j in range(L - 2):
            tpm, tmp = nnn[j]
            H = H + tpm * site_op(Z.SP, j, L) @ site_op(Z.SM, j + 2, L) + tmp * site_op(Z.SM, j, L) @ site_op(Z.SP, j + 2, L)
    return H.tocsr()


def ed_sector(H, L, P):
    idx = np.array([i for i in range(2 ** L) if bin(i).count('1') % 2 == P])
    Hs = H[idx][:, idx]
    e, v = spla.eigsh(Hs, k=1, which='SA', tol=1e-13)
    psi = np.zeros(2 ** L, dtype=complex)
    psi[idx] = v[:, 0]
    return e[0], psi


def ed_obs(psi, L, gwin_lo, gwin):
    S = []
    for l in range(1, L):
        s = np.linalg.svd(psi.reshape(2 ** l, -1), compute_uv=False) ** 2
        S.append(Z.entropy(s))
    Zs = [np.real(np.vdot(psi, site_op(Z.ZZ_, j, L) @ psi)) for j in range(L)]
    Pb = []
    for j in range(L - 1):
        O = 2 * (site_op(Z.SP, j, L) @ site_op(Z.SM, j + 1, L) + site_op(Z.SM, j, L) @ site_op(Z.SP, j + 1, L))
        Pb.append(np.real(np.vdot(psi, O @ psi)))
    Gb = np.zeros((gwin, gwin), dtype=complex)
    Gf = np.zeros((gwin, gwin), dtype=complex)
    for a in range(gwin):
        for b in range(gwin):
            i, j = gwin_lo + a, gwin_lo + b
            Ob = site_op(Z.SM, i, L) @ site_op(Z.SP, j, L)
            Gb[a, b] = np.vdot(psi, Ob @ psi)
            if i < j:
                Of = site_op(Z.SM, i, L)
                for k in range(i + 1, j):
                    Of = Of @ site_op(Z.ZZ_, k, L)
                Of = Of @ site_op(Z.SP, j, L)
            elif i > j:
                Of = site_op(Z.SP, j, L)
                for k in range(j + 1, i):
                    Of = Of @ site_op(Z.ZZ_, k, L)
                Of = site_op(Z.SM, i, L) @ Of     # c+_i c_j for i > j: sm_i Z..Z sp_j with the same string
            else:
                Of = Ob
            Gf[a, b] = np.vdot(psi, Of @ psi)
    return np.array(S), np.array(Zs), np.array(Pb), Gb, Gf


def ed_cp(psi, L, i0):
    def Pop(j):
        return 2 * (site_op(Z.SP, j, L) @ site_op(Z.SM, j + 1, L) + site_op(Z.SM, j, L) @ site_op(Z.SP, j + 1, L))
    P0 = Pop(i0)
    e0 = np.real(np.vdot(psi, P0 @ psi))
    out = []
    for k in range(i0 + 2, L - 1):
        Pk = Pop(k)
        out.append((k - i0, np.real(np.vdot(psi, P0 @ (Pk @ psi))) - e0 * np.real(np.vdot(psi, Pk @ psi))))
    return out


def run_case(label, mpo, gph, h2s_orig, h1s_orig, L, nnn=None, chi=64):
    H = ed_ham(L, h2s_orig, h1s_orig, nnn)
    worst = {}
    for P in (0, 1):
        e_ed, psi_ed = ed_sector(H, L, P)
        S_ed, Z_ed, P_ed, Gb_ed, Gf_ed = ed_obs(psi_ed, L, 2, 6)
        cp_ed = ed_cp(psi_ed, L, 2)
        dt = complex if mpo.dtype is complex else float
        res = []
        for init in ('rand', 'prod'):
            if init == 'rand':
                psi = Z.random_mps(L, P, D0=4, seed=7 + P, dtype=dt)
            else:
                occ = [1, 5] if P == 0 else [3]
                psi = Z.product_mps(L, occ, dtype=dt)
            sw = [dict(chi=16, alpha=1e-3, res_tol=1e-8), dict(chi=32, alpha=1e-5, res_tol=1e-9),
                  dict(chi=chi, alpha=0.0, res_tol=1e-11), dict(chi=chi, alpha=0.0, res_tol=1e-12),
                  dict(chi=chi, alpha=0.0, res_tol=1e-12)]
            psi, hist = Z.dmrg(psi, mpo, sw)
            assert psi.check(), 'parity structure broken'
            E = Z.energy(psi, mpo) + mpo.const
            var, _ = Z.variance(psi, mpo)
            S = Z.entropies(psi)
            m = Z.measure(psi, mpo, gph=gph, gwin=6, cp_i0=2, cp_rmax=L)
            m['Gb_re']
            Gb = np.array(m['Gb_re']) + 1j * np.array(m['Gb_im'])
            Gf = np.array(m['Gf_re']) + 1j * np.array(m['Gf_im'])
            cp = np.array([v for r, v in m['CP']])
            par = Z.parity_expectation(psi)
            errs = dict(E=abs(E - e_ed), S=np.max(np.abs(S - S_ed)), Z=np.max(np.abs(np.array(m['Z']) - Z_ed)),
                        Pb=np.max(np.abs(np.array(m['P_bond']) - P_ed)), Gb=np.max(np.abs(Gb - Gb_ed[:, :])),
                        Gf=np.max(np.abs(Gf - Gf_ed)), CP=np.max(np.abs(cp - np.array([v for r, v in cp_ed]))),
                        par=abs(par - (1 - 2 * P)), var=abs(var), Ehist=abs(hist[-1]['E'] + mpo.const - e_ed))
            res.append((init, errs))
            for k, v in errs.items():
                worst[k] = max(worst.get(k, 0.0), v)
            print(f"{label} P={P} init={init}: E_ed={e_ed:.12f} " +
                  ' '.join(f"{k}={v:.1e}" for k, v in errs.items()), flush=True)
    return worst


if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    t0 = time.time()
    allw = {}
    # 1. phi = 1.53 bond (complex), no gauge
    for name, gauge in (('p153', None), ('rat', None), ('sym40', 'stag'), ('sym40', None)):
        b = Mo.bond_family(name)
        hz = b['hz0']
        mpo, gph, shift, info = Mo.candidate(name, L, hz, gauge=gauge)
        h1o = [hz * Z.ZZ_ + shift * Z.I2] * L
        w = run_case(f"{name}/{gauge}", mpo, gph, [info['h2']] * (L - 1), h1o, L)
        allw[f"{name}/{gauge}"] = w
    # 2. controls
    n_op = np.diag([0.0, 1.0])
    for legs in (1, 2):
        mpo, gph = Mo.xx_control(L, mu=-1.2, legs=legs)
        h2o = [-(np.kron(Z.SP, Z.SM) + np.kron(Z.SM, Z.SP))] * (L - 1) if legs == 1 else [np.zeros((4, 4))] * (L - 1)
        nnn = None if legs == 1 else [(-1.0, -1.0)] * (L - 2)
        w = run_case(f"xx{legs}", mpo, gph, h2o, [1.2 * n_op] * L, L, nnn=nnn)
        allw[f"xx{legs}"] = w
    print('worst errors per case:')
    for k, w in allw.items():
        print(' ', k, ' '.join(f"{a}={v:.1e}" for a, v in w.items()))
    print(f'time {time.time() - t0:.1f}s')
