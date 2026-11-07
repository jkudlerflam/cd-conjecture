"""Zero-temperature (m = infinity) limit of the majorant constraint:
lim_{beta->inf} (L/beta) log Z~(L) = L (L e0 - E_L) =: D(L), so S pi V'/12 >= sup_{even L} D(L).
e0 exact (Bethe ansatz for XXZ, free fermions for the XY family) or extrapolated (gapped models);
E_L by sparse Lanczos on even rings L = 2..LMAX. Writes data/zeroT.json."""
import os, json, numpy as np
from scipy.integrate import quad
from scipy.sparse.linalg import eigsh
from cflib import ring_bond_ops, sz_sectors
from models import MODELS, get

HERE = os.path.dirname(os.path.abspath(__file__))


def e0_xxz(D):
    if D == 1:
        return 0.25 - np.log(2)
    if D < 1:
        g = np.arccos(D)
        f = lambda x: 2 * np.exp(-2 * g * x) * (-np.expm1(-2 * (np.pi - g) * x)) / ((-np.expm1(-2 * np.pi * x)) * (1 + np.exp(-2 * g * x))) if x > 0 else (np.pi - g) / np.pi
        return D / 4 - np.sin(g) * quad(f, 0, np.inf, epsabs=1e-13, limit=400)[0]
    eta = np.arccosh(D)
    n = np.arange(1, 400)
    return D / 4 - np.sinh(eta) * (0.5 + 2 * np.sum(np.exp(-2 * n * eta) / (1 + np.exp(-2 * n * eta))))


def e0_xy(Jx, Jy, g):
    """H = -sum [Jx sx sx + Jy sy sy] - g sum sz (Pauli), per site"""
    lam = lambda k: np.sqrt((g + (Jx + Jy) * np.cos(k)) ** 2 + (Jx - Jy) ** 2 * np.sin(k) ** 2)
    return -quad(lam, 0, 2 * np.pi, limit=400, points=[np.pi / 2, np.pi, 3 * np.pi / 2])[0] / (2 * np.pi)


E0_EXACT = {'XX': e0_xxz(0.0), 'XXZ05': e0_xxz(0.5), 'HEIS': e0_xxz(1.0), 'XXZ2': e0_xxz(2.0),
            'TFIM': -4 / np.pi, 'XYI': e0_xy(0.25, -0.125, 0.125)}


def ground(h, L, u1):
    ops = ring_bond_ops(h, 2, L)
    H = sum(ops).tocsr()
    if u1:
        # original-frame XXZ conserves total Sz; positive frame conserves staggered Sz: use a generic solve
        pass
    if H.shape[0] <= 4096:
        return np.linalg.eigvalsh(H.toarray())[0]
    return eigsh(H, k=1, which='SA', tol=1e-12)[0][0]


def run(LMAX=18):
    out = {}
    for name in MODELS:
        h, mode, M, A, label = get(name)
        Ls = list(range(2, LMAX + 1, 2))
        E = [float(ground(h, L, mode == 'u1')) for L in Ls]
        if name in E0_EXACT:
            e0 = E0_EXACT[name]; src = 'exact'
        else:   # gapped: E_L/L converges exponentially
            e0 = E[-1] / Ls[-1]; src = 'E_%d/%d' % (Ls[-1], Ls[-1])
        D = [L * (L * e0 - EL) for L, EL in zip(Ls, E)]
        out[name] = dict(label=label, M=M, A=A, e0=e0, e0_source=src, L=Ls, E=E, D=D, Dmax=max(D), Dmax_L=Ls[int(np.argmax(D))],
                         Dmax_over_M=max(D) / M, cft_value_over_M=np.pi * A / 6 / M)
        print('%-8s M=%.4f A/M=%.4f  e0=%.8f (%s)  sup_L D(L)/M = %.4f at L=%d ; pi A/(6M) = %.4f ; D(L)/M: %s' % (
            name, M, A / M, e0, src, max(D) / M, Ls[int(np.argmax(D))], np.pi * A / 6 / M, np.round(np.array(D) / M, 4)))
    json.dump(out, open(os.path.join(HERE, 'data', 'zeroT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    run()
