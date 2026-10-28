"""Literal check of leading.py.  Builds the lambda^2 g^2 coefficient
K = -Z^dag R0 C R0 C R0 Z directly in the complete two-cell excited space of
the literal 23-letter cell (contrib/routeA_spatial/ycell.py, imported
read-only), on the complete code x code space (dimension 4^n), and compares
its hopping and onsite blocks with the structured high-precision code.
Also checks that K is exactly quadratic plus a scalar (Hilbert-Schmidt fit).
Usage: python3 -B literal_check.py [nmax]   (n = 4 takes a few minutes)
Writes literal_check.json."""
import sys, os, json, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'routeA_spatial'))
import ycell                                   # noqa: E402
from decimal import Decimal                    # noqa: E402
import hpmat                                   # noqa: E402
from leading import leading, alpha4            # noqa: E402

# generic point: every phase nonzero, loop strong, +- splitting on
PT = dict(t=1.0, tx=1.0, tsw=1.0, Delta=1.0, Dx=2.0, DM=2.5, delta=0.15,
          phs={'0': 0.2, '+': -0.3, '-': 0.25}, phx={'0': 0.7, '+': 0.1, '-': -0.4},
          phsw={'0': 0.35, '+': 0.0, '-': 0.5}, phd={'0': 0.9, '+': 1.3, '-': -0.6},
          td=0.3, bphi=1.1, eps=0.0, musel=0.0)


def literal_K(n, Pd):
    c0 = ycell.build_cell(n, ycell.params(**dict(Pd, gp=0.0, gm=0.0)))
    c1 = ycell.build_cell(n, ycell.params(**dict(Pd, gp=1.0, gm=1.0)))
    E0, G, gap = ycell.code(c0)
    ex = c0['exc']
    De = int(ex.sum())
    He = sp.csr_matrix(c0['H'][np.ix_(ex, ex)] - E0 * np.eye(De))
    Ce = sp.csr_matrix((c1['H'] - c0['H'])[np.ix_(ex, ex)])
    I = sp.identity(De, format='csr')
    A2 = (sp.kron(He, I) + sp.kron(I, He)).tocsc()
    C2 = (sp.kron(Ce, I) + sp.kron(I, Ce)).tocsr()
    lu = splu(A2)
    N = G.shape[1]
    Z = np.zeros((De * De, N * N), complex)
    for (fm, fx, c) in c0['terms']:
        a = (c0['FM'][fm] @ G)[ex]
        b = (c0['FX'][fx] @ G)[ex]
        Z += c * np.kron(a, b)
    X = lu.solve(Z)
    Y = lu.solve(C2 @ X)
    W = lu.solve(C2 @ Y)
    K = -Z.conj().T @ W
    return K, E0, gap, De


def blocks(K, n):
    m = 2 * n
    idx = [2 ** (m - 1 - j) for j in range(m)]
    T = np.array([[K[idx[x], idx[n + y]] for y in range(n)] for x in range(n)])
    kv = K[0, 0]
    muL = np.array([[K[idx[x], idx[y]] for y in range(n)] for x in range(n)]) - kv * np.eye(n)
    muR = np.array([[K[idx[n + x], idx[n + y]] for y in range(n)] for x in range(n)]) - kv * np.eye(n)
    return T, muL + muR, kv


def structured(n, Pd):
    hpmat.setprec(60)
    a = Decimal(Pd['Delta']) + 2 * hpmat.dcos(hpmat.dpi() / (n + 2))
    P = dict(a=a, Dx=Decimal(Pd['Dx']), DM=Decimal(Pd['DM']), delta=Decimal(Pd['delta']),
             t=Decimal(Pd['t']), tx=Decimal(Pd['tx']), tsw=Decimal(Pd['tsw']), td=Decimal(Pd['td']),
             phs={k: Decimal(v) for k, v in Pd['phs'].items()},
             phx={k: Decimal(v) for k, v in Pd['phx'].items()},
             phsw={k: Decimal(v) for k, v in Pd['phsw'].items()},
             phd={k: Decimal(v) for k, v in Pd['phd'].items()},
             bphi=Decimal(Pd['bphi']))
    L = leading(n, P)
    a4 = float(alpha4(n))   # leading.py removes alpha^4, the literal K keeps it
    return a4 * L['T'].tocomplex(), a4 * L['mu'].tocomplex(), a4 * float(L['Kvac'])


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    out = []
    for n in range(1, nmax + 1):
        for label, Pd in [('generic', PT),
                          ('C33-like', dict(PT, Delta=2.0, Dx=6.0, DM=6.0, delta=0.0,
                                            phs={'0': 0., '+': 0., '-': 0.}, phsw={'0': 0., '+': 0., '-': 0.},
                                            phx={'0': np.pi / (2 * n), '+': 0., '-': 0.},
                                            phd={'0': np.pi, '+': np.pi / 2, '-': np.pi / 2},
                                            td=0.05, bphi=np.pi / 2))]:
            t0 = time.time()
            K, E0, gap, De = literal_K(n, Pd)
            Tl, mul, kvl = blocks(K, n)
            e0, h, res = ycell.fit_quadratic(K, 2 * n)
            Ts, mus, kvs = structured(n, Pd)
            sc = np.abs(Tl).max()
            rec = dict(n=n, point=label, pair_dim=De * De, code_dim=4 ** n, cell_gap=float(gap),
                       T_scale=float(sc),
                       T_rel_err=float(np.abs(Ts - Tl).max() / sc),
                       mu_err_over_T=float(np.abs(mus - mul).max() / sc),
                       mu_over_T=float(np.abs(mul).max() / sc),
                       vac_rel_err=float(abs(kvs - kvl.real) / abs(kvl)),
                       quadratic_residual_over_T=float(res / sc),
                       hs_vs_oneparticle_T=float(np.abs(h[:n, n:] - Tl).max() / sc),
                       hs_vs_oneparticle_mu=float(np.abs(h[:n, :n] + h[n:, n:] - mul).max() / sc),
                       seconds=time.time() - t0)
            print(json.dumps(rec), flush=True)
            out.append(rec)
    json.dump(out, open(os.path.join(HERE, 'literal_check.json'), 'w'), indent=1)
