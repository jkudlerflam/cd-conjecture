"""Second-order logical chain of the folded-bus cell at the working point:
hopping matrix T, one-cell matrix mu, non-Gaussian residual, Hermitian part
of the phase-rotated T, bands, Fermi points, centres and resonance margins.
Usage: python3 logical_chain.py n [td ...].  Writes logical_n{n}.json."""
import sys, json
import numpy as np
from ycell import params, build_cell, code, second_order, logical, bands, resonance_margin

IR_PARAMS = dict(delta=0.0, Delta=0.0, Dx=4.0, DM=4.0, gp=0.05, gm=0.05,
                 phd={'0': 0.9, '+': 0.0, '-': 0.0})


def analyse(n, td, eps_frac=0.25, nk=40001, **over):
    P = params(td=td, **dict(IR_PARAMS, **over))
    cell = build_cell(n, P)
    E0, G, gap = code(cell)
    H2, k = second_order(cell, G, E0)
    # Hilbert-Schmidt projection onto quadratic operators: the quadratic part
    # referred to half filling (Hartree-consistent near half filling).  The
    # remainder is the non-Gaussian residual.  Reliable in double precision
    # for n <= 3 at these parameters.
    L = logical(H2, n)
    T, mu = L['T'], L['muA'] + L['muB']
    sc = np.abs(T).max()
    R = T * np.exp(-1j * np.angle(T[0, 0]))
    rs = np.linalg.eigvalsh((R + R.conj().T) / 2)
    # common logical chemical potential from a uniform data field, inside the
    # narrowest band: |eps| = eps_frac * 2 * min|eig(Herm part)|
    eps = eps_frac * 2 * np.abs(rs).min()
    bs = bands(T, mu + eps * np.eye(n), nk=nk)
    roots = [r for b in bs for r in b['roots']]
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    vmin = np.inf
    for r in roots:   # numerical Fermi velocity
        h = 1e-5
        e1 = np.linalg.eigvalsh(mu + eps * np.eye(n) + T * np.exp(1j * (r + h)) + T.conj().T * np.exp(-1j * (r + h)))
        e0 = np.linalg.eigvalsh(mu + eps * np.eye(n) + T * np.exp(1j * (r - h)) + T.conj().T * np.exp(-1j * (r - h)))
        vmin = min(vmin, np.min(np.abs(e1 - e0)[np.argsort(np.abs(e1 + e0))[:1]]) / (2 * h))
    d4, d2 = resonance_margin(roots) if len(roots) > 1 else (None, None)
    out = dict(n=n, td=td, E0=E0, cell_gap=gap, poles=k, T_scale=sc,
               T=[[str(x) for x in r] for r in T], mu=[[str(x) for x in r] for r in mu],
               mu_over_T=float(np.abs(mu).max() / sc), nonGaussian_over_T=float(L['res'] / sc),
               herm_part_eigs_over_T=list(rs / sc), eps_over_T=float(eps / sc),
               roots_per_band=[len(b['roots']) for b in bs],
               centres=[b['centre'] for b in bs], fermi_points=roots,
               delta4=d4, delta2=d2, vmin_over_T=float(vmin / sc))
    return out


if __name__ == '__main__':
    n = int(sys.argv[1])
    tds = [float(x) for x in sys.argv[2:]] or [0.0]
    recs = []
    for td in tds:
        o = analyse(n, td)
        print(json.dumps({k: o[k] for k in ['n', 'td', 'cell_gap', 'poles', 'T_scale', 'mu_over_T',
                                             'nonGaussian_over_T', 'herm_part_eigs_over_T', 'eps_over_T',
                                             'roots_per_band', 'centres', 'delta4', 'delta2', 'vmin_over_T']},
                         default=float), flush=True)
        recs.append(o)
    with open('logical_n%d.json' % n, 'w') as f:
        json.dump(recs, f, indent=1, default=float)
