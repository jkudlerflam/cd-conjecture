"""C111 small checks. NumPy/SciPy only; no network or large diagonalizations."""
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])
I = np.eye(2)
XX = (np.kron(X, X) + np.kron(Y, Y)).real / 4


def spectator(delta, a=.75, lam=.1):
    h = np.zeros((16, 16))
    aa = [4*i+j for i in range(2) for j in range(2)]
    bb = [4*i+j for i in range(2, 4) for j in range(2, 4)]
    ha = (-a*np.eye(4) + delta/2*(np.kron(X, I)+np.kron(I, X))
          + lam*(np.kron(Z, I)-np.kron(I, Z)))
    h[np.ix_(aa, aa)] = ha
    for i in range(2):
        for j in range(2):
            h[bb[2*j+i], bb[2*i+j]] = 1
    for i in range(4):
        for j in range(4):
            if (i < 2) != (j < 2):
                h[4*i+j, 4*i+j] = 1
    return h


def chain(h, d, w, periodic=True):
    n = len(w) if periodic else len(w)+1
    dim = d**n
    assert dim <= 729, 'Local verification is restricted to small matrices.'
    codes = np.arange(dim)
    powers = d**np.arange(n)
    digits = (codes[:, None]//powers) % d
    rows, cols, data = [], [], []
    for x, weight in enumerate(w):
        y = (x+1) % n
        old = d*digits[:, x]+digits[:, y]
        for target in range(d*d):
            vals = weight*h[target, old]
            mask = np.abs(vals) > 1e-15
            to = codes + (target//d-digits[:, x])*powers[x]
            to += (target % d-digits[:, y])*powers[y]
            rows.extend(to[mask])
            cols.extend(codes[mask])
            data.extend(vals[mask])
    return coo_matrix((data, (rows, cols)), shape=(dim, dim)).tocsr()


def ground(mat, k=1):
    v0 = np.cos(np.sqrt(2)*np.arange(mat.shape[0]))
    if np.iscomplexobj(mat.data):
        v0 = v0 + 1j*np.sin(np.sqrt(3)*np.arange(mat.shape[0]))
    ev, vec = eigsh(mat, k=k, which='SA', v0=v0, tol=1e-12)
    order = np.argsort(ev)
    ev, vec = ev[order], vec[:, order]
    residual = max(np.linalg.norm(mat @ vec[:, j]-ev[j]*vec[:, j]) for j in range(k))
    return ev, float(residual)


def spin_deficit(w, delta, lam=.1, periodic=True):
    if periodic:
        right, left = w, np.roll(w, 1)
    else:
        right, left = np.r_[w, 0.], np.r_[0., w]
    u = (right+left)/2
    v = right-left
    # Rationalized expression avoids subtracting almost equal numbers.
    return float(np.sum(lam**2*v*v/(np.hypot(delta*u, lam*v)+delta*u)))


def qform(w, periodic=True):
    return float(np.sum((np.roll(w, -1)-w)**2/w) if periodic
                 else np.sum(np.diff(w)**2/w[:-1]))


def algebra_checks():
    rows = []
    for delta in [.05, .01, .001]:
        h = spectator(delta)
        spectrum = np.linalg.eigvalsh(h)
        assert abs(spectrum[0]+1) < 1e-12 and abs(spectrum[-1]-1) < 1e-12
        for n in [3, 4]:
            ev, residual = ground(chain(h, 4, np.ones(n)), k=2)
            expected = n*(-.75-delta)
            assert abs(ev[0]-expected) < 2e-11
            assert abs(ev[1]-ev[0]-2*delta) < 2e-11
            w = 1+.3*np.cos(2*np.pi*np.arange(n)/n)
            weighted, rw = ground(chain(h, 4, w))
            all_a = (-.75-delta)*w.sum()-spin_deficit(w, delta)
            assert weighted[0] <= all_a+2e-11
            rows.append(dict(delta=delta, n=n, energy=float(ev[0]),
                             expected=expected, gap=float(ev[1]-ev[0]),
                             raw_halfspread=float((spectrum[-1]-spectrum[0])/2),
                             weighted_energy=float(weighted[0]), all_a_trial=all_a,
                             residual=max(residual, rw)))
    # P12+P23 >= -I on the two-state B subspace.
    swap = np.eye(4).reshape(2, 2, 2, 2).transpose(1, 0, 2, 3).reshape(4, 4)
    pair = np.kron(swap, I)+np.kron(I, swap)
    assert abs(np.linalg.eigvalsh(pair)[0]+1) < 1e-12
    # Pair ordering here is (XX left, XX right, spectator left, spectator right).
    full = np.kron(XX, np.eye(16))+np.kron(np.eye(4), spectator(.01))
    es = np.linalg.eigvalsh(full)
    assert abs(es[0]+1.5) < 1e-12 and abs(es[-1]-1.5) < 1e-12
    return rows


def smooth_checks():
    rows = []
    for delta in [.05, .01, .001, .0001]:
        for n in [256, 1024, 4096, 16384]:
            w = 1+.3*np.cos(2*np.pi*np.arange(n)/n)
            deficit = spin_deficit(w, delta)
            q = qform(w)
            rows.append(dict(delta=delta, n=n, shape='positive_periodic',
                             deficit=deficit, Q=q, ratio=deficit/q,
                             predicted_ratio=.01/(2*delta)))
        previous = None
        for n in [256, 1024, 4096, 16384, 65536]:
            x = np.arange(1, n, dtype=float)
            w = x*(n-x)/n
            deficit = spin_deficit(w, delta, periodic=False)
            q = qform(w, periodic=False)
            slope = None if previous is None else (deficit-previous[1])/np.log(n/previous[0])
            rows.append(dict(delta=delta, n=n, shape='quadratic_open',
                             deficit=deficit, Q=q, log_slope=slope,
                             predicted_log_slope=.01/delta))
            previous = n, deficit
    return rows


def model_checks():
    models = [('XX', 2, XX),
              ('Heisenberg', 2, XX+np.kron(Z, Z)/4),
              ('XXZ_half', 2, XX+np.kron(Z, Z)/8)]
    swap3 = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            swap3[3*j+i, 3*i+j] = 1
    models.append(('SU3', 3, swap3))
    rng = np.random.default_rng(111)
    for k in range(4):
        b = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        h = (b+b.conj().T)/2
        ev = np.linalg.eigvalsh(h)
        h = (h-(ev[-1]+ev[0])/2*np.eye(4))/((ev[-1]-ev[0])/2)
        models.append((f'random_{k}', 2, h))
    rows = []
    for name, d, h in models:
        for n in ([6, 8] if d == 2 else [6]):
            # Include a second harmonic, since the first harmonic alone can
            # annihilate the uniform vacuum and give a vanishing response.
            angle = 2*np.pi*np.arange(n)/n
            w = 1+.25*np.cos(angle)+.15*np.cos(2*angle)
            uniform, ru = ground(chain(h, d, np.ones(n)))
            weighted, rw = ground(chain(h, d, w))
            # The finite-ring reference removes the uniform Casimir offset.
            defect = uniform[0]*w.sum()/n-weighted[0]
            q = qform(w)
            ev = np.linalg.eigvalsh(h)
            mraw = (ev[-1]-ev[0])/2
            assert defect >= -2e-10
            rows.append(dict(model=name, d=d, n=n, dimension=d**n,
                             finite_reference_deficit=float(defect), Q=q,
                             M_raw=float(mraw), ratio_raw=float(defect/(mraw*q)),
                             residual=max(ru, rw)))
    return rows


if __name__ == '__main__':
    start = time.perf_counter()
    result = dict(algebra=algebra_checks(), smooth=smooth_checks(), models=model_checks())
    result['elapsed_seconds'] = time.perf_counter()-start
    Path('checks_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print('algebra checks passed:', len(result['algebra']))
    print('largest eigenvector residual:', max(r['residual'] for r in result['algebra']+result['models']))
    for delta in [.05, .01, .001, .0001]:
        rows = [r for r in result['smooth'] if r['delta'] == delta]
        print('smooth', rows[3], rows[-1])
    for row in result['models']:
        print('model', row)
    print('elapsed_seconds:', result['elapsed_seconds'])
