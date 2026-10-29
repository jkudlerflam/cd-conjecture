"""Finite-ring diagnostic of the explicit second-order penalty normal form."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
import numpy as np
from scipy.linalg import expm
from check_vacancy import build


def main():
    t, w, pen, _, states = build(2, 1, 1, True)
    B = (t+w).toarray()
    diff = pen[:, None]-pen[None, :]
    mask = diff != 0
    B0, off = B*(~mask), B*mask
    out = []
    for K in (20., 40., 80., 160.):
        s1 = np.zeros_like(B)
        s1[mask] = B[mask]/(K*diff[mask])
        comm = lambda a, b: a@b-b@a
        r2 = comm(s1, B0)+comm(s1, off)/2
        s2 = np.zeros_like(B)
        s2[mask] = r2[mask]/(K*diff[mask])
        S = s1+s2
        assert np.max(abs(S+S.T)) < 1e-14
        U = expm(S)
        H = K*np.diag(pen)+B
        D = K*np.diag(pen)+B0+r2*(~mask)
        err = np.linalg.norm(U@H@U.T-D, 2)
        rec = dict(K=K, dim=len(states), remainder_norm=float(err),
                   remainder_times_K_squared=float(err*K*K),
                   unitarity_error=float(np.linalg.norm(U@U.T-np.eye(len(states)), 2)))
        out.append(rec)
        print(json.dumps(rec))
    with open('sw_checks.json', 'w') as f:
        json.dump(out, f, indent=2)


if __name__ == '__main__':
    main()
