"""Top eigenvalue u0 of U_sym = K_e^{1/2} K_o K_e^{1/2} at large Trotter number by Lanczos with
tensor-structured matvecs (no sparse matrices), giving the bulk crossed vacuum energy eps0 = -(1/2) log u0
= beta f_m(beta). Used to extrapolate eps0 to m = infinity (Richardson in delta^2).
Usage: python3 eps0_lanczos.py NAME M beta1,beta2,...   -> data/eps0_NAME_m{M}.json"""
import sys, os, json, time
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import LinearOperator, eigsh
from cflib import realign, psd_sqrt
from models import get

HERE = os.path.dirname(os.path.abspath(__file__))


def make_op(h, beta, m):
    d = 2
    G = realign(expm(-(beta / m) * h), d)
    Gh, w = psd_sqrt(G)
    N = 4 ** m

    def pairs(v, A):
        t = v.reshape((4,) * m)
        for k in range(m):
            t = np.tensordot(A, t, axes=([1], [k]))
            t = np.moveaxis(t, 0, k)
        return t.reshape(-1)

    def roll(v, s):
        t = v.reshape((2,) * (2 * m))
        return np.ascontiguousarray(np.moveaxis(t, list(range(2 * m)), [(i - s) % (2 * m) for i in range(2 * m)])).reshape(-1)

    def mv(v):
        v = pairs(v, Gh)
        v = roll(pairs(roll(v, 1), G), -1)   # K_o = P^T K_e P on the rolled pairs
        return pairs(v, Gh)
    return LinearOperator((N, N), matvec=mv, dtype=float), G


if __name__ == '__main__':
    name, m = sys.argv[1], int(sys.argv[2])
    betas = [float(b) for b in sys.argv[3].split(',')]
    h = get(name)[0]
    out = {}
    for b in betas:
        t0 = time.time()
        A, G = make_op(h, b, m)
        v0 = np.ones(4 ** m)
        u, _ = eigsh(A, k=1, which='LA', v0=v0, tol=1e-13)
        out['%g' % b] = dict(u0=float(u[0]), eps0=float(-0.5 * np.log(u[0])), seconds=time.time() - t0)
        print(name, m, b, out['%g' % b], flush=True)
    json.dump(out, open(os.path.join(HERE, 'data', 'eps0_%s_m%d.json' % (name, m)), 'w'), indent=1)
