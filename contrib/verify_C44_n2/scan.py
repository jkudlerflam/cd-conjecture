"""Scan of core parameters for the n=2 projected bond (s/delta -> 0 limit).
Reports the band velocities of the quadratic part, the interaction norm and
the length interaction V, all relative to the fast payload hopping."""
import sys, json
sys.path.insert(0, '../verify_C44')
import numpy as np
from core46 import Params
from precise import kbar_bonds
from ferm import decompose


def bands(T, nk=4001):
    """h(k) = T e^{ik} + T^dag e^{-ik} on one cell, T = <left|h|right> block.
    Returns Fermi momenta and velocities at mu = 0."""
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.array([np.linalg.eigvalsh(T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k)) for k in ks])
    out = []
    for b in range(E.shape[1]):
        e = E[:, b]
        s = np.sign(e)
        idx = np.where(s != np.roll(s, -1))[0]
        for i in idx:
            j = (i + 1) % nk
            k0 = ks[i] - e[i] * (2 * np.pi / nk) / (e[j] - e[i])
            v = (e[j] - e[i]) / (2 * np.pi / nk)
            out.append((b, float(k0), float(v)))
    return E, out


def analyse(n=2, rho=0.05, **pk):
    P = Params(n, **pk)
    K = kbar_bonds(n, rho, P=P)
    K = {k: (v + v.conj().T) / 2 for k, v in K.items()}
    sc, h, rest = decompose(K['11'], n)
    T = h[:n, n:]
    E, fp = bands(T)
    c = {k: np.trace(v).real / 4 ** n for k, v in K.items()}
    V = c['22'] - c['21'] - c['12'] + c['11']
    mixed = max(np.abs((K[k] - c[k] * np.eye(4 ** n)) - (K['11'] - c['11'] * np.eye(4 ** n))).max() for k in K)
    onsite = np.abs(h[:n, :n]).max() + np.abs(h[n:, n:]).max()
    return dict(T=T, fermi=fp, rest=float(np.abs(rest).max()), V=float(V), mixed=float(mixed),
                onsite=float(onsite), tmax=float(np.abs(T).max()))


if __name__ == '__main__':
    for kw in json.loads(sys.argv[1]):
        rho = kw.pop('rho', 0.05)
        r = analyse(rho=rho, **kw)
        vs = sorted(abs(x[2]) for x in r['fermi'])
        print(kw, 'rho', rho, 'tmax %.3e' % r['tmax'], 'v/tmax', ['%.2e' % (v / r['tmax']) for v in vs],
              'rest/vslow %.2e' % (r['rest'] / vs[0]), 'V/vslow %.2e' % (r['V'] / vs[0]),
              'mixed/vslow %.2e' % (r['mixed'] / vs[0]), 'k', ['%.4f' % x[1] for x in r['fermi']], flush=True)
