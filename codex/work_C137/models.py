"""C137 test bonds. Random bonds and the chiral bond are read from C135 models.json (same seed 1352026,
radially normalized to M = 1)."""
import json
import numpy as np
from common import xxz, tfim, C135, kron, PX, PY, PZ, I2
from e0 import e0_xxz, e0_tfim

F = np.eye(4)[[0, 2, 1, 3]]


def c135_models():
    m = json.load(open(C135 / 'models.json'))
    return {e['name']: np.array(e['h']['real']) + 1j * np.array(e['h']['imag']) for e in m}


def get(name, n=None):
    """Return dict(h, e0_exact or None, udirs, reflect)."""
    if name.startswith('XXZ_'):
        d = float(name[4:]); h = xxz(d)
        return dict(h=h, e0=e0_xxz(d), udirs=(), reflect=True)
    if name == 'Heisenberg':
        return dict(h=xxz(1), e0=e0_xxz(1), udirs=(), reflect=True)
    if name == 'XX':
        return dict(h=xxz(0), e0=e0_xxz(0), udirs=(), reflect=True)
    if name == 'TFIM_1':
        return dict(h=tfim(1), e0=e0_tfim(1), udirs=(2,), reflect=True)
    if name.startswith('TFIMc'):          # TFIMc<c>: g = 1 + c/n, depends on n
        c = float(name[5:]); g = 1 + c / n
        return dict(h=tfim(g), e0=e0_tfim(g), udirs=(2,), reflect=True, g=g)
    if name == 'chiral':
        h = c135_models()['chiral']
        # On open chains the pi/2 twist is removable: chiral = Heisenberg - 1/4 exactly, e0 = -ln 2.
        return dict(h=h, e0=e0_xxz(1) - 0.25, udirs=(2,), reflect=True)
    if name.startswith('random_'):
        h = c135_models()[name]
        real = np.max(np.abs(h.imag)) < 1e-14
        sym = np.linalg.norm(F @ h @ F - h) < 1e-12
        return dict(h=h, e0=None, udirs=(0, 2) if real else (0, 1, 2), reflect=bool(sym))
    raise KeyError(name)
