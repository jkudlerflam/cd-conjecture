"""Registry of d = 2 chains with a positive checkerboard transfer (bonds in the reflection-positive frame).
name: (bond h', charge mode, M (gauge optimized), infrared A = sum c v or 0 if gapped, label)."""
import numpy as np
from cflib import bond_xxz_rot, bond_positive

MODELS = {
    'XX':      (bond_xxz_rot(0.0), 'u1', 0.5, 1.0, 'XX, Delta = 0'),
    'XXZ05':   (bond_xxz_rot(0.5), 'u1', 0.5, 3 * np.sqrt(3) / 4, 'XXZ, Delta = 0.5'),
    'HEIS':    (bond_xxz_rot(1.0), 'u1', 0.5, np.pi / 2, 'antiferromagnetic Heisenberg'),
    'XXZ2':    (bond_xxz_rot(2.0), 'u1', 0.75, 0.0, 'XXZ, Delta = 2 (gapped Neel)'),
    'TFIM':    (bond_positive(Jx=4.0, gz=2.0), 'z2', np.sqrt(2), 1.0, 'critical transverse-field Ising, H = -sum sxsx - sum sz'),
    'XYI':     (bond_positive(Jx=1.0, Jy=0.5, gz=0.25), 'z2', None, 0.375, 'anisotropic XY on its Ising line (gamma = 3)'),
    'XYZ':     (bond_positive(1.0, 0.5, 0.3), 'z2', 0.375, 0.0, 'XYZ (1, 0.5, 0.3), gapped'),
    'XXZSTAG': (bond_positive(1.0, 1.0, 0.5, 0.0, 0.3), 'u1', None, 0.0, 'XXZ Delta = 0.5 with staggered field 0.3 (gapped)'),
}
LEGACY = {'XX': 0.0, 'XXZ05': 0.5, 'HEIS': 1.0, 'XXZ2': 2.0}


def get(name):
    h, mode, M, A, label = MODELS[name]
    if M is None:
        from cflib import gauge_M
        M = gauge_M(h)
    return h, mode, M, A, label
