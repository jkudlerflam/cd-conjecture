"""Calibration of est.estimate on points of known central charge (C155, Claude Opus agent, 2026-11-07).
Usage: python3 calibrate.py L1 L2
"""
import sys, json, time
import numpy as np
from est import estimate, bond_from_JB, PAULI, I2, SX, SY, SZ

k = np.kron


def nn(J, f=(0, 0, 0)):
    return bond_from_JB(np.asarray(J, float), np.asarray(f, float))


def r3(*terms):
    h = np.zeros((8, 8), complex)
    for c, a, b, d in terms:
        h += c * k(k(a, b), d)
    return h


Z3 = np.zeros
models = {
    # name: (bond, range, expected c or tag)
    'TFI g=1': (nn(-np.diag([1, 0, 0]), (0, 0, -1)), 2, 0.5),
    'TFI g=0.9': (nn(-np.diag([1, 0, 0]), (0, 0, -0.9)), 2, 'gapped xi~10'),
    'TFI g=1.5': (nn(-np.diag([1, 0, 0]), (0, 0, -1.5)), 2, 'gapped'),
    'XY gamma=0.5 h=1': (nn(-np.diag([0.75, 0.25, 0]), (0, 0, -1)), 2, 0.5),
    'XX': (nn(np.diag([1, 1, 0])), 2, 1.0),
    'XX + DMz 0.7': (nn(np.diag([1, 1, 0]) + np.array([[0, .7, 0], [-.7, 0, 0], [0, 0, 0]])), 2, 1.0),
    'XX h=1 (k_F off pi/2)': (nn(np.diag([1, 1, 0]), (0, 0, 1)), 2, 1.0),
    'Heisenberg AF': (nn(np.eye(3)), 2, 1.0),
    'XXZ D=0.5': (nn(np.diag([1, 1, .5])), 2, 1.0),
    'XXZ D=-0.5': (nn(np.diag([1, 1, -.5])), 2, 1.0),
    'XXZ D=-0.9': (nn(np.diag([1, 1, -.9])), 2, 1.0),
    'XXZ D=2': (nn(np.diag([1, 1, 2])), 2, 'gapped'),
    'XXZ D=1.3': (nn(np.diag([1, 1, 1.3])), 2, 'gapped, xi large'),
    'XYZ 1,0.6,0.3': (nn(np.diag([1, .6, .3])), 2, 'gapped'),
    'XZ compass': (nn(np.diag([1, 0, 0]) + 0 * np.eye(3) + np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]) ), 2, 'classical Ising, gapped'),
    'FM Heisenberg': (nn(-np.eye(3)), 2, 'z=2, degenerate'),
    'XX saturation h=2': (nn(np.diag([1, 1, 0]), (0, 0, 2)), 2, 'z=2'),
    'XX h=2.2 polarized': (nn(np.diag([1, 1, 0]), (0, 0, 2.2)), 2, 'gapped product'),
    'FM XXZ D=-1 h=0.3': (nn(np.diag([1, 1, -1]), (0, 0, 0.3)), 2, 'z=2 boundary'),
    'XXZ Ising-AF h at onset': (nn(np.diag([0.2, 0.2, 1]), (0, 0, 0.0)), 2, 'gapped Neel'),
    # range 3 controls
    '2x TFI (NNN)': (r3((-1, SX, I2, SX), (-0.5, SZ, I2, I2), (-0.5, I2, I2, SZ)), 3, 1.0),
    '2x XX (NNN)': (r3((1, SX, I2, SX), (1, SY, I2, SY)), 3, 2.0),
    'XX + 2x XX NNN': (r3((1, SX, I2, SX), (1, SY, I2, SY), (.5, SX, SX, I2), (.5, SY, SY, I2), (.5, I2, SX, SX), (.5, I2, SY, SY)), 3, 'J1-J2 XX, c=1 or 2'),
    'cluster-Ising l=1': (r3((-1, SX, SZ, SX), (.5, SY, SY, I2), (.5, I2, SY, SY)), 3, 1.5),
    'XX + 2x TFI (NNN)': (r3((1, SX, I2, SX), (1, SY, I2, SY), (-1, SX, SX, I2)), 3, 'mixed'),
}

if __name__ == '__main__':
    Ls = tuple(int(a) for a in sys.argv[1:3]) if len(sys.argv) > 2 else (12, 16)
    res = {}
    for name, (h, r, exp) in models.items():
        t = time.time()
        e = estimate(h, Ls, r=r)
        res[name] = dict(expected=exp, **e)
        print(f"{name:26s} exp={str(exp):22s} c_S={e['c_S']:7.3f} c_fit={e['c_fit']:7.3f} c_cas={e['c_cas']:7.3f} "
              f"c_obj={e['c_obj']:7.3f} v={e['v']:.3f} m0={e['m0']} Lgap={[round(g, 3) for g in e['Lgap']]} ({time.time() - t:.1f}s)",
              flush=True)
    json.dump(res, open(f'calib_{Ls[0]}_{Ls[1]}.json', 'w'), indent=1)
