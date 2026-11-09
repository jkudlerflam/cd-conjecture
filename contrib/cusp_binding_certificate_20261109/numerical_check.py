"""Supplementary banded-matrix checks; all momenta are in the real gauge."""
import json
import platform
import numpy as np
import scipy
from scipy.linalg import eig_banded


def lowest(cutoff, u=1., v=0., g=2/3):
    h = np.zeros((3, cutoff))
    h[0, :2] = [-2+v, -g]
    h[1, :-1] = 2*u
    h[2, :-2] = 2
    e, psi = eig_banded(h, lower=True, select='i', select_range=(0, 1))
    return e, psi[:, 0]


def main():
    convergence = []
    for cutoff in [32, 64, 128, 256, 512, 1024]:
        e, psi = lowest(cutoff)
        convergence.append({
            "cutoff": cutoff, "eigenvalues": e.tolist(),
            "binding": float(-4.5-e[0]),
            "mean_separation": float(np.dot(np.arange(2, cutoff+2), psi**2)),
            "last_ten_weight": float(np.sum(psi[-10:]**2))
        })
    assert abs(convergence[-1]['eigenvalues'][0]+4.519526534666398) < 1e-11
    transition = []
    for u in [0.5, 1., 2., 3.5]:
        vc = u*u/4
        eth = -4-u*u/2
        for factor in [0., 0.8, 1., 1.2]:
            e, _ = lowest(1024, u=u, v=factor*vc)
            if factor >= 1:
                assert e[0] >= eth-1e-12
            if factor == 0:
                assert e[0] < eth
            transition.append({"u": u, "v": factor*vc, "v_critical": vc,
                               "energy_minus_threshold": float(e[0]-eth)})
    print(json.dumps({"python": platform.python_version(), "numpy": np.__version__,
                      "scipy": scipy.__version__, "convergence": convergence,
                      "threshold_checks": transition,
                      "scope": "finite relative-coordinate cutoffs, supplementary to exact proof"}, indent=2))


if __name__ == '__main__':
    main()
