"""Check the archived target fit points and the conditional two-boson metric algebra."""

import json
from pathlib import Path
import re

import numpy as np


HERE = Path(__file__).resolve().parent
TABLES = HERE.parents[1] / "d2_xy_yx_chirality_imps_20261118/tables"
TAG = "tg_d0.18750_h1.78300_rot"


def target_rows(file, prefix):
    text = (TABLES / file).read_text().split("== " + TAG, 1)[1].split("==", 1)[0]
    return [s for s in text.splitlines() if s.strip().startswith(prefix)]


def number(row, name):
    return float(re.search(r"(?:^|\s)" + re.escape(name) + r"=\s*([+\-0-9.eE]+)", row)[1])


def fit(x, y):
    design = np.column_stack([np.ones(len(x)), x])
    coefficients, rss, _, _ = np.linalg.lstsq(design, y, rcond=None)
    se = np.sqrt((rss[0] / (len(y) - 2)) * np.linalg.inv(design.T @ design)[1, 1])
    return dict(slope=float(coefficients[1]), regression_se=float(se))


def estimators():
    infinite = [r for r in target_rows("target_ladders.txt", "chi=") if number(r, "chi") >= 96]
    assert [int(number(r, "chi")) for r in infinite] == [96, 128, 160, 200, 256]
    x = np.log([number(r, "xi") for r in infinite])
    entropy = fit(x / 6, [number(r, "S") for r in infinite])
    chirality = fit(x, -np.log([abs(number(r, "m")) for r in infinite]))

    finite = [r for r in target_rows("open_chain_ladders.txt", "L=")
              if number(r, "L") >= 96 and abs(number(r, "dSbar(chi)")) <= .02]
    lengths = [int(number(r, "L")) for r in finite]
    assert lengths == [96, 128, 192, 256, 384]
    coord = []
    for length in lengths:
        cuts = np.arange(3 * length // 8, 5 * length // 8 + 1)
        coord.append(float(np.log((2 * length / np.pi) * np.sin(np.pi * cuts / length)).mean() / 6))
    means = [number(r, "Sbar") for r in finite]
    return dict(iMPS_c=entropy, iMPS_kappa_effective_x=chirality,
                OBC_lengths=lengths, OBC_c=fit(coord, means), OBC_without_L96=fit(coord[1:], means[1:]))


def square_root(matrix):
    eigenvalues, vectors = np.linalg.eigh(matrix)
    assert eigenvalues.min() > 0
    return (vectors * np.sqrt(eigenvalues)) @ vectors.T


def gaussian_examples():
    # Pick an explicit canonical change of field coordinates, separately from
    # reconstructing G from H. These are examples, not fitted target velocities.
    row = next(r for r in (TABLES / "correlators_odd_even.txt").read_text().splitlines()
               if r.startswith(TAG + "_L384_"))
    p, eta_b, raw_K = (number(row, key) for key in ("pair3", "eta_b", "Kc"))
    d = 4 * eta_b - p
    zero, one = np.zeros((2, 2)), np.eye(2)
    eta = np.block([[zero, one], [one, zero]])
    discrepancies = []
    for t in [0., -.7, 1.5]:
        q = np.array([[p, p*t], [p*t, 1/d + p*t*t]])
        cholesky = np.linalg.cholesky(q)
        transform = np.block([[cholesky, zero], [zero, np.linalg.inv(cholesky).T]])
        canonical_G = transform @ transform.T
        inverse = np.linalg.inv(transform)
        for test_speeds in [(1., 1.), (.6, 1.7)]:
            h = inverse.T @ np.diag([*test_speeds, *test_speeds]) @ inverse
            hs = square_root(h)
            velocities, u = np.linalg.eigh(hs @ eta @ hs)
            hi = np.linalg.inv(hs)
            g = hi @ (u * abs(velocities)) @ u.T @ hi
            discrepancies.extend([np.max(abs(g - canonical_G)), np.max(abs(g @ eta @ g - eta)),
                                  np.max(abs(velocities - np.sort([*test_speeds, *(-np.array(test_speeds))])))])
            for r in range(-3, 4):
                for s in range(-3, 4):
                    ell = 2*np.sqrt(2) * np.array([r, -s, 0, 0])
                    dimension = ell @ g @ ell / 4
                    discrepancies.append(abs(dimension - 2 * (p*(r - t*s)**2 + s*s/d)))
            assert g[2, 2] + 1.e-12 >= 1/p and g[1, 1] + 1.e-12 >= 1/d
    assert max(discrepancies) < 1.e-10
    return dict(eta_b=eta_b, pair_exponent=p, d_if_matched=d, x_quartet_if_matched=2*p,
                lower_bound_on_generic_s_nonzero_x_if_matched=2/d,
                lower_bound_on_matched_density_K=1/p, raw_density_K=raw_K,
                product_test=raw_K*p, growth_factor_96_to_384_if_linearized=4**(2 - 2*p),
                max_error_on_algebra_examples=float(max(discrepancies)))


if __name__ == "__main__":
    print(json.dumps(dict(archived_fits=estimators(), metric=gaussian_examples()), indent=2))
