"""Independent checks of an explicit Gaussian stability certificate.

No imports from the proposed construction. Run with Python and NumPy.
The proof is in GAUSSIAN_CERTIFICATE.md; these finite tests are not its proof.
"""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import json

import numpy as np


def eigh(a, b=None, eigvals_only=False):
    if b is not None:
        inv = matrix_power_spd(b, -.5)
        a = inv @ a @ inv
    if eigvals_only:
        return np.linalg.eigvalsh((a+a.T)/2)
    return np.linalg.eigh((a+a.T)/2)


def block_diag(a, b):
    return np.block([[a,np.zeros((len(a),len(b)))],
                     [np.zeros((len(b),len(a))),b]])


def matrix_power_spd(a, power):
    val, vec = eigh((a + a.T) / 2)
    assert val[0] > 0, val[0]
    return (vec * val**power) @ vec.T


def metric(g):
    """Equal-time dimension metric for [phi,theta] in canonical units."""
    n = len(g) // 2
    j = np.block([[np.zeros((n, n)), np.eye(n)],
                  [np.eye(n), np.zeros((n, n))]])
    root = matrix_power_spd(g, .5)
    inverse = matrix_power_spd(g, -.5)
    val, vec = eigh(root @ j @ root)
    return inverse @ ((vec * np.abs(val)) @ vec.T) @ inverse


def relative_metric_bounds(m, m0):
    val = eigh(m, m0, eigvals_only=True)
    return float(val[0]), float(val[-1])


def generic_matrix_cases():
    rng = np.random.default_rng(20261028)
    rows = []
    # Nonzero phi-theta blocks, unequal velocities, and many components.
    for n in [1, 2, 7, 22, 53]:
        for spread in [1., 10., 100.]:
            vel = np.geomspace(1. / spread, spread, n)
            g0 = np.diag(np.r_[vel, vel])
            root0 = np.diag(np.sqrt(np.r_[vel, vel]))
            assert np.max(np.abs(metric(g0) - np.eye(2*n))) < 1e-10
            for delta in [.05, .25, 5. / 13 - 1e-5, .60]:
                x = rng.normal(size=(2*n, 2*n))
                q, _ = np.linalg.qr(x)
                eig = np.linspace(-delta, delta, 2*n)
                r = (q * eig) @ q.T
                g = root0 @ (np.eye(2*n) + r) @ root0
                m = metric(g)
                lower, upper = relative_metric_bounds(m, np.eye(2*n))
                kappa = np.sqrt((1 - delta) / (1 + delta))
                assert lower >= kappa - 2e-7
                assert upper <= 1/kappa + 2e-7
                # The positive solution of M G M = J G J is unique.
                j = np.block([[np.zeros((n,n)), np.eye(n)],
                              [np.eye(n), np.zeros((n,n))]])
                residual = np.linalg.norm(m @ g @ m - j @ g @ j, 2)
                residual /= np.linalg.norm(g, 2)
                assert residual < 1e-7
                rows.append(dict(modes=n, velocity_ratio=spread**2,
                                 delta=delta, certified_lower=kappa,
                                 observed_lower=lower, observed_upper=upper,
                                 riccati_relative_residual=float(residual)))
    return rows


def sharpness_cases():
    out = []
    for delta in [.05, .25, 5. / 13, .60]:
        # This one-mode squeeze saturates both sides of the bound.
        g = np.diag([1 + delta, 1 - delta])
        kappa = np.sqrt((1 - delta)/(1 + delta))
        expected = np.diag([kappa, 1/kappa])
        error = float(np.max(np.abs(metric(g) - expected)))
        assert error < 1e-12
        out.append(dict(delta=delta, kappa=kappa,
                        six_field_lower_bound=3*kappa, error=error))
    return out


def drift_cases():
    out = []
    for aa, bb in [(1.,1.), (1.25,.75), (.2,7.)]:
        for drift_fraction in [-.99,-.4,0.,.7,.99]:
            cc = drift_fraction*np.sqrt(aa*bb)
            g = np.array([[aa,cc],[cc,bb]])
            expected = np.diag([np.sqrt(bb/aa),np.sqrt(aa/bb)])
            err = float(np.max(np.abs(metric(g)-expected)))
            assert err < 1e-10
            out.append(dict(A=aa,B=bb,C=cc,
                            right_velocity=float(cc+np.sqrt(aa*bb)),
                            left_velocity=float(cc-np.sqrt(aa*bb)),
                            covariance_drift_independence_error=err))
    return out


def physical_coordinate_cases():
    out = []
    for n in [1, 6, 21]:
        f = .25 + np.arange(1, n+1)/(12*(n+1))
        vel = np.r_[1., np.tan(np.pi*f)]
        for b in [4, 100, 10000]:
            a = b + .5
            t = np.eye(n+1)
            t[0,0] = a
            t[1:,0] = -f
            ti = np.linalg.inv(t)
            # phi_phys=Ti phi_sq, theta_phys=T^T theta_sq.
            s = block_diag(ti, t.T)
            si = block_diag(t, ti.T)
            j = np.block([[np.zeros((n+1,n+1)), np.eye(n+1)],
                          [np.eye(n+1), np.zeros((n+1,n+1))]])
            assert np.max(np.abs(s @ j @ s.T - j)) < 1e-10
            g0_sq = np.diag(np.r_[vel,vel])
            g0 = si.T @ g0_sq @ si
            m0 = s @ s.T
            delta = .25
            # A sharp squeeze in canonical coordinates, transformed to physical.
            g_sq = np.diag(np.r_[(1+delta)*vel, (1-delta)*vel])
            g = si.T @ g_sq @ si
            kappa = np.sqrt((1-delta)/(1+delta))
            m = s @ np.diag(np.r_[np.full(n+1,kappa),
                                  np.full(n+1,1/kappa)]) @ s.T
            # Solve only a well-conditioned congruent generalized eigenproblem.
            lo, hi = relative_metric_bounds(si @ m @ si.T,
                                            si @ m0 @ si.T)
            assert abs(lo-kappa) < 1e-12 and abs(hi-1/kappa) < 1e-12
            residual = np.linalg.norm(m @ g @ m - j @ g @ j, 2)
            residual /= np.linalg.norm(g, 2)
            assert residual < 1e-10
            # A physically Euclidean-small interaction would not certify this.
            out.append(dict(payload_modes=n, b=b,
                            marker_K=float(m0[0,0]),
                            stiffness_min=float(eigh(g0, eigvals_only=True)[0]),
                            certified_metric_lower=lo,
                            riccati_relative_residual=float(residual)))
    return out


def resonance_cases():
    """Exact rational enumeration of every girth-six low-dimensional class.

    Species neutrality leaves z=0 through four fields, hence integer u and
    Delta_pay=|u|^2. We enumerate all nonzero u with |u|^2<=2 and h=0,+/-1.
    Fractions express momenta in units of 2*pi. Every such class must detune.
    """
    out = []
    for n in [1, 6, 21, 52, 186]:
        f = [Fraction(1,4) + Fraction(a,12*(n+1)) for a in range(1,n+1)]
        payload = []
        for a in range(n):
            for sign in [-1, 1]:
                payload.append((1, sign*f[a]))
        for a, b in combinations(range(n), 2):
            for sa, sb in product([-1,1],repeat=2):
                payload.append((2, sa*f[a]+sb*f[b]))
        min_pay = min(min(q % 1, (-q) % 1) for dim,q in payload)
        assert min_pay >= Fraction(1,12*(n+1))
        mixed = []
        # Endpoints suffice for this interval because the two disjoint ranges
        # cannot cross an integer; also check the midpoint exactly.
        for nu in [Fraction(5,12), Fraction(1,2), Fraction(7,12)]:
            for dim,q in payload:
                if dim != 1:
                    continue
                for h in [-1,1]:
                    mixed.append(min((q-h*nu)%1, (-q+h*nu)%1))
        min_mixed = min(mixed)
        assert min_mixed >= Fraction(1,12)
        # Explicit large-b conversion D+S=h(b+nu) for pure and mixed resonances.
        conversion_count = 0
        for b in [4,100,10000]:
            a = Fraction(b) + Fraction(1,2)
            for h in [-3,-2,-1,0,1,2,3]:
                for integer_shift in [-2,-1,0,1,2]:
                    # A chosen S resonant at nu=1/2, irrespective of its origin.
                    S = Fraction(h,2) + integer_shift
                    D = h*b-integer_shift
                    assert Fraction(D+S)**2/a**2 == h*h
                    conversion_count += 1
        out.append(dict(payload_modes=n, low_payload_vertices=len(payload),
                        min_payload_mismatch_over_2pi=str(min_pay),
                        min_mixed_mismatch_over_2pi=str(min_mixed),
                        exact_large_b_harmonic_checks=conversion_count))
    return out


def main():
    result = dict(status="independent finite checks; analytic proof is separate",
                  generic_matrices=generic_matrix_cases(),
                  sharpness=sharpness_cases(),
                  one_mode_drift=drift_cases(),
                  physical_coordinates=physical_coordinate_cases(),
                  resonances=resonance_cases())
    target = Path(__file__).with_name("gaussian_bound_results.json")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(dict(generic_matrices=len(result["generic_matrices"]),
                         sharpness_cases=len(result["sharpness"]),
                         one_mode_drift_cases=len(result["one_mode_drift"]),
                         physical_coordinate_cases=len(result["physical_coordinates"]),
                         resonance_sizes=[r["payload_modes"] for r in result["resonances"]],
                         output=str(target)),indent=2))


if __name__ == "__main__":
    main()
