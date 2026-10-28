"""Small independent controls for C49's operator estimates. NumPy/SciPy only."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, expm

rng = np.random.default_rng(491028)


def hermitian(dim):
    x = rng.normal(size=(dim, dim)) + 1j*rng.normal(size=(dim, dim))
    x = (x+x.conj().T)/2
    return x/np.linalg.norm(x, 2)


def tensor(seq):
    out = np.ones((1, 1), complex)
    for x in seq:
        out = np.kron(out, x)
    return out


def embed(x, start, width, count, dim):
    return np.kron(np.kron(np.eye(dim**start), x),
                   np.eye(dim**(count-start-width)))


def norm(x):
    return float(np.linalg.norm(x, 2))


def finite_spectral_checks():
    records = []
    local_p = np.diag([1., 1., 0.])
    local_a = np.diag([0., 0., 1.])
    for count in (2, 3, 4, 5):
        dim = 3**count
        amat = sum(embed(local_a, j, 1, count, 3) for j in range(count))
        nmat = amat.copy()
        pmat = tensor([local_p]*count)
        pidx = np.flatnonzero(np.diag(pmat).real > .5)
        qidx = np.flatnonzero(np.diag(pmat).real < .5)
        dmat = np.zeros((dim, dim), complex)
        omat = np.zeros_like(dmat)
        dnorm = np.zeros(count)
        onorm = np.zeros(count)
        # Diagonal terms can move defects within a bond's complementary block.
        for width in (1, 2):
            px = tensor([local_p]*width)
            qx = np.eye(3**width)-px
            for j in range(count-width+1):
                raw = hermitian(3**width)
                dx = .012*(px@raw@px + qx@raw@qx)
                raw = hermitian(3**width)
                ox = .004*(px@raw@qx + qx@raw@px)
                dmat += embed(dx, j, width, count, 3)
                omat += embed(ox, j, width, count, 3)
                dnorm[j:j+width] += norm(dx)
                onorm[j:j+width] += norm(ox)
        jnorm = float(max(dnorm))
        epsilon = float(max(onorm))
        delta = 1-2*jnorm
        hmat = amat+dmat+omat
        kmat = hmat[np.ix_(pidx, pidx)]
        kvals = eigh(kmat, eigvals_only=True)
        hvals, hvecs = eigh(hmat)
        energy = float(kvals[0])
        bmat = hmat[np.ix_(qidx, qidx)]
        bmin = float(eigh(bmat, eigvals_only=True, subset_by_index=[0, 0])[0])
        leakage = norm(omat[np.ix_(qidx, pidx)])
        assert leakage**2 <= count*epsilon**2+1e-13
        assert count*epsilon**2 <= delta**2/8
        assert bmin >= energy+delta/2-2*count*epsilon**2/delta-1e-12
        bound = 8*count*epsilon**2/delta
        low = np.flatnonzero(kvals <= energy+delta/8)
        shifts = hvals[low]-kvals[low]
        assert np.max(shifts) <= 1e-12
        assert np.min(shifts) >= -bound-1e-12
        gap_errors = np.abs((hvals[low]-hvals[0])-(kvals[low]-energy))
        assert np.max(gap_errors) <= bound+1e-12
        assert hvals[0] >= energy-count*epsilon**2/delta-1e-12
        # Matrix form version of the Young bound, not just random expectations.
        lower = omat+(delta/2)*nmat+(2*count*epsilon**2/delta)*np.eye(dim)
        form_min = float(eigh(lower, eigvals_only=True, subset_by_index=[0, 0])[0])
        assert form_min >= -1e-12
        max_mixing_ratio = 0.
        for k in low:
            vec = hvecs[:, k]
            ratio = np.linalg.norm(vec[qidx])/np.linalg.norm(vec[pidx])
            max_mixing_ratio = max(max_mixing_ratio, float(ratio))
            assert ratio <= 8*math.sqrt(count)*epsilon/delta+1e-12
        records.append(dict(cells=count, dimension=dim, low_levels=len(low),
                            epsilon=epsilon, delta=delta, leakage_norm=leakage,
                            leakage_bound=math.sqrt(count)*epsilon,
                            spectral_bound=bound, max_level_shift=float(max(abs(shifts))),
                            max_gap_error=float(max(gap_errors)), form_min=form_min,
                            max_mixing_ratio=max_mixing_ratio))
    return records


def one_conjugation_check():
    # A one-support instance tests the exact cancellation and convergent BCH bound.
    amat = np.diag([0., 0., 1., 1.7])
    pmat = np.diag([1., 1., 0., 0.])
    qmat = np.eye(4)-pmat
    raw = hermitian(4)
    dmat = 1e-4*(pmat@raw@pmat+qmat@raw@qmat)
    raw = hermitian(4)
    omat = 1e-5*(pmat@raw@qmat+qmat@raw@pmat)
    inv = np.diag([0., 0., 1., 1/1.7])
    smat = inv@omat@pmat-pmat@omat@inv
    comm_error = norm(smat@amat-amat@smat+omat)
    umat = expm(smat)
    transformed = umat@(amat+dmat+omat)@umat.conj().T
    tmat = transformed-amat-dmat
    a = 1.
    d_weighted = math.exp(a)*norm(dmat)
    o_weighted = math.exp(a)*norm(omat)
    lam = d_weighted+o_weighted
    loss = 64*math.e*lam
    assert loss < a
    zeta = 4*math.e*math.exp(a)*norm(smat)/loss
    actual = math.exp(a-loss)*norm(tmat)
    majorant = zeta/(1-zeta)*lam
    assert comm_error < 1e-14
    assert actual <= majorant+1e-14
    assert actual <= o_weighted/7+1e-14
    # Independently sum the BCH expression in (10).
    ad_d, ad_o = dmat.copy(), omat.copy()
    series = np.zeros_like(dmat)
    for k in range(1, 9):
        ad_d = smat@ad_d-ad_d@smat
        ad_o = smat@ad_o-ad_o@smat
        series += ad_d/math.factorial(k)+(k/math.factorial(k+1))*ad_o
    series_error = norm(tmat-series)
    assert series_error < 1e-14
    return dict(commutator_error=comm_error, series_error=series_error,
                weighted_remainder=actual, majorant=majorant,
                required_contraction_bound=o_weighted/7, locality_loss=loss)


def product_rotation_checks():
    delta, t = 1., .4
    h = np.array([[0., t], [t, delta]])
    values, vectors = eigh(h)
    e_minus = float(values[0])
    alpha = (1+delta/math.sqrt(delta**2+4*t*t))/2
    records = []
    for count in (1, 2, 4, 6, 8):
        htot = sum(embed(h, j, 1, count, 2) for j in range(count))
        psi = vectors[:, 0]
        for _ in range(count-1):
            psi = np.kron(psi, vectors[:, 0])
        p = alpha**count
        qpsi = psi.copy()
        qpsi[0] = 0.
        qpsi /= np.linalg.norm(qpsi)
        ground = count*e_minus
        trial_shift = float(np.real(qpsi.conj()@htot@qpsi)-ground)
        formula = -ground*p/(1-p)
        qmin = float(eigh(htot[1:, 1:], eigvals_only=True,
                          subset_by_index=[0, 0])[0])
        assert abs(trial_shift-formula) < 1e-12
        assert qmin-ground >= -1e-12
        assert qmin-ground <= formula+1e-12
        records.append(dict(cells=count, overlap_squared=p,
                            complementary_gap=qmin-ground, trial_upper_bound=formula,
                            formula_error=abs(trial_shift-formula)))
    return records


def schedule_check():
    a, delta, lam, epsilon = 1., 1., 1e-4, 1e-7
    loss = 64*math.e*lam/delta
    steps = math.floor(a/(2*loss))
    remainder = math.ldexp(epsilon, -steps)
    upper = 2*epsilon*math.exp(-math.log(2)*a*delta/(128*math.e*lam))
    assert a-steps*loss >= a/2
    assert remainder <= upper
    # Check the recurrence for a full range of initial o/Lambda ratios.
    for ratio in np.linspace(0., 1., 101):
        zeta = ratio/8
        majorant = zeta/(1-zeta)  # in units of Lambda
        assert majorant <= ratio/7+1e-15
        assert 1-ratio+2*majorant <= 1+1e-15
    return dict(steps=steps, final_weight=a-steps*loss,
                remainder=remainder, exponential_bound=upper)


def main():
    result = dict(status='all checks passed',
                  finite_spectral=finite_spectral_checks(),
                  conjugation=one_conjugation_check(),
                  product_rotation=product_rotation_checks(),
                  schedule=schedule_check())
    Path(__file__).with_name('checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
