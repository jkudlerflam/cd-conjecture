"""Small controls for C51. Not a simulation of the microscopic wrapper."""
import json
import math
from pathlib import Path
import numpy as np


def diagonalization():
    rng = np.random.default_rng(51)
    records = []
    for n in (1, 2, 3):
        s = np.diag([1.0] * n + [-1.0] * n)
        a = rng.normal(size=(2*n, 2*n))
        v = a.T @ a + 0.3 * np.eye(2*n)
        eig, basis = np.linalg.eigh(v)
        root = (basis * np.sqrt(eig)) @ basis.T
        invroot = (basis * (1 / np.sqrt(eig))) @ basis.T
        d, o = np.linalg.eigh(root @ s @ root)
        order = np.r_[np.flatnonzero(d > 0), np.flatnonzero(d < 0)]
        d, o = d[order], o[:, order]
        assert np.count_nonzero(d > 0) == n
        u = np.abs(d)
        r = invroot @ o @ np.diag(np.sqrt(u))
        errors = [np.max(np.abs(r.T @ s @ r - s)),
                  np.max(np.abs(r.T @ v @ r - np.diag(u))),
                  np.max(np.abs(s @ v @ r - r @ np.diag(d)))]
        assert max(errors) < 2e-11
        records.append(dict(N=n, speeds=u.tolist(), max_errors=errors))
    # Exact diagonal input also checks that equality and extreme speed ratios
    # are permitted by the algebra, independently of mixed-matrix conditioning.
    u = np.array([1e-7, 1.0, 1.0, 2e-7, 0.6, 9.0])
    s = np.diag([1., 1., 1., -1., -1., -1.])
    assert np.allclose(s @ np.diag(u), np.diag(np.diag(s)*u))
    return records


def fock_states(max_level):
    states = []
    def rec(mode, remaining, entries):
        if mode > max_level:
            if remaining == 0:
                states.append(tuple(entries))
            return
        for count in range(remaining // mode + 1):
            rec(mode + 1, remaining - count * mode, entries + [count])
    for level in range(max_level + 1):
        rec(1, level, [])
    return states


def virasoro():
    cutoff = 10
    states = fock_states(cutoff)
    index = {state: j for j, state in enumerate(states)}
    levels = np.array([sum((k+1)*a for k, a in enumerate(s)) for s in states])
    test = np.flatnonzero(levels <= 3)
    def act(state, signed_mode):
        entries = list(state)
        mode = abs(signed_mode)
        count = entries[mode-1]
        if signed_mode > 0:
            if count == 0:
                return None, 0.0
            entries[mode-1] -= 1
            factor = math.sqrt(mode * count)
        else:
            entries[mode-1] += 1
            factor = math.sqrt(mode * (count+1))
        out = tuple(entries)
        if out not in index:
            return None, 0.0
        return out, factor
    reports = []
    eye = np.eye(len(states))[:, test]
    for charge in (0.0, 0.37, math.sqrt(2)):
        operators = {0: np.diag(levels + charge*charge/2)}
        # Modes up to 6 needed on the right of commutators with |m|,|n|<=3.
        for m in range(1, 7):
            matrix = np.zeros((len(states), len(states)))
            for j, state in enumerate(states):
                out, factor = act(state, m)
                if out is not None:
                    matrix[index[out], j] += charge*factor
                for k in range(1, m):
                    mid, f1 = act(state, k)
                    if mid is not None:
                        out, f2 = act(mid, m-k)
                        if out is not None:
                            matrix[index[out], j] += 0.5*f1*f2
                for k in range(1, cutoff-m+1):
                    mid, f1 = act(state, m+k)
                    if mid is not None:
                        out, f2 = act(mid, -k)
                        if out is not None:
                            matrix[index[out], j] += f1*f2
            operators[m] = matrix
            operators[-m] = matrix.T
        worst = 0.0
        for m in range(-3, 4):
            for n in range(-3, 4):
                lhs = operators[m] @ operators[n][:, test] - operators[n] @ operators[m][:, test]
                rhs = (m-n)*operators[m+n][:, test]
                if m+n == 0:
                    rhs = rhs + (m**3-m)/12*eye
                worst = max(worst, float(np.max(np.abs(lhs-rhs))))
        assert worst < 1e-11
        reports.append(dict(charge=charge, maximum_commutator_error=worst))
    # Partition counts are an independent integer recurrence.
    p = [1] + [0]*cutoff
    for mode in range(1, cutoff+1):
        for level in range(mode, cutoff+1):
            p[level] += p[level-mode]
    assert [int(np.count_nonzero(levels == level)) for level in range(cutoff+1)] == p
    return dict(dimension=len(states), tested_levels_through=3, reports=reports,
                partition_numbers=p)


def shifted_sectors():
    a = np.array([0.7, 1.3])
    # The chosen diagonal quadratic form is only an abstract compact-sector
    # fixture, not a claim about the wrapper compactification lattice.
    speeds = np.array([0.8, 1.6])
    casimir = -math.pi/12*speeds.sum()
    raw, corrected = [], []
    for length in range(128, 513):
        alpha = np.mod(length*np.array([math.sqrt(2), math.sqrt(3)]), 1)
        qmin = np.sum(a*(alpha-np.rint(alpha))**2)
        coefficient = 2*math.pi*qmin + casimir
        raw.append(coefficient)
        corrected.append(coefficient - 2*math.pi*qmin)
    error = float(np.max(np.abs(np.array(corrected)-casimir)))
    assert error < 1e-14 and max(raw)-min(raw) > 0.5
    beta = 2.3
    theta_records = []
    for length in (256, 1024, 4096):
        t = 2*math.pi*beta/length
        alpha = np.mod(length*np.array([math.sqrt(2), math.sqrt(3)]), 1)
        cutoff = math.ceil(math.sqrt(48/(t*a.min())))+3
        integers = np.arange(-cutoff, cutoff+1)
        logtheta = sum(math.log(float(np.exp(-t*ai*(integers+shift)**2).sum()))
                       for ai, shift in zip(a, alpha))
        lower = -t*a.sum()
        upper = sum(math.log(2+math.sqrt(math.pi/(t*ai))) for ai in a)
        assert lower <= logtheta <= upper
        theta_records.append(dict(L=length, log_theta=logtheta,
                                  pressure_magnitude=logtheta/(beta*length)))
    assert theta_records[-1]['pressure_magnitude'] < theta_records[0]['pressure_magnitude']
    return dict(raw_casimir_coefficient_range=[min(raw), max(raw)],
                corrected_coefficient=float(casimir), maximum_corrected_error=error,
                theta_records=theta_records)


def shell_bounds():
    omega = 0.4
    reports = []
    for p in (2.25, 2.4, 4.0):
        rows = []
        for k in (6, 10, 14, 18):
            # This computes the scalar error divided by T^2 with T=2^-k,
            # directly, to avoid very small floating-point quantities.
            h = np.arange(k+80, dtype=float)
            weight = np.exp2(np.where(h >= k, -2*(h-k), (p-2)*(h-k)))
            relative = float(np.sum(weight * 2.0**(-omega*h)))
            predicted = 2.0**(-k*min(omega, p-2))
            if abs(p-(2+omega)) < 1e-12:
                predicted *= 1 + k*math.log(2)
            rows.append(dict(k=k, error_over_T_squared=relative,
                             ratio_to_bound=relative/predicted))
        assert rows[-1]['error_over_T_squared'] < rows[0]['error_over_T_squared']
        assert max(row['ratio_to_bound'] for row in rows) < 20
        reports.append(dict(p=p, omega=omega, rows=rows))
    slow = []
    for k in (8, 32, 128, 512):
        h = np.arange(k+100, dtype=float)
        weight = np.exp2(np.where(h >= k, -2*(h-k), 2*(h-k)))
        relative = float(np.sum(weight/(h+1)))
        slow.append(dict(k=k, error_over_T_squared=relative))
    assert slow[-1]['error_over_T_squared'] < 0.03*slow[0]['error_over_T_squared']
    return dict(power_law=reports, slow_convergence=slow,
                p_equals_2_fixed_shell_error_over_T_squared=1.0)


def main():
    result = dict(diagonalization=diagonalization(), virasoro=virasoro(),
                  shifted_sectors=shifted_sectors(), shell_bounds=shell_bounds())
    Path(__file__).with_name('checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print('All C51 controls passed.')
    print('Maximum Virasoro error:', max(r['maximum_commutator_error'] for r in result['virasoro']['reports']))
    print('Maximum mixed canonical error:', max(max(r['max_errors']) for r in result['diagonalization']))
    print('Raw moving-sector coefficient range:', result['shifted_sectors']['raw_casimir_coefficient_range'])
    print('Corrected moving-sector coefficient error:', result['shifted_sectors']['maximum_corrected_error'])
    print('Saved checks.json. These are algebra and asymptotic controls, not microscopic central-charge data.')


if __name__ == '__main__':
    main()
