#!/usr/bin/env python3
"""Finite-volume nested SU(N) Bethe equations with ideal k-strings.

Conventions and derivation
==========================
Put z_j=x+i(k+1-2j)/2, j=1,...,k. Multiplying the bare source
factors and canceling numerator/denominator factors gives

    product_j e_k(z_j) = product_{q=1}^k e_{2q-1}(x).

Multiplying scattering factors between two strings gives

    product_{j,l} e_2(z_j-w_l)
        = e_{2k}(x-y) product_{q=1}^{k-1} e_{2q}(x-y)^2,
    product_{j,l} e_{-1}(z_j-w_l)
        = product_{q=1}^k e_{2q-1}(x-y)^(-1).

These are rational identities after cancellations. Internal scattering
cancels pairwise, e_2(u)e_2(-u)=1, before taking the ideal-string limit.
Define theta_n(x)=2 atan(2x/n), theta'_n=n/(x*x+n*n/4), and

    S(x) = sum_{q=1}^k theta_{2q-1}(x),
    T(x) = 2 sum_{q=1}^{k-1} theta_{2q}(x)+theta_{2k}(x).

Since e_n(x)=-exp(-i theta_n(x)) for real x and positive n, the
continuous logarithmic equations are

    L delta_{r1} S(x_ra) - sum_{b != a} T(x_ra-x_rb)
        + sum_{s=r-1,r+1} sum_b S(x_ra-x_sb) = 2 pi I_ra,
    n_r=L(N-r)/N, I_ra=a-(n_r-1)/2, a=0,...,n_r-1.

Only levels 1,...,N-1 exist. The asymptotic counting-function endpoint
is (n_r+2k-1)/2, using n_{r-1}+n_{r+1}=2n_r at inner levels
and L+n_2=2n_1 at level 1. The above I's have the required integer
or half-integer parity. At k=1,N=2 this is L theta_1-sum theta_2=2 pi I.

For the energy, logarithmically differentiate the source identity:

    sum_j k/(z_j*z_j+k*k/4)
        = sum_{q=1}^k (2q-1)/(x*x+(2q-1)^2/4) = S'(x).
    E_ideal = L H_k - (1/2) sum_a S'(x_1a).

Thus the fused expression sums all member energies without complex
rounding. For one bare root, lambda=(k/2)cot(p/2) gives
Delta E=-(1-cos(p))/k, fixing the energy normalization.

Newton uses an analytic dense Jacobian and backtracking. Starting
points are quantiles of the thermodynamic centre density, not a fit
to energies. To obtain that density, write q=exp(-abs(omega)/2).
The Fourier kernels obey
    S'/(2 pi) -> q(1-q^(2k))/(1-q^2),
    1+T'/(2 pi) -> (1+q^2)(1-q^(2k))/(1-q^2).
Canceling the common factor in the density equations gives
    2 cosh(omega/2) rhohat_r-rhohat_{r-1}-rhohat_{r+1}=delta_{r1},
    rhohat_r=sinh((N-r)omega/2)/sinh(N omega/2).
Consequently rho_r(x)=sin(pi*r/N)/[N(cosh(2*pi*x/N)-cos(pi*r/N))].
Its integral is (N-r)/N. Inverting its cumulative distribution gives
the initial_centres function below. This does not assume the proposed
bulk energy or the proposed conformal finite-size coefficient.

Validation stops production if the ideal-string errors are already
of order 1/L in total energy. The supplied ED numbers are rounded,
so agreement can only be tested at their displayed decimal precision.
Only Python's standard library, NumPy and SciPy are required.
"""

import json
import numpy as np
from scipy.linalg import solve


# Keep strings to preserve each ED target's stated decimal precision.
TARGETS = {
    (2, 1): [(8, '2.3489065911'), (12, '3.6126090826'), (16, '4.8577036394')],
    (2, 2): [(4, '1.6492189406'), (8, '3.8393150668'), (12, '5.8948058277')],
    (2, 3): [(8, '4.9267655562'), (12, '7.5551867405')],
    (2, 4): [(8, '5.7817491256'), (10, '7.3280616238')],
    (3, 1): [(6, '0.6972243623'), (9, '1.2102821205'), (12, '1.6877365977')],
    (3, 2): [(6, '1.2020467253'), (9, '2.0739788158')],
    (3, 3): [(6, '1.5976165898')],
    (4, 1): [(8, '0.53949513'), (12, '0.9445081373')],
    (4, 2): [(8, '0.9558401418')],
}


def phase(x, terms, derivative=False):
    out = np.zeros_like(x, dtype=float)
    for n, weight in terms:
        if derivative:
            out += weight * n / (x*x + n*n/4.0)
        else:
            out += weight * 2.0 * np.arctan(2.0*x/n)
    return out


def initial_centres(N, L):
    blocks = []
    for r in range(1, N):
        n = L*(N-r)//N
        u = (np.arange(n) + 0.5)/L - (N-r)/(2.0*N)
        blocks.append(N/np.pi * np.arctanh(
            np.tan(np.pi*u) * np.tan(np.pi*r/(2*N))))
    return np.concatenate(blocks)


def check_fusion():
    """Check fusion and the energy sum against the unfused expressions."""
    worst = 0.0
    for k in range(1, 5):
        shifts = 0.5j*(k+1-2*np.arange(1, k+1))
        x, d = 0.37, 0.83

        def e(n, z):
            return (z+0.5j*n)/(z-0.5j*n)

        source = np.prod(e(k, x+shifts))
        fused_source = np.prod([e(n, x) for n in range(1, 2*k, 2)])
        pair_differences = d+shifts[:, None]-shifts[None, :]
        same = np.prod(e(2, pair_differences))
        fused_same = e(2*k, d)*np.prod([e(n, d)**2 for n in range(2, 2*k, 2)])
        adjacent = np.prod(e(-1, pair_differences))
        fused_adjacent = np.prod([e(n, d)**(-1) for n in range(1, 2*k, 2)])
        members = np.sum(k/((x+shifts)**2+k*k/4.0))
        fused_energy = sum(n/(x*x+n*n/4.0) for n in range(1, 2*k, 2))
        p = 0.71
        lam = (k/2.0)/np.tan(p/2.0)
        defect_error = -0.5*k/(lam*lam+k*k/4.0)+(1-np.cos(p))/k
        worst = max(worst, abs(source-fused_source), abs(same-fused_same),
                    abs(adjacent-fused_adjacent), abs(members-fused_energy),
                    abs(defect_error))
    if worst > 1e-12:
        raise RuntimeError(f'Fusion or energy normalization check failed: {worst}')
    return worst


def solve_centres(N, k, L, tol=2e-13):
    if N < 2 or k < 1 or L < N or L % N:
        raise ValueError('Require N>=2, k>=1, L>=N, and L divisible by N')
    counts = [L*(N-r)//N for r in range(1, N)]
    ends = np.cumsum([0] + counts)
    slices = [slice(ends[r], ends[r+1]) for r in range(N-1)]
    source = [(n, 1) for n in range(1, 2*k, 2)]
    same = [(n, 2) for n in range(2, 2*k, 2)] + [(2*k, 1)]
    quantum = np.concatenate([np.arange(n)-(n-1)/2.0 for n in counts])

    def equations(x, jacobian=True):
        f = -2*np.pi*quantum.copy()
        J = np.zeros((len(x), len(x))) if jacobian else None
        for r, sl in enumerate(slices):
            xr = x[sl]
            diagonal = np.zeros(len(xr))
            if r == 0:
                f[sl] += L*phase(xr, source)
                if jacobian:
                    diagonal += L*phase(xr, source, True)
            d = xr[:, None] - xr[None, :]
            f[sl] -= phase(d, same).sum(axis=1)
            if jacobian:
                kernel = phase(d, same, True)
                # The self term T(x_a-x_a) has derivative zero.
                np.fill_diagonal(kernel, 0.0)
                diagonal -= kernel.sum(axis=1)
                J[sl, sl] = kernel
            for s in (r-1, r+1):
                if 0 <= s < N-1:
                    ss = slices[s]
                    d = xr[:, None] - x[ss][None, :]
                    f[sl] += phase(d, source).sum(axis=1)
                    if jacobian:
                        kernel = phase(d, source, True)
                        diagonal += kernel.sum(axis=1)
                        J[sl, ss] = -kernel
            if jacobian:
                indices = np.arange(ends[r], ends[r+1])
                J[indices, indices] += diagonal
        return f, J

    x = initial_centres(N, L)
    for iteration in range(80):
        f, J = equations(x)
        residual = np.max(np.abs(f))
        if residual < tol:
            break
        step = solve(J, -f, assume_a='sym')
        alpha = 1.0
        for _ in range(40):
            trial = x + alpha*step
            ordered = all(np.all(np.diff(trial[sl]) > 0) for sl in slices)
            if ordered and np.max(np.abs(equations(trial, False)[0])) < residual:
                x = trial
                break
            alpha *= 0.5
        else:
            raise RuntimeError(f'Newton line search failed: {(N,k,L)}, {residual=}')
    else:
        raise RuntimeError(f'Newton did not converge: {(N,k,L)}')
    symmetry = max(np.max(np.abs(x[sl] + x[sl][::-1])) for sl in slices)
    if symmetry > 1e-11:
        raise RuntimeError(f'Lost reflection symmetry: {symmetry}')
    Hk = np.sum(1.0/np.arange(1, k+1))
    energy = L*Hk - 0.5*phase(x[slices[0]], source, True).sum()
    return float(energy), float(residual), float(symmetry)


def main():
    rows = []
    fusion_error = check_fusion()
    print('Ideal k-string validation: Delta = E_ideal - E_ED')
    print(' N  k   L        E_ideal             E_ED           Delta')
    for (N, k), targets in TARGETS.items():
        for L, ed_text in targets:
            energy, residual, symmetry = solve_centres(N, k, L)
            ed = float(ed_text)
            delta = energy-ed
            rounding = 0.5*10.0**(-len(ed_text.split('.')[1]))
            rows.append(dict(N=N, k=k, L=L, E_ideal=energy, E_ED=ed,
                             ED_text=ed_text, delta=delta, residual=residual,
                             reflection_error=symmetry, ED_rounding=rounding))
            print(f'{N:2d} {k:2d} {L:3d} {energy:18.12f} {ed:16.10f} {delta:+15.9e}')

    exact = [r for r in rows if r['k'] == 1]
    passed = all(abs(r['delta']) <= r['ED_rounding']+2e-12 for r in exact)
    print('\nk=1 agrees within ED rounding plus 2e-12:', passed)
    print(f'Maximum fusion / member-energy / one-defect error: {fusion_error:.3e}')
    print(f'Maximum logarithmic-equation residual: {max(r["residual"] for r in rows):.3e}')
    print(f'Maximum reflection error: {max(r["reflection_error"] for r in rows):.3e}')

    print('\nFinite-size errors for every family with two or more sizes:')
    print(' N  k   L       L*Delta        L^2*Delta    p_eff')
    print('p_eff = -ln(|Delta(L2)/Delta(L1)|)/ln(L2/L1)')
    for N, k in TARGETS:
        group = [r for r in rows if r['N'] == N and r['k'] == k]
        if k == 1 or len(group) < 2:
            continue
        previous = None
        for r in group:
            L, delta = r['L'], r['delta']
            p = ('    n/a' if previous is None else
                 f'{-np.log(abs(delta/previous["delta"]))/np.log(L/previous["L"]):7.4f}')
            print(f'{N:2d} {k:2d} {L:3d} {L*delta:14.9f} {L*L*delta:14.9f} {p}')
            previous = r

    print('\nAll listed systems have ordered real centre solutions.')
    print('For N=2,k=2 the largest supplied comparison is L=12.')
    print('The observed total-energy error is consistent with 1/L,')
    print('hence an energy-density error of order 1/L^2.')
    print('These few small sizes do not prove an asymptotic exponent.')
    print('The error contaminates the conformal 1/L term itself.')
    print('VALIDATION FAILED for finite-size ground-state extrapolation.')
    print('Production and A_infinity fits were not run.')
    with open('bethe_strings_validation.json', 'w', encoding='utf-8') as stream:
        json.dump(dict(status='validation_failed_no_production',
                       k1_agrees_at_ED_precision=passed, data=rows), stream, indent=2)
        stream.write('\n')
    print('Wrote bethe_strings_validation.json (validation energies only).')
    if not passed:
        raise RuntimeError('k=1 validation failed even at the supplied ED precision')


if __name__ == '__main__':
    main()
