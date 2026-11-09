"""Atom self-energy for the complete facilitated extension, at order r^2.

Original C173 gauge: nearest hop i, next-nearest hop -1, chemical
potential mu. G=2 sum_j (n_(j+2)-n_(j-2)) X_j is the coefficient of r.
G maps a Bloch atom of momentum K to a separation-two pair with amplitude
4i sin K. This script independently checks that source on finite rings,
then evaluates the infinite-line limit using convergent banded resolvents.
No many-body phase or finite-r binding threshold is certified.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eig_banded, solveh_banded


def source_check(L=13):
    pairs = [(a, b) for a in range(L) for b in range(a+1, L)
             if (b-a) not in (1, L-1)]
    index = {pair: i for i, pair in enumerate(pairs)}
    worst = 0.
    for m in range(L):
        k = 2*np.pi*m/L
        source = np.zeros(len(pairs), dtype=complex)
        predicted = np.zeros_like(source)
        for x in range(L):
            atom = np.exp(1j*k*x)/np.sqrt(L)
            source[index[tuple(sorted((x, (x-2) % L)))]] += 2*atom
            source[index[tuple(sorted((x, (x+2) % L)))]] -= 2*atom
            pair = tuple(sorted((x, (x+2) % L)))
            predicted[index[pair]] += 4j*np.sin(k)*np.exp(1j*k*(x+1))/np.sqrt(L)
        worst = max(worst, float(np.max(np.abs(source-predicted))))
    assert worst < 1e-12
    return {"L": L, "momenta": L, "max_source_error": worst}


def atom_shift(k, mu, cutoff):
    eps = 2*np.sin(k)-2*np.cos(2*k)
    bands = np.zeros((3, cutoff))
    bands[0, :] = -eps-mu
    bands[0, :2] += [-2, -2/3]
    bands[1, :-1] = 2*np.sin(k/2)
    bands[2, :-2] = -2*np.cos(k)
    src = np.zeros(cutoff)
    src[0] = 1
    solution = solveh_banded(bands, src, lower=True)
    gap = eig_banded(bands, lower=True, eigvals_only=True,
                     select='i', select_range=(0, 0))[0]
    assert gap > 0
    coefficient = -16*np.sin(k)**2*solution[0]
    return {"K": float(k), "mu": mu, "cutoff": cutoff,
            "bare_atom_energy": float(eps),
            "cutoff_resolvent_gap": float(gap),
            "boundary_resolvent": float(solution[0]),
            "two_particle_weight_over_r_squared": float(16*np.sin(k)**2*np.dot(solution, solution)),
            "energy_shift_over_r_squared": float(coefficient)}


def main():
    k_a = -np.arcsin(.25)
    k_b = np.pi+np.arcsin(.25)
    rows = []
    tunings = []
    for mu in [-2.4, -2.3, -2.27, -2.26]:
        shifts = []
        for k in (k_a, k_b):
            convergence = [atom_shift(k, mu, cutoff)
                           for cutoff in [64, 128, 256, 512, 1024]]
            assert abs(convergence[-1]['energy_shift_over_r_squared']
                       - convergence[-2]['energy_shift_over_r_squared']) < 1e-9
            rows.extend(convergence)
            shifts.append(convergence[-1]['energy_shift_over_r_squared'])
        s_a, s_b = shifts
        # Add an effective real hop q*r^2 after canceling the direct F^2 hop.
        q = -(s_a-s_b)/(4*np.cos(k_a))
        tuned_a = s_a+2*q*np.cos(k_a)
        tuned_b = s_b+2*q*np.cos(k_b)
        assert abs(tuned_a-tuned_b) < 1e-12
        tunings.append({"mu": mu, "atom_a_shift": s_a, "atom_b_shift": s_b,
                        "extra_real_hop_over_r_squared": float(q),
                        "common_tuned_atom_shift": float(tuned_a),
                        "two_atom_threshold_shift": float(s_a+s_b)})
    result = {"scope": "order-r^2 isolated-branch perturbation, numerical resolvent limits",
              "source_checks": [source_check(L) for L in [9, 13, 20]],
              "resolvent_rows": rows, "valley_retuning": tunings}
    path = Path(__file__).with_name('facilitated_atom_results.json')
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({"source_checks": result['source_checks'],
                      "valley_retuning": tunings}, indent=2))


if __name__ == '__main__':
    main()
