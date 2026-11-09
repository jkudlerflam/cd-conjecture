"""Exploratory small-u scaling of the facilitated binding derivative.

This computes derivatives at alpha=0, not finite-alpha unbinding.  The
largest gap box has 998,001 states, represented sparsely.  No raw vectors
are retained.  Requires the same NumPy/SciPy setup as the other scripts.
"""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded
from scipy.sparse import eye
from scipy.sparse.linalg import LinearOperator, cg

from facilitated_pair_shift import build, pair_source


def atom_derivative(u, mu, k, cutoff=8192):
    energy = -2 - u * u / 4
    a, b = 2 * u * np.sin(k / 2), -2 * np.cos(k)
    bands = np.zeros((5, cutoff))
    bands[2, :] = -energy - mu
    bands[2, 0] -= 2
    bands[2, 1] -= 2 / 3
    bands[1, 1:] = a
    bands[3, :-1] = a
    bands[0, 2:] = b
    bands[4, :-2] = b
    source = np.zeros(cutoff)
    source[0] = 4 * np.sin(k)
    solution = solve_banded((2, 2), bands, source)
    return {
        "k": float(k),
        "energy_derivative": float(-source @ solution),
        "admixture_derivative": float(solution @ solution),
        "relative_cutoff": cutoff,
    }


def run():
    output = {
        "scope": "derivatives at alpha=0; no uniformly controlled unbinding result",
        "chemical_potential": "epsilon_min-0.1*u^4, before alpha-dependent retuning",
        "runs": [],
    }
    for u, cutoffs in ((1.0, (240,)), (0.5, (160, 240, 360, 600)), (0.25, (240, 400, 600, 1000))):
        mu = -2 - u * u / 4 - 0.1 * u**4
        ka, kb = -np.arcsin(u / 4), np.pi + np.arcsin(u / 4)
        atoms = [atom_derivative(u, mu, k) for k in (ka, kb)]
        atom_sum = sum(v["energy_derivative"] for v in atoms)
        result = {"u": u, "mu": mu, "atoms": atoms, "pair_boxes": []}
        for cutoff in cutoffs:
            info, source, direct = pair_source(cutoff, u)
            h3 = build(3, np.pi, 1j * u, cutoff)
            denominator = h3 + (-mu - info["energy"]) * eye(h3.shape[0])
            inverse_diag = 1 / denominator.diagonal()
            preconditioner = LinearOperator(
                denominator.shape, matvec=lambda x: inverse_diag * x, dtype=complex
            )
            solution, flag = cg(
                denominator, source, M=preconditioner, rtol=1e-10, atol=0, maxiter=3000
            )
            if flag:
                raise RuntimeError(f"CG did not converge: u={u}, R={cutoff}, flag={flag}")
            pair_sigma = float(direct - np.vdot(source, solution).real)
            binding_sigma = atom_sum - pair_sigma
            formal_alpha = float(np.sqrt(-info["binding"] / binding_sigma))
            row = {
                "cutoff_each_gap": cutoff,
                "dimension": denominator.shape[0],
                "pair_derivative": pair_sigma,
                "pair_admixture_derivative": float(np.vdot(solution, solution).real),
                "binding_derivative": binding_sigma,
                "binding_derivative_over_u_squared": binding_sigma / u**2,
                "original_binding": info["binding"],
                "formal_quadratic_zero_alpha": formal_alpha,
                "formal_alpha_over_u": formal_alpha / u,
                "linear_solve_residual": float(np.linalg.norm(denominator @ solution - source)),
            }
            result["pair_boxes"].append(row)
            print(f"u={u}, R={cutoff}: binding derivative={binding_sigma:.12f}, formal alpha/u={formal_alpha/u:.8f}", flush=True)
        output["runs"].append(result)
    return output


if __name__ == "__main__":
    target = Path(__file__).with_name("small_u_scaling.json")
    target.write_text(json.dumps(run(), indent=2) + "\n")
    print(target)
