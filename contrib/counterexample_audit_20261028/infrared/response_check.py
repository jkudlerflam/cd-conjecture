"""Independent phase-model response and slow/fast Gaussian checks.

This uses exact diagonalization of the gap fluid, not the d=24 model.
No campaign code is imported.  Run with numpy and scipy available.
"""
import argparse
import itertools
import json
import math
import os
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import scipy.integrate as integrate
import scipy.sparse as sparse
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent


def phase_model(n):
    states = []
    for cuts in itertools.combinations(range(2*n-1), n-1):
        separators = (-1,) + cuts + (2*n-1,)
        states.append(tuple(separators[j+1]-separators[j]-1 for j in range(n)))
    index = {s: j for j, s in enumerate(states)}
    row, col = [], []
    for j, state in enumerate(states):
        for donor in range(n):
            if state[donor] == 0:
                continue
            for receiver in ((donor+1) % n, (donor-1) % n):
                target = list(state)
                target[donor] -= 1
                target[receiver] += 1
                row.append(index[tuple(target)])
                col.append(j)
    h = sparse.coo_matrix((-np.ones(len(row)), (row, col)),
                          shape=(len(states), len(states))).tocsr()
    eig, vectors = sla.eigsh(h, k=1, which="SA", tol=1e-12,
                             v0=np.ones(len(states)))
    e0, psi = float(eig[0]), vectors[:, 0]
    if psi.sum() < 0:
        psi = -psi
    error = float(np.linalg.norm(h @ psi-e0*psi))
    exact_e0 = -2 / math.sin(math.pi / (2*n))
    assert abs(e0-exact_e0) < 1e-10
    assert error < 1e-9
    return h, np.array(states), e0, psi, error


def response(h, diagonal, e0, psi, momentum, omega):
    n = diagonal.shape[1]
    mode = np.exp(-2j*math.pi*momentum*np.arange(n)/n)/math.sqrt(n)
    source = (diagonal @ mode)*psi
    source -= psi * np.vdot(psi, source)
    # Resolvent at imaginary frequency.  Twice its real part is the
    # Euclidean two-sided spectral response.  At omega=0 use a ground
    # projector to remove the exact zero eigenvalue.
    def action(vector):
        return h @ vector-e0*vector + psi*np.vdot(psi, vector)
    if omega == 0:
        operator = sla.LinearOperator(h.shape, matvec=action, dtype=complex)
        result, info = sla.cg(operator, source, rtol=1e-10, atol=1e-12)
        residual = np.linalg.norm(action(result)-source)
    else:
        operator = sla.LinearOperator(h.shape,
            matvec=lambda x: action(action(x))+omega*omega*x,
            dtype=complex)
        rhs = action(source)
        result, info = sla.cg(operator, rhs, rtol=1e-10, atol=1e-12)
        residual = np.linalg.norm(operator @ result-rhs)
    assert info == 0, info
    assert residual < 1e-8, residual
    return float(2*np.vdot(source, result).real), float(residual)


def run(n):
    started = time.monotonic()
    h, occupation, e0, psi, residual = phase_model(n)
    contact = (occupation == 0).astype(float)
    rows = []
    for mode in sorted(set((0, 1, n//2))):
        q = 2*math.pi*mode/n
        for omega in (0., q) if mode else (0.,):
            grr, err_r = response(h, occupation, e0, psi, mode, omega)
            gaa, err_a = response(h, contact, e0, psi, mode, omega)
            ll = 4/math.pi*q*q/(omega*omega+q*q) if q else 0.
            rows.append(dict(mode=mode, q=q, omega=omega,
                             G_density=grr, G_contact=gaa,
                             smooth_density_prediction=ll,
                             residual=max(err_r, err_a)))
    result = dict(material_sites=n, holes=n, physical_sites=2*n,
                  dimension=len(psi), ground_energy=e0,
                  ground_residual=residual,
                  contact_mean=float((psi*psi) @ contact.mean(axis=1)),
                  seconds=time.monotonic()-started, responses=rows)
    print(json.dumps(result), flush=True)
    return result


def analytic_checks():
    a = -(0.25+1/math.pi**2)
    coeff = a*a/(math.pi/2)
    integrals = []
    # Quantum-current contraction: explicit integral checks its finite
    # spatial first moment and its absence of an inverse slow velocity.
    for speed in (1., .1, .01, .001):
        f = lambda tau: ((1-tau*tau)*(1-speed*speed*tau*tau)
            / ((1+tau*tau)**2*(1+speed*speed*tau*tau)**2))
        value = 2*integrate.quad(f, 0, np.inf, epsabs=1e-12)[0]
        exact = 2*math.pi*speed/(1+speed)**3
        assert abs(value-exact) < 1e-10
        integrals.append(dict(payload_velocity=speed, integral=value,
                              analytic=exact))
    modes = []
    v_internal = 1e-8
    for J in (1e-12, 1e-11, 1e-10, 1e-9, 1e-8):
        u, g, vc, kc = J*v_internal, J, 1., 4.
        # B density-stiffness; A dual-current stiffness.
        # Compute the small eigenvalue using determinant/lambda_max to
        # avoid catastrophic cancellation at the deliberately huge ratio.
        trace = vc*vc+u*u
        determinant = vc*vc*u*u-vc*kc*u*g*g
        vmax2 = (trace+math.sqrt(trace*trace-4*determinant))/2
        vmin2 = determinant/vmax2
        modes.append(dict(J=J, u=u, forward_coefficient=g,
                          relative_stiffness_squared=kc*g*g/(vc*u),
                          slow_speed_squared_ratio=vmin2/(u*u),
                          stable=bool(vmin2>0)))
    cell_a = lambda length: 2*(1-math.cos(math.pi/(length+1)))
    cell_delta = lambda length: ((length+1)*cell_a(length)
                                -length*cell_a(length+1))
    selector = (cell_delta(24)+cell_delta(23))/2
    cell_density = (cell_a(24)-selector)/24
    long_cell_cutoff = math.ceil(2*selector/abs(cell_density))
    return dict(alpha=0.5-2/math.pi**2, alpha_prime=a,
                e_second=math.pi/4, K=4, velocity=1,
                density_response_bound=4/math.pi,
                static_energy_kernel_coefficient=coeff,
                current_contraction=integrals, two_speed_examples=modes,
                long_cell=dict(delta_selector=selector,
                    selected_energy_density=cell_density,
                    cutoff=long_cell_cutoff,
                    margin_per_site=abs(cell_density)/2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", type=int, nargs="+", default=[4, 6, 8, 10])
    args = parser.parse_args()
    output = dict(scope="Exact phase fluid, not the full d=24 model",
                  analytic=analytic_checks(), finite=[])
    for n in args.sizes:
        output["finite"].append(run(n))
        (HERE/"response_results.json").write_text(json.dumps(output, indent=2)+"\n")
