"""Bounded complete-sector check of a number-changing escape candidate.

No central charge, thermodynamic stiffness, or compact CFT is inferred.
"""
from __future__ import annotations
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import json
import time
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy import linalg
from threadpoolctl import threadpool_limits
from normal_map import expected, A


def inspect(L, a, zeta, mu):
    states, h = expected(L, a=a, zeta=zeta, mu=mu)
    energy, vectors = linalg.eigh(h.toarray())
    multiplicity = int(np.sum(energy-energy[0] < 1e-9))
    grounds = vectors[:, :multiplicity]
    probs = np.mean(abs(grounds)**2, axis=1)
    number = np.array([int(s).bit_count() for s in states])
    nmean = float(probs@number)
    nvar = float(probs@(number**2)-nmean**2)
    polarized = np.array([len({j % 3 for j in range(L) if int(s) >> j & 1}) <= 1 for s in states])
    multi_probability = float(probs@(number>=2))
    conditional_single = (float(probs@(polarized & (number>=2)))/multi_probability
                          if multi_probability > 1e-12 else None)
    density_mode = np.array([sum(np.exp(2j*np.pi*j/3) for j in range(L) if int(s) >> j & 1) for s in states])
    elements = vectors[:, multiplicity:].conj().T@(number[:,None]*grounds)
    gap = energy[multiplicity:]-energy[0]
    susceptibility = float(2*np.sum(abs(elements)**2/gap[:,None])/L/multiplicity)
    ground_number = grounds.conj().T@(number[:,None]*grounds)
    commutator = h@sparse_number(number)-sparse_number(number)@h
    return dict(L=L, dimension=len(states), a=float(a), zeta=float(zeta), mu=float(mu),
                energy=float(energy[0]), first_gap=float(energy[1]-energy[0]),
                ground_multiplicity=multiplicity,
                first_distinct_gap=float(min(gap)),
                density=nmean/L, number_variance=nvar,
                chemical_response_per_site=susceptibility,
                ground_number_nonscalar_error=float(np.max(abs(ground_number-nmean*np.eye(multiplicity)))),
                single_residue_support_probability=float(probs@polarized),
                at_least_two_particles_probability=multi_probability,
                single_residue_probability_given_at_least_two=conditional_single,
                period_three_density_structure_per_site=float(probs@abs(density_mode)**2/L),
                number_commutator_max=float(max(abs(commutator.data), default=0)))


def sparse_number(number):
    from scipy import sparse
    return sparse.diags(number,dtype=float)


def physical_point(epsilon):
    """Literal rational six-Pauli-term representative, including odd jets."""
    e=Fraction(epsilon)
    u=2*e**2-34*e**4-256*e**5-Fraction(1012,3)*e**6
    t=-e**2+8*e**4+32*e**5-Fraction(5968,75)*e**6
    p=3*e**2+16*e**3
    coefficients=dict(ZZ=Fraction(1,4), Z=-(1+u)/2, XZ=e/2, ZX=-e/2,
                      XX=(t+p)/2, YY=(t-p)/2)
    spread=max(map(abs,coefficients.values()))/min(map(abs,coefficients.values()))
    identity=np.eye(2); x=np.array([[0,1],[1,0]])
    y=np.array([[0,-1j],[1j,0]]); z=np.diag([1,-1])
    paulis={'X':x,'Y':y,'Z':z}
    canonical=sum(float(v)*np.kron(paulis[k[0]],paulis[k[1]])
                  for k,v in coefficients.items() if len(k)==2)
    canonical+=float(coefficients['Z'])/2*(np.kron(z,identity)+np.kron(identity,z))
    scalar=float(Fraction(1,4)+u/2)*np.eye(4)
    telescope=float(e)/2*(np.kron(identity,x)-np.kron(x,identity))
    literal=A.bond(v=1,u=float(u),t=float(t),p=float(p),a=float(e),b=-float(e))
    error=float(np.max(abs(literal-canonical-scalar-telescope)))
    assert error < 1e-14
    return dict(epsilon=str(e),u=str(u),t=str(t),p=str(p),
                pauli_coefficients={k:str(v) for k,v in coefficients.items()},
                coefficient_spread_exact=str(spread),coefficient_spread=float(spread),
                pauli_identity_error=error,emergent_sixth_scale=float(Fraction(128,3)*e**6))


if __name__ == '__main__':
    start=time.monotonic()
    with threadpool_limits(1):
        # Two declared chemical tests at three bounded sizes; no scan.
        result = dict(candidate=dict(a=2, zeta="1/100", mu="1/2"),
                      low_density_control=[inspect(L, 2, 1/100, 3/100) for L in (12, 15, 18)],
                      checks=[inspect(L, 2, 1/100, 1/2) for L in (12, 15, 18)],
                      physical_points=[physical_point(e) for e in ('1/20','1/5')])
    result['runtime_seconds']=time.monotonic()-start
    out=Path(__file__).with_name('escape_results.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
