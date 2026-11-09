"""Literal 4x4 nearest-neighbour bonds for the facilitated cusp search.

These are candidate test points, not measured or certified counterexamples.
The q values balance dressed valley energies only through order r^2.
"""
import json
from pathlib import Path
import numpy as np


def candidate(V, r, mu_eff, q, Delta=1.):
    identity = np.eye(2)
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.])
    n = (identity-z)/2
    create = (x-1j*y)/2
    destroy = create.conj().T
    kron = np.kron
    lam = r*Delta
    scale = Delta**2/(2*V)
    w = -lam**2/V+scale*(q*r*r+1j)
    mu = (3*Delta**2-2*lam**2)/V+scale*mu_eff
    jx, jy = (Delta+w.real)/2, (w.real-Delta)/2
    dz, dy, jz = w.imag/2, lam/2, V/4
    hz = (mu-V)/2
    bond = (jx*kron(x, x)+jy*kron(y, y)+jz*kron(z, z)
            + dz*(kron(x, y)-kron(y, x))
            + dy*(kron(x, z)-kron(z, x))
            + hz/2*(kron(z, identity)+kron(identity, z)))
    physical = (Delta*(kron(create, create)+kron(destroy, destroy))
                + w*kron(destroy, create)+w.conjugate()*kron(create, destroy)
                + V*kron(n, n)-mu/2*(kron(n, identity)+kron(identity, n))
                + lam*(kron(n, x)-kron(x, n)))
    # Scalar and one-site divergence disappear from periodic-ring dynamics.
    scalar = V/4-mu/2
    divergence = lam/2*(kron(identity, x)-kron(x, identity))
    error = np.max(np.abs(physical-bond-scalar*np.eye(4)-divergence))
    herm = np.max(np.abs(bond-bond.conj().T))
    assert error < 1e-12 and herm < 1e-12
    return {"V": V, "Delta": Delta, "r": r, "mu_eff": mu_eff, "q": q,
            "effective_energy_unit": scale,
            "couplings": {"Jx": jx, "Jy": jy, "Jz": jz,
                          "Dz_XY_minus_YX": dz, "Dy_XZ_minus_ZX": dy,
                          "hz_total_field": hz},
            "w_bare": [w.real, w.imag], "mu_bare": mu,
            "bond_basis": ["00", "01", "10", "11"],
            "bond_real": bond.real.tolist(), "bond_imag": bond.imag.tolist(),
            "hermiticity_error": float(herm),
            "spin_occupation_identity_error": float(error)}


def main():
    points = [candidate(V, r, mu, q)
              for V in [40., 80.]
              for mu, q in [(-2.30, -0.364132190384574),
                            (-2.27, -0.4746147911696791)]
              for r in [0., .05, .10, .15]]
    result = {"status": "candidate test family, no central-charge claim",
              "definition": "H=sum_j h_(j,j+1), periodic spin chain; Pauli convention",
              "central_charge": None, "points": points}
    path = Path(__file__).with_name('microscopic_candidates.json')
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"Verified {len(points)} Hermitian 4x4 bonds and occupation-basis identities.")


if __name__ == '__main__':
    main()
