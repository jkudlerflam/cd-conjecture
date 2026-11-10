"""C190 model definitions (Claude Opus 5.5 agent, 2026-11-09).

Qubit chain H = sum_j [J_x X X + J_y Y Y + D (X Y - Y X) + J_z Z Z]_{j,j+1} + h_z sum_j Z_j on an open chain,
with the field on every site (the collaborator's convention).  A constant e_shift per site is added so that
the energy of the empty state is close to zero (better conditioned local solves); E_phys = E - L e_shift.

Fermion form (C173 PROOF eq. 0.1): w = J_x + J_y + 2 i D, Delta = J_x - J_y, V = 4 J_z, mu = 4 J_z + 2 h_z.
Second order (C173 PROOF 3.1): |t_2| = Delta^2 / (2V), mu_eff = (mu - 3 Delta^2 / V) / |t_2|.
"""
import numpy as np
import z2dmrg as Z


def bond_family(name):
    """Return dict(Jx, Jy, Jz, D, hz0) for the named bond.  hz0 is a reference field (tuned later)."""
    if name == 'p153':            # collaborator's phi = 1.53 parity-family bond (unrounded)
        w = np.exp(1.53j) / 80
        return dict(Jx=(1 + w.real) / 2, Jy=(-1 + w.real) / 2, Jz=10.0, D=w.imag / 2, hz0=-19.97589125)
    if name == 'rat':             # collaborator's rational simplification
        return dict(Jx=2001 / 4000, Jy=-1999 / 4000, Jz=10.0, D=1 / 160, hz0=-199759 / 10000)
    if name == 'sym40':           # symmetric line phi = pi/2, |w| = |t2|, V/Delta = 40 (C173 candidate bond)
        return dict(Jx=0.5, Jy=-0.5, Jz=10.0, D=0.00625, hz0=-19.97625)
    if name == 'sym80':           # symmetric line, V/Delta = 80: |t2| = 1/160, D = |t2|/2
        return dict(Jx=0.5, Jy=-0.5, Jz=20.0, D=1 / 320, hz0=-39.988125)
    raise ValueError(name)


def eff_params(b, hz):
    Delta = b['Jx'] - b['Jy']
    V = 4 * b['Jz']
    w = b['Jx'] + b['Jy'] + 2j * b['D']
    t2 = Delta ** 2 / (2 * V)
    mu = 4 * b['Jz'] + 2 * hz
    return dict(Delta=Delta, V=V, t2=t2, w_over_t2=w / t2, absw=abs(w) / t2, phi=float(np.angle(w)),
                mu=mu, mu_eff=(mu - 3 * Delta ** 2 / V) / t2)


def hz_for_mu_eff(b, mu_eff):
    Delta = b['Jx'] - b['Jy']
    V = 4 * b['Jz']
    t2 = Delta ** 2 / (2 * V)
    mu = mu_eff * t2 + 3 * Delta ** 2 / V
    return (mu - 4 * b['Jz']) / 2


def candidate(name, L, hz, gauge=None):
    """MPO and site phases for a candidate bond.  gauge='stag' (only for J_x + J_y = 0) makes H real."""
    b = bond_family(name)
    h2 = Z.pauli_bond(Jx=b['Jx'], Jy=b['Jy'], Jz=b['Jz'], D=b['D'])
    gph = Z.gauge_phases(L, gauge)
    if gauge:
        h2s = [Z.rotate_bond(h2, gph[j], gph[j + 1]) for j in range(L - 1)]
        for h in h2s:
            assert np.isrealobj(h), 'gauge did not make the bond real'
    else:
        h2s = [h2] * (L - 1)
    shift = -(b['Jz'] + hz)
    h1 = hz * Z.ZZ_ + shift * Z.I2
    mpo = Z.MPO(L, h2s, [h1] * L)
    return mpo, gph, shift, dict(bond=b, h2=h2, hz=hz, h1=h1)


def xx_control(L, mu, legs=1):
    """legs=1: H = -sum (Sp Sm + Sm Sp)_{j,j+1} - mu sum n_j  (free fermions, c = 1).
    legs=2: two decoupled XX legs: -sum (Sp_j Sm_{j+2} + h.c.) - mu sum n_j (c = 2), no nearest-neighbour term."""
    n_op = np.diag([0.0, 1.0])
    h1 = -mu * n_op
    if legs == 1:
        hop = -(np.kron(Z.SP, Z.SM) + np.kron(Z.SM, Z.SP))
        return Z.MPO(L, [hop] * (L - 1), [h1] * L), np.ones(L, dtype=complex)
    zero = np.zeros((4, 4))
    return Z.MPO(L, [zero] * (L - 1), [h1] * L, nnn=[(-1.0, -1.0)] * (L - 2)), np.ones(L, dtype=complex)
