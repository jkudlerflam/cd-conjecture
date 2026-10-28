# Choosing the core accuracy and the physical coupling in the right order

## Status

[proved] The estimates below give a consistent hierarchy of *finite-order*
errors for the sparse first-order construction. [conditional] Concluding
an infinite-volume Gaussian phase also uses the local low-energy and RG
inputs listed in [IR_AUDIT.md](IR_AUDIT.md). An extensive finite-volume
perturbation estimate is not promoted to a bound on conformal level gaps.

Derived by Codex, GPT-6 family, 2026-10-28, during the research continuation.
No human review is claimed.

The primary construction now uses [exact algebraic preparation](../exact_route/EXACT_AMPLITUDES.md),
so its circuit approximation error is zero. The first section records
an optional approximate implementation; the coupling hierarchy applies
to the exact construction without that additional error.

## 1. Approximation does not have to resolve the clock gap

Fix a prime p, its incidence graph, and the corresponding finite alphabet.
Write N for its number of edges and u for its number of upper vertices.
The desired two-port isometry W has hopping visibility

    eta = (1-1/(2N))/(4uN).

Let W' be an exactly isometric circuit approximation, with
||W'-W||<=delta on the intended input space. This can be realized by a
finite word of exact, charge-preserving elementary unitaries. In the
history parent, these actual unitaries still give an exact zero-energy
code: approximating the readout is not the same operation as perturbing
the parent Hamiltonian or splitting its code.

For every norm-one observable A,

    ||W'^dagger A W' - W^dagger A W|| <= 2 delta.

On two independently encoded cells the corresponding error for any
norm-one inter-port bond is at most 4 delta. This is an operator-norm
estimate on the full bond, so it does not introduce an unnecessary
factor N from summing its matrix elements.

Suppose the history readouts have positive, input-independent weights
r_L,r_R. A link between cells has the common factor r_L r_R on both
the target and its approximation error. The ratio is consequently

    ||bond error|| / target hopping scale <= 4 delta / eta.

The same cancellation holds for an exactly applied common wrapper
visibility factor. Thus a sufficient accuracy requirement is

    delta << eta * epsilon,

where epsilon is the dimensionless interaction budget for the target.
For the projective-plane family, eta is of order p^(-5), and the
payload mismatch is of order p^(-3). Taking delta=p^(-K_accuracy)
with a sufficiently large fixed exponent makes the algebraic error
smaller than both. Its precision-bit count is only O(log p).

This observation avoids a circular resource requirement. The core gap
can be exponentially or even more severely small in its length, while
the preparation circuit still needs only polynomial precision in p.
The physical interaction scale, chosen later, absorbs the small gap.
During repeated compute/readout/uncompute cycles, the inverse of the
*actual* gate word cancels it exactly. The readout error therefore does
not accumulate over the exponentially many clock sweeps.

## 2. Normalize the intended effective bond before choosing its strength

After fixing the core and wrapper, all endpoint probabilities, helix
compensator norms, marker-hopping overlaps and gaps are finite numbers.
They may be extremely unfavorable. Absorb their inverses into a fixed
microscopic perturbation V_p whose first-order projection is the desired
payload and length-fluid Hamiltonian in units of one common scale.
Use the physical family

    H(s) = H_base + s V_p,       s>0.

Let B_p bound the local norm of V_p and Delta_p be the unperturbed gap
above the complete code manifold. B_p is not assumed to be order one.
At second order, connected corrections have local scale at most

    C_p s^2 B_p^2 / Delta_p.

Relative to the retained scale s this is
C_p s B_p^2/Delta_p. For any fixed positive tolerance epsilon_p one can
choose

    0 < s << epsilon_p Delta_p / (C_p B_p^2).

At any fixed higher perturbative order, decreasing s further makes
that order's relative correction smaller. The constants depend on the
fixed cell, not on the total physical ring length. This does not prove
convergence of the full expansion or stability of its eliminated sector.

This order of choices is essential:

1. Fix p and hence N, q and the required Gaussian/resonance tolerance.
2. Fix a finite circuit accuracy satisfying the first section.
3. Build the exact history core and the complete mobile wrapper.
4. Calculate or bound their finite gaps and microscopic visibility norms.
5. Choose one nonzero s satisfying the local hierarchy.
6. Only then take the physical ring length to infinity.

One must not send s to zero with the total ring length, nor infer a
thermodynamic CFT from agreement of a few eigenvalues at each finite L.

## 3. What the graph symmetry buys at nonzero coupling

Every microscopic term commutes with the species torus and marker number.
Effective terms obtained by eliminating gapped degrees of freedom inherit
those symmetries. Hence no arbitrarily small symmetry-breaking mass needs
to be tolerated: it is absent exactly. The graph forbids nonzero neutral
flavour transfers below six fermion fields. Its remaining forward
scattering changes the Gaussian matrices continuously.

At fixed p the number of potentially relevant or marginal oscillatory
channels is finite, their mismatches have a positive minimum, and the
free Gaussian form is strictly positive. If a local weak-coupling
effective action exists and depends continuously on s, these strict
inequalities survive for sufficiently small s. The resulting prediction
is N+1 compact-boson modes. The audit identifies precisely which parts
of that last sentence are standard infrared input rather than theorems
proved by this construction.
