# Controlled quartet amplitudes and the limits of a simple detuning

Date: 2026-11-18. Author and AI assistance: OpenAI Codex, GPT-6 Astra
(Max reasoning requested for this assignment). No independent human or
second-model review is claimed.

## Claim and status

**A quartet source is nonzero in controlled finite examples of the actual
nearest-neighbour qubit family. A special leading cancellation is lifted at
the next microscopic order. Neither calculation yet establishes a stable
two-mode thermodynamic phase or decides V4.**

The physical family is

\[
 H=\sum_j[aX_jY_{j+1}+bY_jX_{j+1}+uZ_jZ_{j+1}-hZ_j],
\]

with microscopic d=2, Hermiticity, range one and one-site translation
invariance. There are four independent real coefficients before removing
one common positive scale, and four literal Pauli monomials in this basis.
No blocking, auxiliary states or L-dependent thermodynamic coupling
schedule is proposed. The desired infrared quantity remains unitary compact
Virasoro c, allowing a sum of Gaussian contributions with distinct
velocities; no value of c or c_eff is established here.

The new results are:

1. **[derived exact finite matching; numerical evaluations]** The full
   interacting spectral sum for a resonant N to N+4 transition is evaluated
   with every intermediate N+2 state in the appropriate momentum and
   symmetry sector. Both generated attractions and the blockade are kept.
   At L=14, kappa=2, N=0 to 4, its dimensionless magnitude is
   **0.06381310620545**, giving an off-diagonal term
   **(8/9) times this number times Delta^4/V^3**. Direct diagonalization of
   the literal qubit ring approaches the predicted coefficient with relative
   errors 6.30%, 2.53%, 1.13% at Delta/V=1/10, 1/16, 1/24. These errors
   decrease as (Delta/V)^2. Nonzero background N=2 to 6 calculations are
   included. These are finite spectral amplitudes, **not a measurement of
   the continuum C4 or its dimension at V4**.

2. **[proved, exact finite perturbation algebra]** At L=12, kappa=1, the
   uniform four-particle same-sublattice state is dark to the leading pair
   source. It remains exactly dark in the truncated leading model at any
   source strength. In the **full microscopic Hamiltonian**, after matching
   the crossing field consistently, its first nonzero off-diagonal coupling
   is

   \[
   \boxed{\frac{6\sqrt{30}}5\frac{\Delta^6}{V^5}
     +O(\Delta^8/V^7).}
   \]

   Exact rational matrix elimination gives the coefficient and the field
   counterterm; independent floating-point qubit spectra approach it.
   This is a finite dark-state control, not an extended liquid or a useful
   symmetry protection for c=2. See [AMPLITUDE.md](AMPLITUDE.md).

3. **[derived restriction; tested leading kinematics]** Changing only
   a,b,u,h preserves the symmetry that fixes the balanced valley sum to pi.
   It cannot detune the quartet momentum. Its leading separation-three
   source cannot be canceled by changing the imaginary hop or the cusp
   chemical potential. Adding just **v XX** breaks that symmetry and
   detunes the bare quartet momentum. A rational five-term example has
   two bare occupied wells and removes the particular momentum-pi dimer's
   energetic advantage over two lower-valley atoms. It also generates a
   separation-two pair source and worsens the literal coefficient spread
   to **4064**. A matched Gaussian matrix, sector stability and the
   finite-density locking thresholds are still missing. This example is
   **not a stabilized c=2 candidate**. See [DETUNING.md](DETUNING.md).

4. **[implemented, calibrated controls]** The connected momentum-pi pair
   diagnostic subtracts the complex disconnected expectation and retains
   phases. On exact two-leg XX correlations it returns pair exponents
   **1.00119, 1.00062, 1.00038** in three stated windows, versus the exact
   value 1. The stability routine exhausts the **entire integer lattice**
   using a proven eigenvalue bound; it catches both a (1,1) and a (3,1)
   mixed instability despite irrelevant pure vertices. No existing V4
   correlation data were loaded or fitted. See [DIAGNOSTICS.md](DIAGNOSTICS.md).

## What the amplitude does and does not fix

The spectral convention uses eigenstates of the *same full leading
interacting Hamiltonian* as the external dressed fields. There is no free
denominator substitution and no bare-band external state. The phase of a
matrix element depends on the phases of those eigenstates; its magnitude
and the avoided-crossing gap do not. An on-shell connected-source identity
and an exactly canceling quadratic-fermion control distinguish this
calculation from disconnected virtual-pair dressing.

Finite matching requires tuning a crossing field separately for each
finite calculation. Those fields are **measurement devices**, not a
proposal to change a microscopic Hamiltonian with chain length. The table
does not extrapolate at fixed density or extract a cosine normalization.
In particular the quartet's connected pi-pair exponent at V4 remains
unestablished. The complete symmetry-allowed lattice from the
[locking audit](../v4_locking_operators_20261118/README.md) remains necessary.

## Reproducing and checking

Dependencies: numpy, scipy, threadpoolctl; sympy for the exact certificate.
From this directory, with one worker:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B quartet.py --output quartet.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B microscopic.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B exact_dark.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B detuning.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B diagnostics.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B checks.py
```

Each script writes the correspondingly named JSON. The full bounded suite
takes a few seconds on the working machine; per-script runtimes are saved.
The largest microscopic matrix is 330 by 330 (a symmetry sector of a
14-site ring); the largest matching endpoint sector is 561 by 561.
The two-body check uses a 512-site *relative-coordinate* tridiagonal-plus-
second-diagonal matrix. No DMRG, MPS checkpoint, background process,
production job or large parameter scan was used. Battery-saving limits
remain in force.

The derivations and independently coded checks were performed by the same
author, so they are not labeled independently [checked]. Source files in
the earlier contributions were read and left unchanged. The coordinator
owns shared status pages; [SHARED_SUMMARY.md](SHARED_SUMMARY.md) supplies
proposed text without editing them.

## Most decisive next calculation

At one **fixed** physical Hamiltonian and a converged finite-density bulk
state, match the connected momentum-pi pair exponent, the other entries of
Q, and the quartet matrix element in the same dressed convention. If the
pair exponent is below one and the quartet survives that matching, it is
relevant. If it is above one, test the full lattice and the relative mode
before claiming stability. For the five-term detuning, measure dressed
momentum/density mismatches and susceptibilities together; bare well
occupancy and a nonzero mismatch alone cannot answer the question.

The immediate deliverable is an explicit amplitude calculation and an
exact next-order cancellation test, not a claimed c=2 phase or a d=2
upper bound.
