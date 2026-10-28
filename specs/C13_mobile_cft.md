# C13: Is the mobile-cell construction a CFT in the strict sense?

You are a long-running research agent. Your working directory is your own folder. If `PROGRESS.md` exists there, read it first and continue from where it stops. Update `PROGRESS.md` after every substantive step, and write `RESULT.md` at the end.

## Background

The repository is two levels above your folder (`../../`). Read `../../STATUS.md` first: it states the c–d conjecture and lays out a candidate construction (sites form long cells; each cell holds a degenerate logical code; neighbouring cells couple through their ports; the logical chain is tuned to an n-flavour fermion target with incommensurate Fermi momenta, c = n).

The problem this run addresses. In the version with frozen cells, the cell markers do not move, so translation is spontaneously broken: on rings with b | L the ground state is b-fold degenerate, and on other ring lengths there are defects (at least L ground states on large prime rings). That is not "the low-energy spectrum is a unitary CFT" in the strict sense the conjecture intends, which includes a unique ground state. The version meant to fix this has mobile cells: the cell crystal melts into a c = 1 "length fluid" and the payload rides on it, as in spin–charge separation, giving c = n + 1 on every ring with a unique ground state.

Material already in the repository (read these; they are Astra drafts, unchecked unless stated):
- `../work_C8/RESULT.md` §12 to §17 and `../work_C8/{ring_obstruction,mobile_selector,equal_hopping_selector,length_fluid,padding_transport,selector_wrapper}.md`: frozen-marker degeneracy, the mobile two-length selector, the exact finite-hopping selector with c = 1 on all rings, and moving logical states with the markers.
- `../work_C8/helix_equalizer.md`: makes the first-order logical bond independent of the lengths of the two adjacent cells (local dimension 96q + 48). Its check script was rerun by the coordinator and passes.
- `../work_C8/joint_fluid.md`: the complete first-order spectrum of length fluid plus an arbitrary translation-invariant payload, with the two sectors glued by one momentum constraint, checked against direct diagonalization to 1e-14 and giving c = 2.0000 (Casimir) for an XX payload. Rerun by the coordinator and passes. Its §4 states what it does not cover: corrections at finite coupling.
- `../../contrib/ir_control/README.md` and `../../contrib/ir_control/astra_A/RESULT.md`: infrared analysis of the frozen-cell n-flavour target (operator classification, incommensurate detuning, Sidon-set margins, the on-site Z₂ that acts as logical fermion parity).
- `../../contrib/minimal_d0/README.md` and `../../contrib/programmable_cell/README.md` for the cell machinery and local dimensions.

## The question

At fixed small coupling and on every sufficiently large ring length L, is the low-energy spectrum of the mobile-cell chain (cells of two adjacent lengths b and b + 1, a payload of n logical fermion flavours with distinct incommensurate fillings) that of n + 1 compact bosons with z = 1, a unique ground state, conformal towers, and c = n + 1? The standard is theoretical physics: standard bosonization and RG power counting, controlled Schrieffer–Wolff perturbation theory, and numerics calibrated with at least two estimators. Rigor is not required, but every step must be shown.

## Tasks

1. **Effective theory at finite coupling.** Starting from the first-order reduction of `joint_fluid.md`, list every correction that finite marker hopping, leakage out of the code and higher orders generate: residual dependence of logical terms on cell lengths, payload dependence of the marker hopping, three-cell terms, and so on. Bosonize the (n + 1)-component theory. Using the microscopic symmetries that actually survive (one-site translation, the on-site Z₂, any conserved marker or flavour numbers), classify the local operators and decide which are relevant or marginal at the Gaussian point. In particular find the resonance conditions that involve the length fluid's Fermi momentum together with the payload 2k_F's, and show (or refute) that an open set of parameters avoids all of them. Extend the Sidon-set margins of `ir_control` to include the length mode. Identify the quadratic mixing terms and the positivity condition that keeps the theory Gaussian with n + 1 gapless modes.
2. **Towers and ground-state uniqueness.** Generalize equations (3) and (4) of `joint_fluid.md` to the bosonized theory: give the zero-mode lattice including the momentum gluing. Show whether the vacuum is unique for every large L, with whatever conditions on L mod b, the parities of r and n, and the chemical potentials are needed. Predict the finite-size spectrum and the Casimir energy with distinct velocities, and say how the total c is read off when the velocities differ.
3. **What could gap the length mode.** Address commensurate locking of the cell lengths to the lattice (umklapp for the length fluid), phase separation between lengths, and ordering of the payload that could feed back on the lengths. Give conditions that avoid each.
4. **A numerical test in your folder.** Take the smallest explicit model you can (for example the joint-fluid model with generic length–payload couplings and interactions added, with n = 1 and then n = 2 incommensurate flavours) and show, with ED or DMRG and two calibrated estimators, whether c = n + 1 survives and whether the ground state is unique on rings of every length in a range that includes primes. Keep each computation to minutes on one machine.
5. **The accumulation point.** As n grows at fixed local dimension, the steps h_n accumulate somewhere in the compact space of bonds. By the c-theorem, c is locally bounded near any point whose IR is a CFT with finite c, so the accumulation point must be gapless but not conformal. Identify it in this construction (we expect the threshold where cells become infinitely long, with z = 2 and logical energy scales vanishing) and check that this is consistent.

## Rules

- Tag every claim as proved, numerical (with sizes and estimator), heuristic, or assumed. A precisely located obstruction counts as progress.
- Do not put probabilities, odds or percent-confidence estimates in any file. State evidence and status instead.
- Python (numpy and scipy) only, inside your folder. No network. Never modify files outside your folder.
- Write prose plainly, without em-dashes.

## Deliverables (in your folder)

`PROGRESS.md` (resumable state), `RESULT.md` (verdict on the question above, with the effective theory, the operator classification, the tower structure, the numerics with reproduction instructions, and ranked next steps), and code and data.
