# C14: Prove the Casimir bound Σ c_a v_a ≤ 2π‖h‖, in general or in the broadest classes you can

You are a long-running research agent. Your working directory is your own folder. If `PROGRESS.md` exists there, read it first and continue from where it stops. Update `PROGRESS.md` after every substantive step, and write `RESULT.md` at the end.

## Background

The repository is two levels above your folder (`../../`). Read `../../STATUS.md` first, then `../../contrib/refined_conjecture/README.md` (the source of this conjecture, with its bounds, its SU(d) numerics and its heuristic argument about emergent cells), its Astra subfolders `astra_fsum/` and `astra_ULS/`, and `../../REPORT.md` §3.1 to §3.4 and §7 (why a general proof of c ≤ d − 1 was blocked).

The class. H = Σ_i h_{i,i+1} on (ℂ^d)^⊗L with periodic boundary conditions, one fixed Hermitian h, one-site translation invariance, and a low-energy theory that is a set of unitary CFTs with central charges c_a and velocities v_a. Write A = Σ_a c_a v_a, the coefficient of the Casimir energy E_0(L) = e_0 L − πA/(6L) + o(1/L). Define ‖h‖ as half the spread of the spectrum of h, minimized over h → h + o⊗1 − 1⊗o (which leaves H unchanged).

The conjecture. A ≤ 2π‖h‖ for every chain in the class. It contains no d. The SU(d) chain h = P has A = 2π(1 − 1/d)‖h‖ (v = 2π/d), so the constant would be sharp. It is proved (at the physics level) for number-conserving free fermions (A ≤ π‖h‖) and for Sugawara-complete current algebras built from on-site generators (A ≤ 8πC₂(R)‖h‖/(k + h^∨)). Together with the entropy bound c ≤ 3ℓ* log d/(π(1 − ε)) it would give c² ≤ 6R log d/(1 − ε) with R = ‖h‖/T*.

## Tasks

1. **A general argument.** Look for a sum rule or variational inequality that bounds Σ c_a v_a by a local norm of h with no dependence on d. Candidates include the thermal f-sum rule for the energy current (the thermal Drude weight of a CFT is proportional to Σ c_a v_a T, while the f-sum side is a local double commutator), the equal-time or dynamical structure factor of the energy density at small q, a variational bound on the Casimir energy using twisted or boosted trial states, and Lieb–Robinson-type bounds on the energy current. For each, either complete it at the physics standard or locate exactly where it fails, and say whether the failure is technical or a sign that the conjecture is false.
2. **Broader classes.** Extend the proof beyond the two known classes as far as you can. Priorities: quadratic fermion chains without U(1) (Majorana chains, any period's worth of orbitals per site); multi-component Luttinger liquids (using vK = πD and v/K = 1/(πκ), with Drude weights bounded by f-sum rules, and whatever inequality couples the modes); Bethe-solvable chains (gl(d) and gl(m|n) graded swaps, spin-S Takhtajan–Babujian, SO(n) and Sp(2n) integrable points); and chains with an on-site symmetry whose IR is a WZW model times cosets.
3. **Stress tests.** Search for chains that come close to A = 2π‖h‖ or exceed it. Natural suspects: graded swaps on ℂ^{m|n} viewed as bosonic chains (the supersymmetric t–J model and relatives), SU(d) chains in other representations, decoupled or weakly coupled stacks (check additivity of both sides), and the emergent-cell construction of STATUS.md (whose velocities shrink with the cell length). Compute A and ‖h‖ explicitly where you can, analytically or with small exact diagonalization and finite-size fits calibrated on known cases.
4. **Consequences.** If you prove the bound in some class, state exactly what it implies for c in terms of d and the conformal window, and whether it gives any c ≤ f(d) in that class.

## Rules

- The standard is theoretical physics (standard RG and sum rules, controlled perturbation theory, numerics calibrated with at least two estimators). Show every step. Tag every claim as proved, numerical (with sizes and estimator), heuristic, or assumed. A precisely located obstruction counts as progress.
- Do not put probabilities, odds or percent-confidence estimates in any file.
- Python (numpy and scipy) only, inside your folder. No network. Never modify files outside your folder. Keep each computation to minutes on one machine.
- Write prose plainly, without em-dashes.

## Deliverables (in your folder)

`PROGRESS.md` (resumable state), `RESULT.md` (verdict, every proof with its steps, the numerical survey with reproduction instructions, and ranked next steps), and code and data.
