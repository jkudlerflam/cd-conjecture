# C82 result

[Proved construction; conditional infrared conclusion] F5 has a local
repair with d=17. Replace the single vacancy by V_0,V_1, carrying the
grammar state of the preceding material letter. Penalize incompatible
vacancy-letter and vacancy-vacancy pairs, and let vacancy hopping update
the color. The occupied 15-letter bond is unchanged. The complete
289 by 289 bond is in [bond17.py](bond17.py).

[Proved] Every illegal material pair now forces a local penalty, even
when arbitrarily many vacancies separate its letters. On the legal
space the coloring is unique, and the hopping reduces to ordinary
vacancy hopping with amplitude -t. At finite K the full Hamiltonian
still mixes legal and illegal words. Their suppression, not exact
finite-K invariance, is the repaired claim.

[Numerical] On the referee's two n=3 cells plus one vacancy at J=t, the
uncolored illegal material weight tends to 0.0256. With the repair it is
1.07551e-4, 1.08213e-6, and 1.08279e-8 at K=100,1000,10000. The repaired
closure has dimension 4608. A literal full-bond test with nonzero port
and conversion, dimension 1344, also gives K^-2 suppression. The largest
eigenpair residual in these tests is below 2.2e-11. These are leakage
tests, not central-charge measurements.

## The two proof gaps

[Proved, using C75's all-cell inequality] F3 is closed for the frozen
chain at nonzero lambda. Apply the parent's boundary-energy lemma to
complete finite-K cells, including their illegal components. With
B=sqrt(2)|lambda| and gamma=delta_selector-1/(K-4), a compatible ring
with k wrong-length cells satisfies

    E_wrong-E_selected >= k(gamma-5B)>0.

Marker-free rings are excluded by C75's separate density bound. This
selects the full marker sectors. It does not treat bad-pair labels as
conserved or claim a finite-K ground state with zero illegal weight.

[Proved, finite-order local normal form] F4 has an additional
Schrieffer-Wolff layer for the commuting integer penalty. Its denominators
are kappa times nonzero integers, with kappa=K in the frozen chain and
kappa=JK in the repaired lift. At order r,

    U_r H U_r^dagger = kappa N_bad+D_r+R_r,
    [D_r,N_bad]=0,
    ||R_r||_local,exponential <= C_r nu (nu/kappa)^r.

Here nu is the unpenalized interaction strength per site and C_r is
independent of chain length. [PROOF.md](PROOF.md), Section 4, gives the
local decomposition, both first generators, the effective second-order
term, and the connected-support remainder estimate. It does not assert
an exact all-order diagonalization or a globally separated many-body
band. A 96-dimensional matrix check verifies the second-order signs
and remainder scaling.

## Dimension and candidate central charge

[Proved, new computer-assisted certificate] The extra vacancy color
costs one state relative to C81. C75's original n=15 payload would give
the conditional count c=16=d-1 at d=17. To recover a violating candidate,
this folder also certifies a new n=16 member using attributed copies of
C75's integer/rational routines.

| Finite input | Result | Tag |
| :-- | :-- | :-- |
| Local alphabet | 15 material letters and 2 vacancies, d=17 | Proved |
| Material cell | n=16 data sites, b=17 | Definition |
| Fermi roots | 32 simple roots, exhaustive | Proved, computer-assisted |
| Minimum root separation | greater than 0.0145298937 | Proved, computer-assisted |
| Pair-sum separation | greater than 5.92937e-7, all 528 pairs tested | Proved, computer-assisted |
| Finite-g,K error budget | below one eighth of 2^-35 | Proved, computer-assisted |
| Frozen nonzero-port selector | gamma>2^-18, gamma>5sqrt(2)lambda | Proved |
| Payload infrared theory | 16 bosons after matching and stability | Assumed |
| Vacancy reference | one charge boson at half occupancy | Assumed, inherited phase-model result |
| Finite-J matching and dressed defect comparison | still required | Assumed |
| Total c | c_candidate=17>d-1=16 | Heuristic, conditional on the preceding infrared inputs |

[Definition] The new certificate uses D=1000, Delta=4-2cos(pi/18),
phi_A=pi/32, phi_B=0, psi_A=pi, psi_B=pi/2, and

    beta=1/d_17,   d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1),
    t_d=beta 2^-92,    g=2^-288,    K>=2^1280,
    0<lambda<=2^-704,  field=lambda^2 alpha^4 beta^2 g^2 2^-149,
    alpha=sqrt(2/18)sin(pi/18).

The b=17 midpoint selector has mu_sel=1.912016042804466... . The raw
root routine's finite_g_power=-264 is an intermediate sufficient
threshold. The final choice is -288, recorded in budget_n16.json.

[Proved scale requirement] In the mobile model the penalty scale is JK,
while the color-changing hopping has strength t. The new virtual terms
are O(nu^2/(JK)). Choose the material parameters, lambda and positive J/t,
then increase K until JK>>nu and C_2 nu^2/(JK)<J tau, with tau the desired
material accuracy, including the smallest logical scale. Increasing K
improves the finite certificate. The lower bound 2^1280 alone does not
certify a mobile matching window.

[Assumed, remaining] A positive infrared tolerance tau, volume-uniform
finite-J matching, dressed material-defect stability, and multichannel
Gaussian stability are not supplied by the leakage tests or the band
certificate. The d=17, c=17 statement is a conditional candidate, not
an established microscopic counterexample. No d=16 repair or universal
minimality theorem is claimed.

## Files and reproduction

Run from this folder with Python 3, NumPy and SciPy. All scripts write
only here. The computations use small matrices and one BLAS thread.

    python3 check_vacancy.py
    python3 check_full.py
    python3 check_parent.py
    python3 check_sw.py
    python3 certify_n16.py
    python3 budget_n16.py

The outputs are vacancy_checks.json, full_checks.json, parent_checks.json, sw_checks.json,
root_certificate_n16.json, root_certificate_exact_n16.json, and
budget_n16.json, together with the intermediate integer matrices.
The scalar routines copied from C75 are identified in their first lines.
The new n=16 certificate does not overwrite or rerun C75's n=15 outputs.

All work, including the local memory in STATE.md, is confined to work_C82.
C81 and project-level memory were read but not modified. No network
commands, cluster jobs, or subagents were used.
