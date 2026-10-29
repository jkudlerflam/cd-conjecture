# C103 result

2026-10-29.

[Open] I did not obtain a general physics-standard proof of

    A = sum_a c_a v_a <= pi M log_2 d.

The bound survives the supplied solvable-family checks and the sharp free-fermion theorem. Standard infrared CFT assumptions and RG reasoning do not close any of the proposed routes developed here without an additional microscopic inequality. No crossover scale has been hidden in the conclusion.

## The best positive argument

[Proved conditional implication] Suppose physical conserved on-site currents exhaust the infrared sectors, and their stiffness matrix is

    D = (1/pi) sum_a v_a G_a.

If a positive matrix W satisfies tr(W G_a)>=c_a for each sector, the optical sum rule gives

    A <= pi tr[W(K-R)],

where K is the diamagnetic response and R>=0 is the regular optical weight. An independently proved microscopic budget

    tr[W(K-R_lower)] <= M log_2 d,   0<=R_lower<=R

would prove the desired bound with no matching length or temperature. All assumptions and the gauge qualification are in [PROOF.md](PROOF.md), Section 2. C101 supplies genuine current-complete examples of this mechanism, including balanced Hubbard and symmetric supersymmetric t-J.

[Proved] A useful improvement over discarding the regular weight follows from its positive spectral measure. If s_0 and s_1 are its zeroth and first moments, Cauchy-Schwarz gives

    tr(W R) >= 2 s_0^2/s_1,
    A <= pi [tr(W K)-2 s_0^2/s_1].

These moments can be expressed as a ground-state total-current covariance and a double commutator. This inequality is scale-free, but its right side has not been bounded universally by pi M log_2 d.

[Proved objection, using standard exact Heisenberg data] The simplest budget, obtained by discarding R, fails at Heisenberg. In spin units M=1/2, the physical spin twist has

    K=(2/3)(ln 2-1/4),   D=1/4,   W=2,
    W K=0.5908629074... > M.

The missing regular weight is R=0.0454314537.... Recovering the sharp coefficient by this route requires all of it at this saturator. More generally, neutral sectors need not carry the conserved on-site currents, so current completeness is another substantive assumption. Neither gap is resolved by ordinary spectral positivity.

## Two monotonicity shortcuts fail

[Proved] The raw Casimir estimator need not decrease toward the infrared. Start with a fermion dispersion

    epsilon(k)=-2t_1 cos k-2t_3 cos 3k,
    t_1/27<t_3<t_1/3.

Blocking three auxiliary sites gives a literal nearest-neighbor chain with d=8. In physical site units,

    A=2(t_1-3t_3)/3,
    A_ring(L)=A+[7pi^2/(810L^2)](t_1-27t_3)+O(L^-4).

It approaches A from below. The derivation uses an exact finite-ring energy on compatible vacuum sectors. The model satisfies A<=2M and therefore also the logarithmic conjecture.

[Proved] The normalized thermal energy-current variance also fails the proposed direction of monotonicity, already in the nearest-neighbor XX chain. For hopping t,

    R_J(T)=3 Var(J_E)/(pi L T^3)
          =2t-7pi^2 T^2/(20t)+O(T^4/t^3).

It increases as temperature decreases. Its high-temperature limit is zero, while its infrared limit is A=2t=2M. The energy current is conserved in this example, so current relaxation is not the cause.

[Open qualification] These results disprove identifying a weighted c-function with these observables along a general lattice flow. They do not disprove every possible velocity-weighted c-theorem, nor endpoint ordering for flows between continuum CFTs. A common fixed velocity and a relativistic c-theorem do give A_IR<=A_UV in the restricted setting, but the required lattice ultraviolet budget is still absent.

## Checks and remaining work

[Assumed inherited continuum data; proved comparisons in C101] XXZ gives Q=sin(gamma)/gamma, with Q=A/(pi M log_2 d). Heisenberg saturates. The symmetric SU(N), SO(N), conditional Sp(2r), Hubbard, t-J and critical XYZ results remain consistent. Symmetric SU(8), k=2 has Q=0.406195448..., despite violating A<=2pi M. Thus the dimension dependence cannot be discarded. Every quasi-free chain has Q<=2/(pi log_2 d), using the sharp theorem A<=2M. The full comparison table and completeness qualifications are in PROOF.md, Section 8.

[Numerical] All checks in [checks.py](checks.py) pass. Exact blocked-ring energies in both parity sectors agree with independently summed levels to 5.68e-14, and the blocked bands agree with the unfolded dispersion to 1.39e-15. Independent thermal-current integrals agree to 2.22e-16. Small XXZ-ring calculations also verify the optical-moment identity and the twist curvature. [CHECKS.md](CHECKS.md) contains the values. These are algebra and normalization checks, not new evidence for a general lattice-to-CFT theorem.

[Open] The strongest remaining target is a microscopic positive response budget that includes neutral stress-tensor sectors and bounds the surviving infrared weight. The optical formulation shows what such a budget accomplishes and why its sharp coefficient needs more than a local norm. LR plus SSA and spectral matching at frequency of order M do not supply it. The C99 and C100 matching-scale obstructions remain in place.

[Proved work record] Only work_C103 was modified. No large computation, cluster job, network shell command, delegation, commit or push was used. A few primary literature pages were consulted through the browsing tool. The parent STATE.md was left unchanged because writes outside this folder were forbidden. The local [STATE.md](STATE.md) records the handoff.
