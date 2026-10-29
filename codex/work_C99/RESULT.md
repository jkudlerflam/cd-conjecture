# C99 result

2026-10-29. Natural logarithms unless a base is displayed. Full derivations are in [PROOF.md](PROOF.md). Numerical diagnostics are in [CHECKS.md](CHECKS.md), reproduced by `python3 checks.py`.

[Open] No length-independent bound A<=C M log d is proved here. No counterexample to that bound is constructed. The investigation does establish a general obstruction to combining SIE-controlled rates, gives three explicit counterexamples to auxiliary inequalities that could otherwise appear to remove the length, and finds a positive local-quench integral that measures A.

## 1. What combining protocols can prove

[Proved conditional theorem] Suppose conformal response rates at lengths ell_j>=ell_* have the form

    r_j(t) = A f_j(t)/ell_j,      |f_j(t)|<=b n_j,

where n_j is the weighted number of cut bonds. SIE gives |r_j|<=18 n_j M log d. For any signed, absolutely integrable combination

    Q = sum_j int w_j r_j dt = A beta,
    D = sum_j n_j int |w_j| dt,

one has

    A <= 18 M log d D/|beta|,
    D/|beta| >= ell_*/b.

Thus weighting, repeating, integrating, or subtracting these rates cannot remove the smallest controlled length by the separate SIE estimates. Divergent D gives no finite bound. Cancellations could help only through an additional joint inequality. This theorem concerns the stated response class, not every conceivable entropic observable.

For C69's boundary-state plateaus, b=pi/12. The controlled result remains

    A <= 216 ell_* M log d/(pi alpha),

with relative lower matching alpha. No uniform ell_* is supplied by an infrared CFT identification.

## 2. A state-dependent correlation length fails without a restriction on the states

[Proved, microscopic auxiliary counterexample] Fix a d=4 chain consisting of an XX critical qubit and an on-site gapped spectator qubit. Its Hamiltonian and Casimir coefficient are fixed, with

    A=2J,       J<=M<=J+Delta/2.

Choose an active state with nonzero instantaneous cut rate

    r_0=(sqrt(3)/2)J log 3.

Independently prepare the spectators in a pure coherent classical-Ising state, whose local Z correlations are (tanh k)^r. Its correlation length is xi_Z=-1/log(tanh k), which diverges with k. Spectator evolution is on-site and changes no cut entropy. Therefore

    sup_psi xi_Z(psi)|dot S_psi| = infinity

for this one fixed Hamiltonian. All states have finite cut entropy. The example refutes the unrestricted state-length supremum, while satisfying the proposed Casimir bound. A length tied specifically to the excitations producing the rate, with an appropriate restriction on states, remains a possible direction. It needs a new estimate.

## 3. Local joins give a positive integral, but SIE cannot bound it

[Proved from the regulated continuum join] With common spatial smoothing epsilon,

    Delta S(t)=sum_a (c_a/6)log[1+(v_a t/epsilon)^2],
    A=(6epsilon/pi) int_0^infinity dot S(t) dt/t.

The integral restores the velocity weighting even when the individual rate peaks occur at different times. However, inserting SIE gives int dt/t, divergent at both endpoints. Restricting to t_-<t<t_+ gives only

    (6epsilon/pi) int_(t_-)^(t_+) dot S(t)dt/t
        <= (108epsilon/pi) M log d log(t_+/t_-).

A controlled lattice version needs an additional weighted-rate bound and continuum matching. The universal late-time join rate is c/(3t), which loses velocities. The smoothing parameter and the time when all sectors become visible cannot be assumed microscopic uniformly over chains.

## 4. Other proposed routes

| Route | Proven obstruction or missing input |
|:--|:--|
| Arbitrary states and local unitary control | Net entropy production is at most E(h) times interaction time. A lower bound E(h)>=kappa A is missing. |
| Inhomogeneous or Floquet driving | SIE charges the integrated norm of the actual crossing interaction. Faster clocks and repeated cycles rescale this cost too. A response lower bound in terms of the original H is additional input. |
| Distant mutual information | SIE gives a separation-independent rate bound. Existing long-distance entanglement can be transferred locally into the observed interval, so an all-state exponentially suppressed MI rate outside the LR cone is false. LR alone bounds velocities, not their central-charge-weighted sum. |
| Modular flow | The full modular unitary factorizes across its defining cut and produces no entropy there. Across another cut, a local CFT modular generator has coefficients of order R/v. No bound by the physical M for a general microscopic modular generator is supplied. |
| Euclidean evolution | An ancillary boundary-only entropy-rate bound is false, even for a generator acting strictly on one side. |

[Proved Euclidean counterexample] Take Schmidt weights p and N equal weights (1-p)/N. Correlate the first group with |0> of one left qubit and the second group with |1>. Normalized filtering by exp(-tau M Z) gives

    p'=-4M p(1-p),
    S=H_2(p)+(1-p)log N,
    S'(0)=M log N at p=1/2.

The acted-on dimension is two, its norm is M, and N is arbitrary. The filter reweights preexisting entanglement. This proves that analytic continuation of SIE is invalid; it does not refute a state-restricted Euclidean method with other estimates.

[Proved MI diagnostic] If a neighbor E of A is initially Bell-entangled with distant B, a local XX interaction gives I(A:B)=H_2(sin^2 Jt), independent of the separation. No new signal travels to B. The identity also holds in a uniform nearest-neighbor XX-Pauli chain with suitable product states on the remaining sites.

## 5. What remains open

[Proved sufficient implication] Let E(h) be the ancillary instantaneous entangling capacity across one microscopic bond. If

    E(h)>=kappa A

holds with dimension-independent kappa>0, then

    A<=(18/kappa) M log d.

This is the concrete escape from the current obstruction. SIE supplies E(h)<=18M log d; the required inequality goes in the other direction. A sufficiently strong restricted weighted-rate theorem could provide an alternative escape.

[Proved limitation of the negative result] The scalar rate, entropy-density, and LR constraints allow arbitrarily large formal A if the onset of conformal response is postponed to a sufficiently large length. Those formal response functions do not construct a finite-dimensional local Hamiltonian. The actual counterexamples above invalidate auxiliary inequalities, not the A conjecture. The campaign's slow-mode constructions likewise do not violate the velocity-weighted bound.

[Heuristic counterexample mechanism only] To violate the conjecture, a microscopic family would have to combine many infrared modes with an unsuppressed sum c_a v_a, while delaying their conformal visibility. Neither SIE nor LR rules that scenario out, but no such family is constructed here. A theorem excluding every entanglement-rate route would require substantially more than these obstructions.

## 6. Verification and scope

[Numerical] All small checks passed. A four-spin calculation verifies the nonzero active rate to 1.83e-10. The coherent Ising correlations and invariance of spectator entropy agree to 1.22e-15 and 2.29e-16. The Euclidean filter derivative agrees to 1.66e-10. A seven-spin uniform chain verifies the distant-MI identity to 1.3e-15. Numerical integration of the multisector local-join moment recovers A to 8.88e-16. These checks verify the displayed counterexamples and factors, not a general microscopic CFT theorem.

Primary sources checked through the browsing tool were [the ancillary SIE theorem](https://arxiv.org/pdf/1304.5931), [the conformal local join](https://arxiv.org/pdf/0708.3750), and [the CFT interval modular generator](https://arxiv.org/pdf/1102.0440). No network shell commands, large computation, cluster jobs, delegation, commits, or pushes were used. Only work_C99 was modified. The parent project memory was left unchanged because the task forbids writes outside this folder; local STATE.md records the handoff.
