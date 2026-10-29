# C90 result

[Proved, main conclusion] C48's abstract massive-matching theorem, C58's rooting and ensemble lemmas, and C51's Gaussian central-charge identification transfer to the vacancy constructions with the hypotheses stated in [PROOF.md](PROOF.md). The complete d0=4464 microscopic application does not transfer merely by replacing the alphabet or the length fluid. The supplied notes leave finite-J vacancy initialization open. A thermal defect comparison and Gaussian completeness must also be distinguished from the narrower ground-energy and RG assertions.

## What transfers

| Ingredient | Transferred statement | Application boundary |
| :-- | :-- | :-- |
| C48 | Spectral resummation of the entire massive auxiliary, all-degree retarded kernels, leakage-linear error, massive thermal winding and controlled source derivatives | Requires a verified mobile representation, dressed gap, strong coefficient bounds and smallness. These are not established for C88 by its frozen certificate. |
| C58 | Fixed-density root pinning, bounded seams, pressure zero, canonical duality, and preservation of uniform activated or o(T^2) pressure errors | The wrapper's invariant mismatch flags and all-string onsite proof do not transfer to the colored bond. A mobile defect trace comparison is needed. |
| C51 | c_R=c_L=n+1 for a complete positive compact Gaussian theory, arbitrary velocity ratios, and the thermal scalar summation | Gaussian completeness remains necessary for equality of the full Virasoro central charge. TP remains necessary for the lattice thermal statement. |

[Proved, routine transfers supplied] PROOF sections 2, 3, 5, 6, 7, 9 and 10 give the coloring argument, charge normalization, conditional C48/source transfer, root-pinning and ensemble proofs, a sufficient defect-activity summation, the scalar-shell summation, and the Gaussian central-term calculation. The existing material certificates and finite-order grammar elimination are imported with attribution rather than relabeled as new proofs.

## What changes in the mobile chains

[Proved] The charge variable is an unbounded gap occupation with norm-one shifts E^+/-, not C58's binary length bit. At nonzero occupancy, the colored legal basis has one coloring for each squeezed material word and occupancy pattern. Colors add no independent legal oscillator. The all-vacancy sectors are separate competitors.

[Proved reference thermodynamics; conditional Gaussian dictionary] With rho=N_occ/L,

    e_particle(rho)=-(2t/pi)sin(pi rho),
    v_gap=2t rho sin(pi rho), K_gap=rho^-2.

At half occupancy K_gap=4, while K_particle=1 and the physical charge speed is 2t. These are different coordinate and compactification conventions. Treating the phase-model charge as a small perturbation of a K=1 fermion in the same gap variables is unjustified. The microscopic contact/harmonic dictionary and its finite-J continuation belong to initialization.

[Imported, proved data; conditional counts] The constructions have the following selected lengths and coordinate factors.

| Chain | n | Material length b | Physical length per cell at rho=1/2 | Conditional c |
| :-- | --: | --: | --: | --: |
| C88, d=13, two colors | 13 | 15 | 30 | 14 |
| C83, d=14, two colors | 13 | 14 | 28 | 14 |
| C85, d=14, three colors | 13 | 15 | 30 | 14 |
| Original d=24, one vacancy | 23 | 24 | 48 | 24 |

[Proved, finite-K qualification] C88's complete leading grammar correction is

    -(t^2 C+J^2 A+2J^2 lambda^2 Z)/(J K_pen).

The t^2/(J K_pen) term is a bit-dependent contact field, with no lambda^2 suppression. Its size relative to the material scale J grows as t^2/(J^2 K_pen). Thus taking J smaller at fixed K_pen does not automatically improve the colored chain. PROOF (8) gives sufficient finite-penalty inequalities. A joint schedule must also control the K_pen dependence of the separate matching constants. A fixed local remainder cannot be discarded as o(T^2).

[Proved, two representation issues] Contact factors can vanish on retained gap configurations, so the frozen material gap does not imply a uniform onsite massive gap in the mobile representation. Also, Q=2N_payload on the small-alphabet code, but odd-Q auxiliary states make exp(i pi Q/2) have eigenvalues +/-i. The audited separate-Jordan-Wigner lemma cannot use it as a full-space fermion parity. The original d=24 frozen matching verifies a different grading; this does not solve mobile matching for d=13 or d=14.

[Proved, locality qualification] A fixed squeezed interval can contain arbitrarily many physical vacancies. The wrapper's bounded-code-length proof of physical source locality therefore does not transfer. Bounded seams suffice for pressure, while a physical field identification also needs control of the fluctuating coordinate map.

## Conditional d=13 theorem

[Assumed, complete input list] Fix one C88 member with n=13 and all microscopic parameters independent of volume and temperature. The theorem requires the following named inputs, stated quantitatively in PROOF.

1. I, source-uniform mobile initialization. This includes finite-J low-frequency vacancy matching, a justified massive/grammar elimination and grading representation, strong all-degree endpoint bounds, the physical harmonic dictionary, exact scalar normalization and source matching, and a jointly admissible finite parameter schedule.
2. D_E and D_T, dressed physical defect comparison. D_E selects the good ground-state phase, including C88's pin defects, wrong lengths and marker-free/all-vacancy competitors. D_T supplies a uniform activated trace or pressure comparison. Long wrong cells need length-weighted costs. C88's ground-energy S_pin alone does not imply D_T.
3. N_mob, the multichannel comparison, counterterm and source theorem for the initialized payload plus K_gap=4 charge action. This must preserve the positive Gaussian neighborhood and nonresonance margins and give the required density response. C79's source-free moment estimate does not close its comparison and source parts.
4. TP_mob, the uniform ground-subtracted scalar estimate for the actual thermal histories and a trace-normalized interacting reference. A sufficient shell bound is C q_h^2 eta_h min(1,(T/q_h)^p), with p>2 and eta_h tending uniformly to zero.
5. G_mob, C51's complete compact Gaussian identification, including charged sectors and exhaustion of the physical infrared theory. This can be included as an explicit conclusion of a stronger N_mob, but is not implied by selected two-point correlations.

[Proved conditional theorem] Under I, D_E, D_T, N_mob and TP_mob, the full physical chain, in the tuned phase or its compatible fixed half-density ensemble, obeys

    f(T)=f(0)-(pi/12) sum_(a=1)^14
                    [1/u_(a,R)+1/u_(a,L)] T^2+o(T^2).

The u's are dressed physical speeds. At half occupancy u=30 v for speeds v in C88 cell coordinates. With G_mob as well,

    c_R=c_L=14, c_cd=14>12=d-1.

The proof retains the exact scalar, performs the thermodynamic limit before T tends to zero, and does not differentiate an uncontrolled thermal remainder. No common Lorentz speed is required in this collection-of-sectors convention.

[Proved, interpretation of the expected three inputs] The shorthand "N, TP, dressed-defect comparison" is valid if N is explicitly expanded to include I and G_mob, and the defect comparison includes both D_E and D_T. With the narrower existing lemmas it omits substantive requirements. This audit does not close those requirements or certify one numerical finite-J infrared point.

[Proved, limits] No o(1/L) spectrum, finite-ring tower multiplicities, interval-entanglement asymptotic, or unconditional c-d violation follows. C88's finite material certificate is unchanged and is not contradicted by these application gaps.

## Verification and files

[Numerical] `python3 -B checks.py` passes. The 54 number-sector pinning tests have compression-identity residual 1.12e-16 or smaller, with largest matrix dimension 64. Independent free-fermion quadrature approaches the half-density charge thermal coefficient pi/(12t), with relative error 1.71e-5 at T=0.01 and t=1.3. Gap-energy finite differences recover K_gap=4 to 1.01e-9 at half density and check K_gap=rho^-2 at five densities. Synthetic penalty-budget checks and the C88 physical momentum-margin arithmetic pass. The synthetic constants are not a determination of C88's matching constants.

[Record] [checks.py](checks.py) and [checks.json](checks.json) contain these controls. [PROOF.md](PROOF.md) contains the derivations and dependency ledger; [PROGRESS.md](PROGRESS.md) and [STATE.md](STATE.md) record the handoff. The C88 certificate and microscopic diagnostics were not rerun. No microscopic central charge was computed. No network commands, subagents or writes outside work_C90 were used.
