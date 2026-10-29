# C101 result

The bound survives every completed test. The work extends C35 with a stiffness proof for balanced Hubbard at all fillings, a stiffness proof for flavor-symmetric supersymmetric t-J at all fillings in its conformal phase, exact field-dependent bond widths, a conditional Sp(2r) comparison using its actual multiple velocities, and a small adversarial ladder search.

It does not prove the requested statement for all integrable chains. Arbitrary finite-density SU(N|M) phases and general several-sea dressed-charge matrices remain unresolved. The Sp(2r) continuum assignment is itself a published analytical/numerical identification, not a rigorous all-rank theorem. Heavy ED did not run because SSH to Typhon timed out. These limitations are part of the result.

## What is proved, and what is assumed

Every M below is gauge optimized over arbitrary Hermitian on-site coboundaries. The proof is reflection averaging, with a parity reduction for Hubbard. Continuum-dependent results assume the stated complete Bethe/Luttinger identification, not merely the presence of a few known modes.

| Family and regime | Exact M | Result for Q=A/(pi M log_2 d) | Status |
| :-- | :-- | :-- | :-- |
| SU(N) fundamental, h=P | 1 | <=1, equality N=2 | Inherited C35 proof conditional on Bethe/CFT input |
| SU(N), symmetric rank k | H_k in the 2H-channel normalization | <=1; every k>=2 has Q<3/4 | Inherited C35 proof |
| Spin S TB, h=sum H_j Pi_j | H_{2S}/2 | Same N=2 symmetric-family inequality | Inherited C35 proof |
| SO(N) vector, N>=3 | N/(N-2) | 1/log_2 N | Inherited C35 proof |
| Zero-field XXZ and critical XYZ equivalents | XXZ M=1; general XYZ M=(two largest absolute J's)/4 in spin units | sin(gamma)/gamma<=1 | Proved from the known critical-locus identification |
| Sp(2r), r>=2 | 1 | [2(r-1)/r+1/(r+1)]/log_2(2r)<1 | Conditional on Martins's multivelocity sector assignment |
| Sp(2) | 1/2 | 1 | Rescaled Heisenberg |
| Balanced repulsive Hubbard, all metallic fillings | Half the range of the explicit 16-level list in PROOF.md; always >=2t | <=3/(2pi)=0.477465 | New stiffness proof, using K_c>=1/2 and complete charge/spin sectors |
| Balanced attractive Hubbard | Same exact spectral formula | <=1/(2pi)=0.159155 | Conditional on spin gap and standard K_c>=1 input |
| Supersymmetric t-J, q=2 equal flavors | Explicit five-value max/min formula, always >=t | <=3/(2 log_2 3)=0.946395 | New stiffness proof plus exact rational four-site certificate |
| Supersymmetric t-J, q>=3 equal flavors | Same formula | <=2/log_2(q+1)<=1 | New stiffness proof using the known exponent range |
| Arbitrary grading and diagonal fields | Exact spectrum given in PROOF.md Section 7 | No general A comparison proved | Local optimization solved; continuum problem remains |

The t-J upper bounds are not predictions for the actual maxima. In particular, the q=3 upper bound of one is not an identified saturator. The only established saturation in these families is Heisenberg and its equivalent normalizations or unitary images. XXZ approaches it continuously. The discrete symmetric-SU(N) runner-up is fundamental SU(3), Q=0.8412396714.

## The new mechanism

The positive optical sum rule bounds the stiffness matrix by its diamagnetic expectation. For balanced Hubbard this gives

    2v_c K_c/pi <= -<T>/L,
    2v_s/pi     <= -<T>/L.

The free-hopping variational bound is -<T>/L<=(4t/pi)sin(pi n/2). Together with K_c>=1/2, it gives A<=6t sin(pi n/2). The exact Hubbard bond has M>=2t, already more than sufficient for d=4. This proof does not require an elementary formula for either velocity at generic filling.

For SU(q) t-J, use a positive weighting of the charge and spin stiffnesses. Kawakami's electron exponent range implies K_c>=1/q in the normalization of PROOF.md. The weighted diamagnetic bond has largest eigenvalue 2t, proving A<=2pi t for q>=3. At q=2, SU(2) averaging and an exact four-site operator certificate improve the estimate to A<=(3/2)pi t. The certificate uses rational LDL factorizations of 15 blocks, total dimension 81. No floating-point assertion enters that inequality.

The stiffness mechanism relies on complete current-carrying Luttinger sectors and a lower bound on the relevant Luttinger parameter. It is not a universal proof and cannot count an extra neutral sector for free.

## Completed adversarial tests

The 66-point repulsive Hubbard quadrature scan found maximum sampled Q=0.3168709574 at U/t=0.5, n=0.9904793267. There v_c=2.0712087014 and v_s=1.9184787897. In particular, v_c can exceed the free single-particle maximum 2t; a proof based on v_c<=2t would be false. The largest final change under quadrature refinement was 1.01e-6 in the checked velocities/density/chemical potential. Four half-filled spin-velocity checks agree with the independent Bessel formula to 4.2e-15. The scan is numerical; the all-filling argument is the stiffness proof.

The d=4 ladder/composite screen used

    h=a P_spin+b P_orbital+c P_spin P_orbital

at 14 coupling points and L=4,6,8, with sector dimensions at most 4900. The exact M is half the range of a s+b t+c s t, s,t=+/-1. No two-size Casimir diagnostic exceeded one. The largest fitted Q was 0.8745965622 at SU(4), whose exact Q is 0.75. Decoupled Heisenberg legs gave 0.5656260606 against the exact 0.5. The largest eigenpair residual was 1.15e-11. For unfamiliar couplings these are fixed-magnetization sector screens, not certified global ground states or thermodynamic CFT determinations.

Products of m identical Heisenberg chains have d=2^m, M=m/2, A=m pi/2, and Q=1/m. More generally, products of verified factors preserve the conjecture. Increasing the number of independent fast species this way cannot violate it. Composite Bethe branches, unselected graded phases with quadratic modes, and longer-range integrable densities cannot be counted as additional admissible nearest-neighbor CFT sectors.

## Should the dimension dependence change?

No completed test requires changing pi M log_2 d. Heisenberg forces the coefficient at d=2. The symmetric SU(N), k=N sequence already established in C35 forces any universal replacement to grow at least as log d/log log d along that sequence. It does not force a full log d, so the asymptotically optimal dimension dependence remains open. A dimension-independent replacement is excluded.

## Files and remaining work

[PROOF.md](PROOF.md) contains the spectra, gauge proofs, current normalizations, stiffness derivations, assumptions, and primary references. [NUMERICS.md](NUMERICS.md) records the finite-size conventions and reproducibility commands. `bond_checks.py`, `tj_certificate.py`, `hubbard_quadrature.py`, and `adversarial_ed.py` produced the corresponding JSON records. `typhon_ed.slurm` prepares the larger ladder screen and refuses login-node execution through the Python guard. It was not submitted.

The unresolved work is (i) an admissible-phase classification and velocity/stiffness-matrix bound for general Hermitian SU(N|M) chains with flavor fields, (ii) a uniform justification of the full Sp(2r) Casimir content, and (iii) larger commensurate ring checks for unfamiliar ladder points. No general proof or counterexample is claimed. Only this task folder was modified; the project-level STATE.md was left untouched under the explicit confinement rule.
