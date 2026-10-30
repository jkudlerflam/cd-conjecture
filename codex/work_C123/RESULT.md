# C123 result

2026-10-29.

[Proved] The bounded-cost reconnection lemma holds for XX at every fixed a>0, uniformly in even L>=4 and beta>=0. With N=max(1,ceil(beta M/(aL))) and M=1/2,

    G <= 2 log 2 + [2 log 2+(10pi/3)a]N.

This follows from C118's full-spin partition-ratio theorem and C121's inequality G<=R_L+2N log 2. It extends to independent XX stacks with log 2 replaced by log d. It does not improve the already known free-fermion Casimir bound A<=2M.

[Proved] There is no finite largest admissible spacing under the stated criterion. If the general lemma holds at one a_0, it holds at every finite a after changing the constants. Equation (3) of PROOF.md gives the conversion. The meaningful optimization is of (3 log d+b_1(a))/a. Constants uniform over arbitrarily large a are impossible already on XX.

[Proved and numerical] Derived and implemented exact finite-grid multi-pinch XX formulas, including both Jordan-Wigner sectors and finite bordered determinants for the zero-mode parity traces. The individual Gaussian histories give finite sums of determinants. The nonselective history sum is not generally one Gaussian determinant. All 51 direct Hilbert-space calibrations pass, with maximum relative error 1.96e-14. Another 32 charge-block XX comparisons pass with maximum absolute error in G of 6.75e-13.

[Numerical] Computed G itself, including the doubly pinched denominator, in 156 large-ring XX cases through L=8192, doubled length 16384, with N=1,2,4,6 and a=tau M/L between 0.125 and 16. Representative values are

| a | G/N, N=1 | G/N, N=2 | G/N, N=4 | G/N, N=6 |
| :--- | ---: | ---: | ---: | ---: |
| 0.25 | -0.531044 | -0.376830 | -0.348434 | -0.323519 |
| 1 | 0.685705 | 0.839919 | 0.876721 | 0.877622 |
| 4 | 5.586324 | 5.590035 | 5.590038 | 5.590038 |

[Numerical correction to the starting datum] At L=8192,beta=4L,N=1, the actual G is 2.372985925. C121's 0.986691567 was log(S_N/Z_(2L)), a lower bound on G. The difference comes from the two-seam pinching loss in the denominator. Both quantities remain bounded at this fixed spacing.

[Proved] For each fixed a>0 and fixed N on even XX rings,

    G=R_L-N log 2+O_(a,N)(L^(-2)),
    lim_(L->infinity)G/N=r_XX(2aN)/N-log 2.

Here r_XX is the inherited full-spin scaling function. Its large-argument behavior gives pi a/2-log 2 after taking L to infinity and then N to infinity. The constants in the error can depend on N. Uniformity in N follows from the first theorem, not from exchanging these limits.

[Numerical] Full-spectrum interacting ED covers Heisenberg and XXZ at Delta=-0.9,-0.5,0.5 on L=4,6 with doubled rings through length 12, and SU(3) on L=3,4 with doubled rings through length 8. All sectors and eigenvalues are retained. There are 192 ED cases in total, including the XX calibration, at a=0.25,1,2,4 and N=1,2,4,8. The largest diagonalized block has dimension 924. At a=1,N=8,

| Model | Smaller-ring G/N | Larger-ring G/N |
| :--- | ---: | ---: |
| XXZ, Delta=-0.9, L=4 to 6 | -0.111314 | -0.303410 |
| XXZ, Delta=-0.5, L=4 to 6 | 0.356419 | 0.341737 |
| XXZ, Delta=0.5, L=4 to 6 | 1.574825 | 1.446196 |
| Heisenberg, L=4 to 6 | 2.035604 | 1.889414 |
| SU(3), L=3 to 4 | 2.979480 | -5.071306 |

[Numerical limitation] The SU(3) L=4 ring is not vacuum-compatible. Its negative value is a finite-channel check, not evidence for convergence of the vacuum sequence. The interacting sizes do not establish a uniform bound or determine its optimal cost. A full SU(3) length-6 versus length-12 comparison and larger qubit calculations were not run. The task forbids network commands, so no Typhon connection or submission was attempted.

[Proved] For every finite-dimensional Hamiltonian, the numerator and denominator have exact positive temporal-transfer representations

    S_N=<u_S,T_S^(N-1)u_S>,
    Q_D=<u_D,T_D^(N-1)u_D>.

The maps and endpoint vectors are specified in PROOF Section 6. If lambda_j is the top eigenvalue and w_D is the squared overlap of u_D with the full top eigenspace, then

    G <= log(||u_S||_2^2/w_D)
           +(N-1)log(lambda_S/lambda_D).

This identity and bound pass 12 independent finite-dimensional checks. They identify a spectral-rate estimate and a boundary-overlap estimate that would suffice for the general theorem.

[Open precise residual] A sufficient completion is

    log(lambda_S/lambda_D) <= b_1(d),
    log(||u_S||_2^2/w_D) <= b_0(d)+b_1(d),

uniformly on the macroscopic time-slab band used by the ceiling prescription, together with the N=1 bound at shorter times. The first estimate is necessary for a uniform all-N bound at tau=aL/M. The second is sufficient but can be stronger than necessary. These transfer operators still retain the full quantum state of the rings.

[Open] Neither Araki-type locality nor a gap of order v/L supplies these estimates. C122's full operator-order obstruction applies, and a gap alone leaves uncontrolled spectral multiplicity. No general interacting reconnection theorem or new universal bound on A was proved. The exact residual is stated rather than assumed.

[Work record] Derivations are in [PROOF.md](PROOF.md), numerical details and reproduction in [CHECKS.md](CHECKS.md), and the handoff in [STATE.md](STATE.md). All files and computations remained in work_C123. No large local calculation, cluster job, external write, or independent agent review was performed.
