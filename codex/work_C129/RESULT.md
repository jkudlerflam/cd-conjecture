# C129 result

2026-11-07. Derivations are in [PROOF.md](PROOF.md), numerical scope and commands in [CHECKS.md](CHECKS.md).

[Assumed, unresolved] The full theorem C_U1<infinity was not proved. No counterexample was found in the completed small-size checks. The proposed partial success (a)+(b)+(e) was also not established in its stated uniform form.

[Proved, literature audit] The key premise in TOOLBOX is incorrect at half filling. In [Kozlowski, arXiv:1508.05741v2](https://arxiv.org/html/1508.05741), Proposition 3.4 gives the all-order expansion in the finite-Fermi-boundary regime. Theorem 3.5 gives convergence, without a finite-size rate, at rho=1/2 for -1<Delta<=1. The 1/L energy term below half filling is already stated in Corollary 4.3. The paper explicitly identifies missing uniform energy control near half filling after Lemma 4.1. Thus piece (a), including its half-filled boundary, is not a corollary of the cited results.

[Proved, dependency audit] Density uniformity at fixed Delta is stated for the leading estimates in Theorem 3.2. Uniformity in Delta is not. PROOF Sections 12-13 identify the root bounds, contour widths, inverse-density denominators, resolvent norms, and size thresholds that must be controlled. Linear integral-equation continuity alone does not make the discrete-root threshold uniform. No explicit compact-region constant is extracted here. The printed complex strip also needs narrowing in the negative-Delta regime, as shown directly from the scattering-kernel poles.

## Elementary results

[Proved] Steps 1-3 hold for every L>=2, including odd lengths and the two-oriented-bond convention at L=2. With J=sqrt(J_perp^2+D^2), the DM term produces flux Ltheta, the field cancels in the selected-sector comparison, and

    e_L(n,Phi)-min_phi e_L(n,phi)
       <= (JL/2)[1-cos(pi/L)] <= Jpi^2/(4L).

The trigonometric estimate is exact as an upper bound. Its explicit fourth- and sixth-order terms are in PROOF Section 5. A minimizing flux is Lpi modulo 2pi for every sector and anisotropy, proved by a configuration-basis diamagnetic inequality. No numerical flux optimization or vacuum-compatible subsequence is needed for these reductions.

[Proved] The optimized physical M is half the range of

    J_z/4 +/- B/2,  -J_z/4 +/- J/2.

At zero field M_0=J/2 for |Delta|<=1 and M_0=J(1+|Delta|)/4 otherwise. For arbitrary field, M(B)>=M_0. Keeping this distinction repairs the normalization in the proposed reduction.

## Explicit partial bounds

[Proved] For J_z=0, every field, DM angle, and L>=2,

    L(E_(2L)-2E_L) <= [8-4sqrt(2)+pi^2/4] M(B).

The twist-optimized sector constant is 8-4sqrt(2). The free-fermion proof includes both Jordan-Wigner parities and open-shell ties.

[Proved] For J_z<=-J the unrestricted energy is

    E_L=L[J_z/4-|B|/2],

so the defect is zero. At Delta=-1 every twist-optimized sector has energy -JL/4. For J=0 the classical repetition argument gives a nonpositive defect. The defect is also zero for J_z>=-J and |B|>=J+J_z by the saturation bound. These results close parts of the original class without the stronger droplet-sector lemma.

[Proved] For Delta>=-1 and 0<=n<=L/2,

    L[F_(2L)(2n)-2F_L(n)]/M_0 <= 2pi^2 n^3/L.

Consequently n^3<=L has sector constant 2pi^2, uniform in Delta. For Delta=-1+epsilon>=-1 there is a second bound

    L[F_(2L)(2n)-2F_L(n)]/M_0 <= epsilon L^2.

Thus epsilon L^2<=A has sector constant A. Add pi^2/4 to either constant for the original physical flux and field when its minimizing sector satisfies the stated condition. These are new elementary variational estimates in this run, not Bethe-ansatz consequences. No claim of literature priority is made.

## Checks and remaining work

[Numerical] The logarithmic Bethe solver agrees with ED in all 250 sector cases for L=2 through 10 on the ten requested Delta values. Maximum energy discrepancy is 7.24e-14. Paired Bethe rings reach 20 sites, and the largest sampled ratio is 4 at L=2, n=1, Delta=1. There are no solver failures in this test.

[Numerical] Twist checks comprise 4,675 samples through L=10 and a separate 3,564-sample run through L=12, including Delta=-0.999999. All pass. The largest sampled ratio to Jpi^2/(4L) is 0.810569469, at L=2. The 36 supplied field/DM points at four lengths give 144 end-to-end checks, all passing and reproducing the printed survey values. The new variational bounds pass 72 sector checks.

[Assumed, unresolved] The gaps are (i) a uniform half-filled finite-size estimate for -1<Delta<=1, (ii) parameter-uniform control beyond the proven ferromagnetic crossover band and beyond the proven dilute band, (iii) uniform half-filled and large-Delta control on the antiferromagnetic side, and (iv) a checked extension of the published even-length asymptotics to the minimizing odd-length flux. The stronger sectorwise droplet lemma for Delta<-1 is unproved here but unnecessary for the original unrestricted target on that region. A fixed-parameter asymptotic statement and a finite-L patch do not give the requested global explicit constant.

[Work record] [run_typhon.slurm](run_typhon.slurm) prepares the requested subsampled Bethe grid through L=4096 and paired length 8192, plus Delta+1=c/L^2 families. It was not submitted. No Typhon connection was attempted. All writes stayed in work_C129, and [STATE.md](STATE.md) is the local handoff because the task prohibits modifying the project-level STATE.md. No independent referee or sub-agent review was performed.
