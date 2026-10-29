# C115 result. Counting fast sectors

2026-10-29.

[Open] The general bound A<=C M f(d) remains unproved. No counterexample to it was found. Finite-size counting needs a velocity-resolved infrared estimate which does not follow from the available locality bounds.

## What is established

[Proved from finite-level CFT matching] The proposed tower inclusion has its inequality reversed. A state of conformal dimension x in a sector of speed v costs 2 pi v x/L. The window 2 pi u x/L includes the tower only through dimension u x/v. Fast sectors cost more energy, while slow sectors populate the window more densely.

[Proved counterexample to the unfiltered trace bound] Two independent XX chains with hoppings 1/(1+epsilon) and epsilon/(1+epsilon) have

    d=4,   M=1,
    v_f=2/(1+epsilon),   v_s=2epsilon/(1+epsilon),
    N(1)=1,   A=2=2M.

For any fixed tau>0, slow-sector oscillator states give

    liminf_L Tr exp[-tau L(H_L-E_0(L))]
       >= [1-exp(-2 pi tau v_s)]^(-1).

This diverges as epsilon tends to zero at fixed d, M and fast multiplicity. The number of eigenstates in a fixed positive O(1/L) window diverges too. The proof takes L to infinity at fixed epsilon, and uses only finite collections of conformal levels before summing. Thus the obstruction survives working strictly in the infrared. It does not rely on a uniform crossover assumption.

[Proved] A local gate's operator Schmidt rank is at most d^2, but a spatial transfer matrix for a quantum Gibbs operator crosses all imaginary-time layers. Its uncompressed dimension grows with the Trotter number. A generic fixed-d^2 quantum transfer matrix does not follow from nearest-neighbor locality. Also, filtering on a many-body ratio E/|P| does not select fast sectors: counterpropagating slow excitations can have arbitrarily large such ratios.

## A precise sufficient lemma

[Definition] Put V=4eM and C_j=sum c_a over V 2^(-j-1)<v_a<=V 2^(-j).

[Assumed missing lemma] If one proves, uniformly over all chains and all j,

    C_j <= K(d) 2^(j/2),

with an explicit function K(d), then no symmetry or crossover assumption is needed beyond the stated infrared matching.

[Proved conditional conclusion]

    A <= sum_j V 2^(-j) C_j
      <= 4e(2+sqrt(2)) M K(d).

The coefficient is about 37.1232. A different sufficient estimate C_j<=K(d)(j+1) gives A<=16e M K(d) and a squared-logarithmic tail for N(u). The function K(d) has not been bounded for general interacting chains.

[Proved spectral lower bound, under full modular state-content assumptions] For a complete modular-invariant CFT factor in one velocity window, let Q_j=Tr exp[-pi sum(h+hbar)] after the infrared limit. S-invariance and the dual vacuum contribution give

    log Q_j >= pi C_j/4.

An upper bound log Q_j<=B(d)2^(j/2) would therefore imply

    A <= [16e(2+sqrt(2))/pi] M B(d).

This counts central charge rather than one vacuum tower per sector. The missing input is the lattice upper bound and, for projected or nonfactorized sectors, a spectral selection with the required lower bound. A trace over vacuum descendants alone is insufficient. Full dimension traces also see compactification radii, so their proposed upper bound is stronger than the central-charge bound. PROOF Section 5 states these qualifications and the unequal-chiral-speed version explicitly.

## Dynamical route

[Proved limitation] The Lieb-Robinson estimate bounds propagation speed but supplies no lower bound on channel bandwidth or independent microscopic access. Even granting a packet protocol of width ell, small incremental entangling only gives C_w<=const (M/w) ell log d for a window w<v<=2w. Constant ell already leaves a nonintegrable 1/w envelope if each window is bounded separately.

[Proved obstruction to that packet construction] C106's fixed-speed pockets require arbitrarily long packets to remain within the individual linear infrared branches. A simultaneous protocol with uniformly bounded packet or preparation scale would suffice, but it is an extra lemma, not a consequence of speed. The same missing scale identified by C69 and C99 remains.

## Tests

[Inherited theorem; proved consequence] Every nearest-neighbor quasi-free chain with d=2^m and a pure conformal infrared obeys N(u)<=min(m,2M/u). It therefore satisfies the sufficient window estimate with K(d)=log_2 d, and its stronger known bound is A<=2M.

[Proved] XX stacks have M=sum t_a, v_a=2t_a and A=2M. Their counting functions and full compact-boson dimension traces satisfy the stated sufficient estimates. Equal stacks also require any proposed power-tail prefactor to obey g(2^m)>=2^(1-delta)m^delta.

[Proved conditional test] If n payload modes have maximum normalized speed at most K exp(-kappa n), then their fast count is at most max(0,log(K M/u)/kappa), plus at most one length mode. It stays bounded for u of order M even as c grows. A bound on the smallest speed alone would not establish this.

[Inherited numerical test with proved gauge normalization] The campaign's spatial-cell table records extremely small largest payload speeds through n=23. The actual vacuum-token hop changes both site letters, so its amplitude -t is unaffected by every bond coboundary and proves M>=t. Restoring the physical cell length b=n+1 bounds v_max/M by b v_max,cell/t, whose reported value at n=23 is approximately 24 times 10^(-205.2) lambda^2 g^2 alpha_n^4. The tabulated frozen effective theories have no payload sector with speed of order M. This retains the campaign's perturbative and infrared qualifications; the printed numerical values are not certified enclosures, and no new all-n theorem is claimed. Comparable slow velocities also warn against assuming a uniform number of sectors per dyadic window.

[Numerical] All small checks pass. The largest matrix is 16 by 16. See CHECKS.md and checks.py. No heavy numerics or cluster jobs were needed.

[Proved work record] PROOF.md contains the derivations and source links. PROGRESS.md and STATE.md give the handoff. Only work_C115 was modified. The project STATE.md and REPORT.md remain unchanged under the folder-only rule.
