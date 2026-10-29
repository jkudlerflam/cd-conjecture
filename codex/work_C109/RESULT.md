# C109 result

[Proved, scope of result] No general bound A_neutral<=C M f(d) for an arbitrary collection of neutral or coset sectors is established here. The three proposed routes have distinct missing inputs, stated below. No counterexample to the desired bound was found or claimed.

## Positive bounds

[Proved, conditional on dynamical CFT matching] A nearest-neighbor Lieb-Robinson estimate with the optimized bond norm gives v_a<=4 e M. Consequently an individual Ising factor obeys cv<=2 e M, an individual unitary Virasoro minimal-model factor obeys cv<4 e M, and an individual standard Z_N parafermion factor obeys cv<8 e M. These hold for interacting lattice realizations as well as solvable ones. They require nonzero local real-time stress response in the scaling limit, not merely a Casimir fit.

[Proved, limitation] With r_min minimal-model factors and r_pf other standard parafermion factors, the corresponding sum is at most 4 e M(r_min+2r_pf). This is crossover-free but is not a bound in terms of d without a multiplicity theorem. No such theorem is supplied. The inherited stronger result A<=2M still applies to whole quasi-free chains.

## Discrete twists

[Proved] A one-bond on-site symmetry seam gives |delta E_g(L)|<=2M. Under the topological-defect CFT identification,

    delta E_g(L) = (2 pi/L) sum_a v_a x_a(g)+o(L^-1).

The norm estimate therefore leaves a factor L. What is missing is a bound on limsup L delta E_g(L), as well as coverage of all sectors by the available twists. A discrete symmetry does not justify distributing its twist as a continuous conserved-charge flux.

[Proved, conditional on the standard seam identification] For Z_N parafermions, c_N=2(N-1)/(N+2) and the disorder dimensions are x_k=k(N-k)/[N(N+2)]. The strongest twist satisfies c_N<=8x_k, while averaging all twists gives c_N<12 x_average. Thus an averaged seam budget

    limsup L average_g delta E_g(L) <= K(d) M

would imply A_neutral<=(6/pi)K(d)M for products with full cyclic symmetry images and product twisted vacua. The factor 12 is independent of parafermion orders and multiplicity. The missing lattice budget remains unproved.

[Proved] A universal twist exists after replication: cyclically joining n copies gives one ring of length nL. Exactly,

    delta E_cycle(L)=E_0(nL)-nE_0(L).

It detects every neutral sector, but for n=2 its limiting coefficient is pi A/4. It is precisely the C98 ring-gluing problem, not an independent bound.

## Vacuum currents

[Proved within the assumed CFT] The equal-space vacuum energy-current tail is

    <j(x)j(0)>_stress = [sum_a c_a v_a^4]/(4 pi^2 x^4).

The energy-density tail instead contains sum c_a v_a^2. The lattice bounds ||h||<=M and ||j||<=2M^2 control the correlation values, not the coefficients of their asymptotic tails. PROOF.md gives a positive-definite covariance kernel with fixed local variance and an arbitrarily large x^-4 coefficient, demonstrating this inference gap.

[Proved counterexample to a velocity floor] The d=2 critical anisotropic XY family has

    M_gamma=sqrt(1+gamma^2),
    c=1/2,  v_gamma=2gamma,  0<gamma<=1.

It has a neutral Ising infrared for every positive gamma, yet v_gamma/M_gamma tends to zero. Thus a universal positive lower bound on v/M cannot repair the proposed moment argument. This example satisfies the desired upper bound on A with increasing slack.

[Proved accounting] Even a bound on sum c v^2 or sum c v^4 plus an upper speed limit does not control sum c v without an additional constraint on multiplicity at small speeds. The exact relation is A=integral_0^(v_LR) N(u) du, where N(u)=sum c_a 1_(v_a>=u). Higher-moment bounds give nonintegrable upper envelopes at u=0.

## Charged extensions

[Proved conditional reduction] A lift to a chain with full microscopic WZW currents would suffice if it simultaneously obeyed d'<=g(d), M'<=K(d)M, and A'_WZW>=eta(d)A_neutral with eta(d)>0. C104 would then give

    A_neutral <= [4 pi K(d)/(3 eta(d))] g(d)[g(d)^2-1] M.

[Proved limitations] A c<1 minimal model cannot acquire a continuous current while keeping the same CFT. Adding a decoupled U(1) spectator does not charge it. The continuum parafermion-plus-boson WZW extension requires appropriate locality and matched velocities, and does not supply the stated lattice controls. Bounded perturbation cost alone gives no comparison between the two Casimir coefficients. A lift only to U(1) also retains C104's stiffness-weight gap.

[Numerical] All small checks passed, with matrices at most 64 by 64. They verify the XY spectrum and gauge optimization identities, Ising twist normalization, the full replica-seam spectrum, parafermion dimension arithmetic and the positive covariance example. See CHECKS.md and run `python3 checks.py` to reproduce them. No heavy numerics or cluster jobs were run.

[Proved, deliverables] PROOF.md contains the derivations, hypotheses and primary-source links. PROGRESS.md and local STATE.md record the outcome. Only files inside work_C109 were written; the parent project memory and REPORT.md were not modified.
