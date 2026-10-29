# C108 result

[Proved, lattice theorem] The chemical-potential range is bounded by the bond scale. If delta_- and delta_+ are the gaps from the minimum and maximum eigenvalues of the on-site charge to their nearest distinct neighbors, then every interior thermodynamic subgradient obeys

    -4M/delta_- <= mu(n) <= 4M/delta_+.

For an integer-valued charge this gives a range of length at most 8M. The proof replaces every nonmaximal charge state by a maximal charge state with a local trace-preserving channel. The energy cost is at most 4M per replaced site. It works with degenerate charge eigenspaces and avoids dividing by a small matrix element of a charge-raising operator. PROOF.md Sections 2 and 3 give the derivation.

[Proved, averaged bound] The distributional curvature nu=e'' is a positive measure and

    nu((q_-,q_+)) <= 4M(1/delta_- + 1/delta_+) <= 8M

in integer normalization. Its atoms include Mott chemical-potential jumps. On smooth one-mode Luttinger density branches J, the supplied stiffness bound gives

    integral_J v(n)^2 dn
       <= 4 sigma_q^2 M^2(1/delta_- + 1/delta_+).

The mean curvature over an interval of length ell is bounded by the curvature mass divided by ell. No uniform pointwise bound follows from this averaging. Under charge rescaling the chemical potential scales inversely, so a bound proportional to the charge spread alone needs a fixed normalization. The endpoint-gap formula has the correct covariance.

[Proved, conditional pointwise theorem] A charged Gaussian mode appearing with nonzero weight in the real-time scaling limit of a local conserved density satisfies

    v_a <= 6 exp(1) M.

An explicit one-dimensional Lieb-Robinson path count gives this constant. Its commutator cone passes to any polynomially normalized distributional scaling limit. A Gaussian density has a nonzero delta-derivative commutator on its propagation ray, which cannot lie outside that cone. The assumption is local real-time Gaussian matching, stronger than a static stiffness dictionary or a Casimir fit alone. The argument needs no uniform crossover length.

[Proved, consequence under the scalar response dictionary] For a non-drifting Luttinger sector,

    D e'' = v^2 <= 36 exp(2) M^2,
    D <= 36 exp(2) M^2/e''    when e''>0.

Thus large smooth curvature does force small stiffness under these infrared assumptions. This is a locality proof, not a joint twist-density variational proof. The separate static bounds permit abstract smooth convex spikes with arbitrarily large D e'', so another input is necessary. No fixed-d nearest-neighbor example of unbounded smooth e''/M for a fixed integer charge was found, and no bound on that curvature alone was proved.

[Proved from the stated Mott-edge assumptions] The ordinary Mott doping edge does not produce large e''. If v tends to zero and K tends to a finite positive value, e''=pi v/K tends to zero, while compressibility diverges. At the plateau itself e(n) has a cusp, not a finite Luttinger curvature. The exact d=3 atomic example in PROOF.md has curvature measure 2U delta_1 and D=0. One cannot use the Luttinger product identity on that atom. The conventional edge behavior is also recorded in [Cazalilla et al., Section VI.A.3](https://arxiv.org/html/1101.5337v1#S6.SS1.SSS3).

[Proved, conditional multi-mode theorem] For non-drifting Gaussian modes,

    D=Z^T V Z/pi,  kappa=Z^T V^(-1) Z/pi.

If the charge-coupling rows are linearly independent, restrict to the active charge subspace and set E=kappa^(-1). The eigenvalues of E^(1/2) D E^(1/2) are the squared velocities. The number m of modes then obeys m<=r<=d-1, so

    A_U(1)=sum_a v_a <= 6 exp(1) M(d-1).

This proves the requested relaxed bound for Gaussian infrareds independently resolved by microscopic charges. The full proof, including the gauge optimization and the matching assumptions, is in PROOF.md Sections 6 through 8.

[Proved, limitation] Having nonzero coupling to a charge does not imply row independence. With one charge, m equal-speed modes and Z_a=1/sqrt(m), both D=v/pi and kappa=1/(pi v) are independent of m, whereas A=m v. This is a continuum obstruction to an inference from the response matrices, not a fixed-d lattice counterexample. For general rectangular Z the ratio D/kappa measures weighted squared velocities; it does not count modes. The general bound on the sum over arbitrarily many modes sharing the microscopic charges remains open. Independent neutral sectors require separate control.

[Proved, convention qualification] If left and right speeds differ, thermodynamic relaxation of the current gives D e''=v_R v_L in the drifting scalar Gaussian model, rather than the square of their arithmetic mean. Locality bounds both chiral speeds, so the velocity and conditional aggregate bounds survive. The simple matrix eigenvalue identity above assumes no drift.

[Numerical] Small deterministic checks passed for 30 local replacement-channel cases, the matrix response identity, the one-charge multiplicity example, 99 free-chain densities, and 101 XXZ anisotropies. Residuals are below 1e-12. See CHECKS.md and check.py. These checks do not test infrared matching. No heavy numerics, cluster jobs, or writes outside work_C108 were performed.
