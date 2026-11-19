# Neither (P) nor (D)

[proved] The missing ferromagnetic strip of the uniform U(1) qubit chain admits a length-uniform doubling bound. Put \(J=\sqrt{J_\perp^2+D^2}>0\) and \(J_z/J=-\cos\eta\). For \(0<\eta\le\pi/4\), every uniform field, every hopping phase and every \(L\ge2\),

\[
 \frac{L(E_{2L}-2E_L)}{M}\le 2\pi\sin\eta+\frac{\pi^2}{4}.
\]

[proved] At zero hopping phase relative to negative real exchange, let \(F_L(n)\) denote the zero-field ground energy at \(n\) down spins. The stronger estimate \(L|F_{2L}(2n)-2F_L(n)|\le\pi J\sin\eta\) holds in every sector. An interpolation that splits each Bethe quantum number into two makes the energy derivative an alternating sum of dressed velocities. Their total variation is at most \(4\pi J\sin\eta\) throughout the interpolation.

[proved] In the strip the Casimir lower bound holds with \(K=(2/3)(2\pi\sin\eta+\pi^2/4)\), by dyadic summation. An exact vacuum convention at prescribed down-spin density \(\rho\) requires \(L\rho\) integral and trivial accumulated hopping phase. The estimates hold for every integer \(L\ge2\).

[assumed] Using the complementary bounds recorded in sections 9.91, 9.94, 9.104 and 9.106 of the supplied report, this closes the entire uniform U(1) qubit class. An explicit common constant is the maximum of \(2\pi+\pi^2/4\), \(6\cdot10^{15}(1-1/\sqrt 2)^{-25/3}\), \(4078\) and \(10781\). The new strip estimate is independent of those inputs.

[proved] A separate bound rules out an unbounded velocity sum in the effective cell Hamiltonian (9) of PROOF.md. Its hypotheses are isolated gapped cells, one microscopic interface bond of strength at most the cell gap, and a sum of the specified same-interface first- and second-order terms that are individually number-conserving quadratic operators on the encoded modes. If \(J_*\) bounds microscopic bond norms in this decomposition, the sum of absolute Fermi velocities in physical lattice units is at most \(C_*d^4 J_*\), with an explicit universal \(C_*\) in the proof. The bound is uniform in the gap, the cell length and the number of encoded modes. A harmonic tower is therefore excluded for this mechanism when \(J_*/M\) is bounded.

[proved] In that cell model the number of chiral crossings of speed at least \(u\) is at most a universal constant times \(d^4[1+\log_+^2(J_*/u)]\), with a universal cutoff in \(u/J_*\). Approximating the virtual energy denominator bounds the number of crossings at a given speed. A locality estimate for adjacent projected interfaces bounds the largest speed when their separation grows.

[numerical] The Bethe endpoint energies agreed with the spin Hamiltonian to within \(7\cdot10^{-13}J\) in 60 comparisons. The spin calculation passed 68 reference checks against XX and the isotropic ferromagnet. The script also tests the doubled energy, the sampled velocity variation and the interpolation derivative. Run with **python3 -B check_xxz_doubling.py**.

[proved] The spin proof requires the ordered counting equation (2). Even within XXZ, its sufficient kernel inequality (6) fails for \(\pi/4<\eta<\pi/2\). The cell proof requires (9), in particular the absence of uncontrolled terms coupling the two ends of a cell. A general interacting bond supplies neither representation. Higher perturbative orders in a proposed encoder would also need uniform control of their contribution to the microscopic Casimir energy. These are the stopping points for the two routes.
