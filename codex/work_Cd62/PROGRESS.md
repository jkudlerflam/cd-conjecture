# Progress

The four requested background files have been read. The target is a uniform upper bound on the positive doubling defect. The infrared velocity bound alone leaves open how many channels can contribute at a given length.

The direct length comparison in the attractive XXZ chain started from a duplication of the rapidities of the shorter ring. Moving each pair of Bethe quantum numbers apart makes the energy derivative pair alternate signs with samples of the discrete dressed velocity.

The kernel estimate holds for all \(0<\eta\le\pi/4\) with \(J_z/J=-\cos\eta\). On each side of a source, the variation from an observation point out to infinity is at most the kernel's value at that point. Differentiating the finite Gaudin equation cancels the kernel contributions after integration. The resulting matched-sector bound is \(L|F_{2L}(2n)-2F_L(n)|\le\pi J\sin\eta\). A convexity argument gives the interpolating rapidities even when pairs coincide. Positivity of the spin ground state and continuation from XX identify the endpoint energies.

The strip proof is in PROOF.md. The independent spin matrices passed both free references and the Bethe endpoint comparisons. The supplied results on the complementary anisotropy ranges give finiteness for the entire uniform on-site U(1) qubit class. The proof uses a scalar Bethe scattering equation with an ordered solution.

For the separate cell construction I obtained an absolute approximation of the second-order resolvent by a sum of products. Its rank grows quadratically with the logarithm of the requested accuracy and is independent of the largest virtual energy. A signed-velocity counting argument turns this into a count of fast Fermi crossings.

For well separated interfaces, spectral filtering makes adjacent projected interface terms nearly commute. In a number-conserving quadratic chain this bounds the largest group velocity, as can be seen by differentiating a two-cell Bloch propagator. Combining the two estimates gives a length-independent sum of velocities for a precise class of first- and second-order effective Hamiltonians. The outstanding issue for the general problem is control beyond this effective class, in particular extra terms mixing the encoded modes inside the cells and higher perturbative orders.

Both derivations and their stopping points are recorded in PROOF.md. The final run of check_xxz_doubling.py passed the reference, sector-energy, doubled-energy, velocity-variation and energy-derivative checks. RESULT.md distinguishes the new proofs from the bounds used from the supplied report. The outcome for arbitrary bonds is neither (P) nor (D).
