# C108 progress

[Proved, scope] Read the project memory, the full background note with both updates, REPORT Sections 9.68 through 9.77, C104's local memory and derivations, and the pairing README. All writes stay in work_C108. No network shell commands or heavy numerical calculations have been run.

[Proved, inherited input] The charge stiffness satisfies D_q <= sigma_q^2 M. C104 corrects the Gaussian susceptibility to kappa = Z^T V^(-1) Z/pi and explains why nonzero charge coupling does not bound the number of modes by the charge rank.

[Heuristic, work in progress] A local channel that replaces nonmaximal charge states by maximal charge states should give a finite saturation field, hence a bound on the total distributional curvature of the convex canonical energy. The charge gap at the spectral endpoints, rather than the charge spread alone, fixes the normalization.

[Heuristic, work in progress] The Lieb-Robinson cone should bound every Gaussian velocity with a nonzero local density response, and thus bound D e'' without an independent bound on e''. A proof must state the real-time infrared matching assumption and distinguish full-rank charge response from many modes sharing one charge.

[Proved, completed lattice derivation] The local replacement channel gives e(n)>=e(q_+)-4M(q_+-n)/delta_+ and the corresponding lower-end inequality. Convexity then bounds the interior chemical potentials and the entire distributional curvature measure. Degeneracy of endpoint eigenspaces causes no difficulty.

[Proved, conditional pointwise result] An explicit nearest-neighbor commutator path count gives propagation speed 6 exp(1) M. Assuming the local density has its Gaussian real-time scaling limit, every nonzero charged-mode row obeys v_a<=6 exp(1) M. Hence D e''<=36 exp(2) M^2 in the non-drifting one-mode theory. The proof uses the scaling limit at each fixed Hamiltonian and requires no uniform matching scale.

[Proved, completed multi-mode analysis] On the full-row-rank Gaussian subclass the eigenvalues of kappa^(-1/2) D kappa^(-1/2) are all the physical v_a^2. This yields A_U(1)<=6 exp(1) M(d-1). Mere nonzero coupling does not bound the number of modes, and the corresponding general aggregate bound remains open.

[Proved, corrections and scope] The Mott doping edge has e'' tending to zero when v tends to zero and K stays finite. A plateau cusp instead contributes an atom to distributional curvature. The scalar identity D e''=v^2 also needs the non-drifting convention; with unequal chiral speeds and relaxed current it gives v_R v_L. These qualifications are explicit in PROOF.md.

[Numerical] check.py passed 30 replacement-channel cases, a positive-matrix generalized eigenvalue check, the one-charge multiplicity example, 99 XX densities, and 101 XXZ anisotropies. Largest channel inequality residual was 9.550e-16. Largest matrix-spectrum residual was 8.882e-16. No heavy calculations were used.

[Proved, completion] RESULT.md, PROOF.md, CHECKS.md, check.py, and local STATE.md record the conclusions and the remaining gaps. Literature was checked only with the browsing tool. No parent report, parent project memory, or other work folder was edited, in accordance with the task's write restriction.
