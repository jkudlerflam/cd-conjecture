# C126 handoff

2026-11-06.

[Proved] Permutation seam equals the NL ring after relabeling. CFT and fixed-L large-N limits agree with the task. Distinguish modifying a seam bond sum from conjugating the whole Hamiltonian.

[Proved for the specified CAR and parity-resolved spin realizations] Continuous free-copy rotations give A<=pi M. Number conservation gives U(N), pairing gives O(N), and odd cycles suffice. Spin tensor permutations are not graded mode permutations. Signed even cycles give boundary sign -p^N and establish the finite-size D(L)<=pi^2 M/(4L) bound for the quadratic spin realization. The inherited sharp all-quasi-free A<=2M is stronger.

[Proved] On the Heisenberg singlet vacuum at even L>=4, for every N, an exact rung-unitary stabilizer path cannot connect the identity to a nontrivial copy permutation. Two-rung correlation eigenvalue t=4e/3 isolates su(2)^N, whose ideals cannot be continuously permuted.

[Proved] Every two-copy fractional-swap profile on this Heisenberg sequence costs at least 1/3. The exact functional is F=sum[a(1-zz')+b ww'] with an antiperiodic seam, a=-e in [3/8,1/2], b=4a^2/3. The stronger lower bound uses F>=2(a-|b|), from multilinearity on the antiperiodic cube. Uniform spread costs L(a+b)/2+O(1/L). The optimized profile cannot give 1/L.

[Proved] Odd-prime-N cyclic logarithms on replicated Ising cat vacua have halfway excess per copy per bond 2J(1-1/N^2)/3. Large copy number does not make this on-site path low energy. The unrestricted defect is nevertheless zero, with aligned polarized phases.

[Numerical] checks_local.json and TABLES.md cover doubled spin sizes through 20 and SU(3) through 12. Fractional and inherited coherent seam ratios increase. Max eigenpair residual 1.06e-11. algebra_checks.json records full-space sector checks, actual seam conjugation, full-vector fractional energy, cyclic-cat checks, and free sewing checks.

[Unresolved] No general finite C(d), no counterexample to it, no theorem excluding approximately ground-preserving or spatially extended paths. PROOF.md section 9 states the integrated seam-screening lemma, equivalent to the target by concavity and spectral perturbation theory.

[Operational] All files are local to work_C126. Project STATE.md and REPORT.md were not modified under the task's outside-write prohibition. run_typhon.slurm is ready but unsubmitted. The pending SU(3) L=9 versus 18 comparison needs a 17,153,136-dimensional balanced sector and was not run locally. No network commands or Typhon attempts were made. Literature used only the browsing tool. No independent review.


[Consolidated additional theorem] COPY_ROTATION_PROOF.md section 3 proves that all products of arbitrary on-site U(2^N) rung unitaries applied to N symmetric Ising cats have minimum twisted excess JN for every N>=2 and L>=3. The true defect is zero. The proof relaxes local squared unitary entries to doubly stochastic matrices, reduces to bijections, and applies the Hamming triangle inequality. This is stronger than failure of one cyclic-logarithm path, but can be evaded by vacuum selection or coherent filtering. The coherent seam trial has cost 4J/(L+1) on these cats.

[Numerical supplement] copy_rotation_checks.py and copy_rotation_validation.json validate seam spectra, one- and two-mode Gaussian signed trial vectors, both parity signs, the all-N Ising orbit obstruction, and the Heisenberg stabilizer. COPY_ROTATION_RESULT.md is the supplementary summary. Independent overlapping interacting calculations agree to 3.78e-14. The primary checks.py was rerun through 20 spin sites and 12 SU(3) sites. The larger SU(3) comparison remains pending.
