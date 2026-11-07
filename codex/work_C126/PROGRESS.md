# C126 progress

2026-11-06.

[Completed] Read the target, all background updates, reconnection note, requested REPORT sections, C125 result, and the relevant free and stiffness proofs. Read the project memory before exploratory work. All writes stayed in this folder.

[Proved] Verified the seam relabeling, its gauge invariance, the fixed-L large-N identity, and the CFT twist weights. Derived A<=pi M from a boundary-matched continuous free-copy twist, with explicit grading and O(N) qualifications.

[Proved] Established disconnected on-site stabilizers for the Heisenberg product vacuum at every replica number. Derived the fractional-swap functional and its L-independent 1/3 lower bound at all even L>=4. Derived the odd-prime large-N Ising-cat barrier and distinguished constrained phase walls from the unrestricted permutation defect.

[Numerical] Completed spin ED through 20 sites and SU(3) through 12. Compared bare gluing, a uniform fractional path, optimized fractional profiles, and the inherited optimal coherent seam family. Independently checked seam conjugation with complex oriented bonds, full-space sector choices on small rings, an eight-spin trial vector, the free phase conventions, and Ising cyclic rotations. Checks and numerical limits are in RESULT.md and TABLES.md.

[Unresolved] General C(d), approximate ground-preserving paths, and more general correlated large-N constructions remain open. The single integrated seam-screening lemma in PROOF.md is an equivalent quantitative remainder, not a new proved reduction.

[Pending external computation] The 18-site SU(3) ED comparison is SLURM-ready and unsubmitted. No connection to Typhon was attempted. Local STATE.md carries the handoff because outside writes are forbidden.


[Consolidated] Concurrent work in this directory was preserved and checked. COPY_ROTATION_PROOF.md supplies the signed-cycle admissibility proof, the arbitrary-rung-unitary Ising minimum JN for every N>=2 and L>=3, and the stronger Heisenberg profile lower bound 1/3. The shared PROOF.md, RESULT.md, and STATE.md now point to the supplement. copy_rotation_checks.py validates explicit Gaussian trial vectors, including odd initial parity and paired two-mode sites. Independent overlapping ED and trial costs agree to 3.78e-14.


[Numerical, prepared larger trial] su3_pending_trial.py computes the L=9 SU(3) ground state and all trial costs without the large E18 calculation. It gives E9=-6.5794357590864685, optimized fractional excess 1.9518332831194254, and coherent excess 0.6769252653558792. The ground-sector dimension is 1680 and residual is 3.06e-12. These are recorded in su3_L9_trial.json. No trial-to-exact ratio is asserted until E18 is computed on the cluster.
