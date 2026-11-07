# C126 result

2026-11-06. The general bound remains open. Full derivations are in [PROOF.md](PROOF.md), and the numerical comparison is in [TABLES.md](TABLES.md).

[Proved] The permutation seam is an exact reformulation. Replacing the closing bonds by sum_a h_((L-1,a),(0,a+1)) gives one NL ring, so its ground-energy shift is E_(NL)-N E_L. This is a bond modification, not conjugation of the entire replicated Hamiltonian. Dividing by N and taking N to infinity at fixed L gives e0 L-E_L. Within the assumed CFT vacuum sector, the twist weights c(N-1/N)/24 give the requested coefficient pi A(N-1/N)/(6L).

[Proved, for the stated quadratic CAR and parity-resolved spin realizations] Spreading the continuous copy rotation gives

    E_(NL)-N E_L <= [pi^2 M/(6L)](N-1/N),
    A <= pi M.

Number-conserving copies have U(N), while general paired copies have O(N). Odd cyclic permutations lie in SO(N) and suffice. Ordinary spin swaps and graded fermion swaps differ. The spin construction uses a signed even cycle so that its unwrapped boundary sign is -p^N, the sign required by the trial parity p^N. This proves admissibility and gives the finite-size bound D(L)<=pi^2 M/(4L) in that quadratic spin class. Omitting these distinctions gives a wrong boundary sector. The spread bound is weaker than the inherited sharp all-quasi-free A<=2M. Its corresponding doubled-ring coefficient pi^2/4 is also different in scope from C119's all-temperature coefficient 10pi/3. The trial is optimal for the half-filled XX chain when its actual hopping expectation is retained.

[Proved] Exact on-site ground-preserving paths fail already on the Heisenberg vacuum, for every number of copies. Its two-rung correlation matrix singles out the direct sum of the copy su(2) algebras. A continuous stabilizer path cannot permute these simple ideals. Thus the identity and a nontrivial copy permutation lie in different components of that stabilizer. This excludes paths preserving the product ground vector exactly. It does not exclude spatially extended or approximately ground-preserving constructions.

[Proved] The entire N=2 isotypic-phase trial family is obstructed. For arbitrary local fractional swaps exp(i theta_x S_x), the Heisenberg trial cost is

    F=sum_x [a(1-z_x z_(x+1))+b w_x w_(x+1)],
    z_x=cos(2theta_x),    w_x=sin(2theta_x),
    (z_L,w_L)=-(z_0,w_0),    a=-E_L/L,    b=4a^2/3.

A uniform spread has F=L(a+b)/2+O(1/L). Every angle profile on every even L>=4 obeys F>=1/3, whereas the true defect is pi^2/(8L)+o(1/L). Optimizing the path therefore cannot give the required 1/L cost. This conclusion does not depend on numerical optimization.

[Proved] Large N does not automatically remove the bulk potential. For N odd prime, the principal cyclic logarithm applied halfway to N copies of the ferromagnetic Ising cat vacuum raises the energy per copy per ordinary bond by

    (2J/3)(1-1/N^2) -> 2J/3.

Direct matrix checks at N=3,5,7 agree to 1.67e-16. The exact permutation defect in this example is zero, since all copies can choose the same polarized phase. This is a counterexample to automatic low-energy preservation of the specified product vacuum, not to the Casimir bound or to every possible large-N trial.

[Proved] Finite-group LSM supplies no missing estimate from replica permutations alone. Their on-site representation is linear and has invariant vectors. An ordinary Ising spin-flip twist costs 2J, but an unrestricted identical-copy permutation twist has zero cost in that model. Fixing opposite copy phases instead forces two walls, costing 4J, even on compatible lengths. The distinction is between the unrestricted ground energy and a fixed sector or an inadequately correlated trial state. For quantum ordered phases, excluding a nonzero asymptotic wall cost requires a finite-size ground-energy input, stated explicitly in PROOF.md rather than inferred from the word "compatible."

[Numerical] ED reaches 2L=20 for Heisenberg and XXZ at Delta=0.5,2. It reaches 2L=12 for SU(3). The coherent seam family from C114 is the best of the evaluated trials at every reported point. Its ratio to the exact defect still increases. Representative ratios are below. "Fractional" is the best multistart profile found, not a certified optimum.

| Model | L range | Fractional / exact | Coherent / exact |
| :--- | :--- | :--- | :--- |
| Heisenberg | 4 to 10 | 2.2876 to 5.7003 | 1.1095 to 1.3728 |
| XXZ, Delta=0.5 | 4 to 10 | 2.4198 to 6.0957 | 1.1460 to 1.4167 |
| XXZ, Delta=2 | 4 to 10 | 2.7190 to 8.7458 | 1.2506 to 1.9055 |
| SU(3) swap | 3 to 6 | 1.4343 to 3.4203 | 1.1474 to 1.4959 |

[Numerical] Maximum eigenpair residual is 1.06e-11. Full-space small-ring comparisons agree to 1.34e-15. The fractional-swap local contraction agrees to 8.89e-16, and an independent eight-spin vector expectation agrees to 1.78e-15. The increasing ratios are evidence against the trial families, not an asymptotic fit. In particular, a ratio to an exponentially small gapped defect can grow even for a useful 1/L upper bound. The analytic Heisenberg obstruction avoids that ambiguity.

[Assumed sufficient lemma, unresolved] PROOF.md section 9 gives one precise missing estimate. For H(s)=H_L^(2)+s(B_sigma-B), the integrated ground-energy curvature must cancel the initial seam slope to accuracy C(d)M/L. In a nondegenerate interpolation it asks for

    2 integral_0^1 (1-s) sum_(n>0)
      |<n_s|B_sigma-B|0_s>|^2/(E_n(s)-E_0(s)) ds
        >= <B_sigma-B>_(s=0)-C(d)M/L.

The curvature-measure formulation covers crossings. This would imply the target, but is equivalent to it after an exact identity. No weaker sufficient lemma was proved. The disconnected stabilizer is a new obstruction to the suggested on-site route. The bulk potential also relates to lesson 2, and the remaining integrated reconnection estimate encounters lesson 4. No counterexample to a finite C(d) was found.

[Work record] [checks.py](checks.py) and [algebra_checks.py](algebra_checks.py) reproduce the completed checks. The 18-site SU(3) sector has 17,153,136 states and was not run locally. [run_typhon.slurm](run_typhon.slurm) prepares that ED comparison with 128 GB inside a SLURM allocation. The sandbox cannot reach Typhon, so no job was submitted. All writes remain in this folder, and [STATE.md](STATE.md) is the local handoff instead of a project-level edit.


[Proved, strengthened ordered obstruction] The consolidated supplement [COPY_ROTATION_PROOF.md](COPY_ROTATION_PROOF.md), section 3, proves more than failure of the cyclic logarithm. For every N>=2 and L>=3, the minimum twisted excess over all products of arbitrary on-site rung unitaries applied to the Ising cat product is JN. The actual defect is zero. Thus the on-site-unitary obstruction persists per copy at large N. This does not exclude choosing an aligned ordered vacuum or coherently correlating the copies.

[Numerical, independent validation] [copy_rotation_checks.py](copy_rotation_checks.py) verifies explicit signed Gaussian trial vectors in both parity sectors, with pairing and with two modes per site. Overlapping independent interacting energies and trial costs agree within 3.78e-14. The current checks.py was rerun through 20 spin sites and 12 SU(3) sites. [COPY_ROTATION_RESULT.md](COPY_ROTATION_RESULT.md) records the supplementary statements and their limitations. The 18-site SU(3) run remains unsubmitted.


[Numerical, prepared larger trial] su3_pending_trial.py computes the L=9 SU(3) ground state and all trial costs without the large E18 calculation. It gives E9=-6.5794357590864685, optimized fractional excess 1.9518332831194254, and coherent excess 0.6769252653558792. The ground-sector dimension is 1680 and residual is 3.06e-12. These are recorded in su3_L9_trial.json. No trial-to-exact ratio is asserted until E18 is computed on the cluster.
