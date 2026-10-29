# C117 progress

2026-10-29.

[Proved work record] Read the project memory, the full bound background including all updates, and REPORT sections 9.68 through 9.78. Read C100 spectral conventions and the C103 and C106 results. All writes stay in work_C117.

[Open] No general crossover-free bound on A has been obtained. The targets are a finite-support obstruction for local marginal SDPs, small explicit dual certificates, and a genuinely infrared filtered spectral formulation.

[Proved work record] The authorized Typhon connection was attempted with LogLevel=ERROR, BatchMode=yes, and ConnectTimeout=10. It exited 255 with "Connection timed out during banner exchange". No remote job ran. Large numerics will be left as a SLURM-ready driver. Local work is limited to small matrices and solver checks. No package installation or delegation.

[Proved conditional implication, being written] A fixed-support translation-invariant local marginal hierarchy has an energy density independent of the circumference once all supports fit. A shorter ring supplies a feasible point. This prevents fixed-support certificates from resolving the asymptotic Casimir correction. The statement does not cover global moment matrices or explicitly circumference-dependent constraints.

[Assumed infrared identification] A low-frequency filtered inverse energy-density moment measures A with coefficient 1/(24 pi), using the momentum-averaged conventions of C100. Its microscopic upper bound remains open.

## Completed analysis and checks

[Proved] PROOF.md establishes the fixed-support obstruction and its conditional sublinear-support extension. Exact three-site optima were obtained for Heisenberg, XXZ at Delta=1/2, and SU(3) swap. It also derives the constant-3 local q^4 coefficient bound and an optimized-gauge massive Ising example with a nonzero q^4 term despite A=0.

[Proved] Fifteen local dual certificates and ten small full-ring intervals passed independent integer checks. The reference energy intervals have width below 3.83e-6. Solver optimality is not claimed for the random examples. The random d=3 k=3 searches hit the iteration cap; their certified lower endpoints remain valid. CHECKS.md records all results and limitations.

[Assumed, not proved] The spectral route needs a bound on an infrared-filtered inverse moment, yielding constant 24pi K. The ring route needs paired upper/lower certificates controlling L(u_2L-2ell_L), yielding constant 4K/pi. The modular route still needs the C106 quantitative field-counting and support lemmas, with the inherited C109 speed constant 4e. No general C M f(d) bound has been obtained.

[Proved work record] RESULT.md, PROOF.md, CHECKS.md, and local STATE.md are written. The larger block search is prepared as run_typhon.slurm and guarded by a SLURM allocation check. It was not submitted because the one SSH attempt timed out. All local computations were small validation calculations. Parent files remain unchanged.
