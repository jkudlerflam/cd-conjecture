# C127 independent handoff

2026-11-06.

[Assumed scope] This folder is a self-contained audit of the C127 crossed-majorization route. It was isolated after concurrent edits changed the parent checks.py and checks.json. Do not combine the schemas or numerical grids of the two runs. Parent edits were preserved. The project-level STATE.md is outside the authorized write scope and was not modified.

[Proved] The source rank per signed Matsubara frequency is <=d^2-1, but the spectral cyclic span is not so bounded. Exact XX density insertions detect infinitely many energies at fixed frequency. The low-frequency support result needs a compatible reflected positive spectral measure, beyond spectral compactness.

[Proved] Exact XX satisfies log Ztilde <= log 3+8 beta/(3 pi L) for even L>=2 and all beta. A scalar boson product follows with S=6,V'=4eM,K=log 3+8/(3e). Momentum-preserving global linear-energy comparison fails on the same lattice model. The diagonal ferromagnetic d-state chain forces K(d)>=log d.

[Proved] Existence of the scalar bosonic majorant with arbitrary K,S is equivalent up to constants to the original scalar trace target. The softened logarithmic ultraviolet reference and Gaussian charge factors admit beta/L bounds once their many-body spectral comparisons are assumed. The missing lattice input is an additive many-insertion energy/counting comparison, not merely source completeness or two-point decay.

[Numerical] checks.py and audit.py pass in isolation. They compute 48 momentum-resolved transfer matrices of dimension <=256, 192 scalar comparisons and 480 folded-sector comparisons. Kstar passes scalar samples, K=log 3 fails two m=1 cases, and 56 momentum-sector comparisons fail. The scripts keep Trotter errors, product tails, and diagonalization residuals separate. Larger m is guarded by SLURM and has not run.

[Assumed or unresolved] No proof or counterexample for the unrestricted interacting scalar target. No minimal sufficient lattice hypothesis has been derived. RESULT.md gives the verdict; PROOF.md gives full arguments; CHECKS.md records reproduction, failures and limits. run_typhon.slurm is prepared but unsubmitted. No outside write, network command, or separate-agent review occurred.
