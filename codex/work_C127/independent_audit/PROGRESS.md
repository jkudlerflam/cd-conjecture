# C127 independent audit progress

2026-11-06.

[Assumed] This run read the target, background, reconnection note, REPORT 9.76 through 9.88, C125 files and project memory. All writes remain inside work_C127.

[Proved] The insertion source rank does not bound its spectral cyclic dimension. Exact XX gives a lattice witness. Its full excitation trace obeys log Ztilde <= log 3 + 8 beta/(3 pi L), and a scalar boson product bound with S=6, V=4eM and K=log 3+8/(3e).

[Proved] The unrestricted scalar product target is equivalent, with changes of constants, to K+kappa M beta/L. A logarithmic high-frequency reference and Gaussian charge factors can each be bounded at this scale once the required spectral comparisons are assumed.

[Numerical] Initial checks covered 48 transfer matrices, including m=1, and 480 folded momentum-sector comparisons. Kstar passed the scalar samples; K=log 3 failed at two m=1 points. One Sz frequency source detected up to ten distinct energies.

[Assumed or unresolved] During verification, checks.py and checks.json in the parent folder were replaced concurrently by an unrelated implementation, and PROGRESS.md gained entries from that run. Those edits are preserved. This subfolder isolates this run and its output schemas. The initial audit rerun failed with KeyError: transfer because the parent JSON schema changed. All isolated checks will be rerun before delivery.

[Assumed or unresolved] No universal interacting scalar product theorem or lattice counterexample to it is established.

[Numerical] Both isolated scripts completed successfully. The 48 transfer matrices and 480 folded-sector comparisons reproduce the original diagnostics. Additional L=2 exact XX checks retain explicit product-tail errors. A diagonal ferromagnetic lattice test verifies the necessary K>=log d mechanism.

[Assumed work record] RESULT.md, PROOF.md, CHECKS.md and STATE.md now give the final limitations and file map. The universal interacting scalar question remains unresolved. A guarded SLURM job is prepared but not submitted.
