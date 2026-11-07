# C127 local handoff

2026-11-06. Only work_C127 was modified. The parent project STATE.md is outside the authorized write scope, so this file carries the substantive handoff.

[Proved] The fixed-frequency insertion map has rank <=d^2-1, but its spectral cyclic span need not. Actual XX density insertions at fixed Matsubara k have infinitely many positive reflected spectral masses. PROOF section 3 gives their energies and weights, with the necessary reflection sign at odd separation. At k=1 the count below fixed E grows as beta sinh(E/2)/(2pi). Do not interpret insertion rank as an eigenstate or species count.

[Proved] The C125 low-frequency cone applies to the detected cyclic subspace if its reflected positive spectral representation exists, with V=4eM and domain 0<|omega|<V. This does not reconstruct all crossed states. At high frequency use the logarithmic LR rate. Exact XX disproves a global momentum-preserving linear Fock majorant even with exp(K); that failure does not refute a scalar trace bound.

[Proved] For XX, all beta and even L>=2: log Ztilde<=log 3+beta[csc(pi/L)-L/pi]<=log 3+(2-4/pi)beta/L. Convex midpoint quadrature of the exact crossed product proves it. Consequently a scalar bosonic product holds with S=2,V'=4eM,K=log 3-2log(1-exp(-1)). Tensor products of independent XX chains inherit explicit bounds. See PROOF section 5.

[Proved] The diagonal ferromagnetic d-state chain forces K(d)>=log d in every universal oscillator-only bound. If a separately assumed Gaussian charge-sector comparison holds with sum rho<=4M, its log trace is <=16M beta/L with no species-count factor. The stiffness budget alone does not establish that spectral comparison.

[Proved] Existence of a scalar bosonic product bound with arbitrary dimension-dependent S,K and V'=cM is equivalent, up to explicit changes in constants, to the scalar trace target. The elementary floor-count lower bound on the oscillator product proves the converse. A logarithmic ultraviolet oscillator dispersion also yields an O(S M beta/L) reference bound. The missing input remains an energy-filtered comparison of the complete many-insertion space to such a reference.

[Numerical] checks.py and audit_checks.py pass their stated assertions. Exact XX has 30 samples. Momentum-resolved transfer tests have 48 cases, m=2,3,4, beta=0.25,1,2,8, Delta=0,0.5,1,2. All 240 scalar samples pass S=6,V=4eM,K=log 3. Sectorwise comparison fails. A continuum-cone check also fails at coarse m=2,beta=8 for Heisenberg, so finite-step data must not be promoted to a continuum cone theorem. The L=2 exact XX crossed product has a controlled 0.00519 tail at beta=64, not machine-precision agreement. CHECKS.md gives all errors and qualifications.

[Assumed or unresolved] No general interacting scalar majorant and no lattice counterexample to it were established. Tasks 2 and 4 remain incomplete for that scalar question. Do not report the rank obstruction as a solution of C(d). No independent agent/referee was used. run_typhon.slurm is prepared for m=5,6, unrun and unsubmitted.

[Work record] RESULT.md is the summary, PROOF.md contains the derivations, CHECKS.md gives the scope, and PROGRESS.md records the session. checks.json and audit.json are generated numerical outputs, not orientation files.
