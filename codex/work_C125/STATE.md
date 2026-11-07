# C125 local handoff

2026-11-06.

[Proved] The main result is a negative implication theorem. Even a positive finite-dimensional crossed spectrum, an entire momentum tilt, the full global cone epsilon-epsilon_0>=|p|/V, the microscopic ratio ceiling, and the full infrared CFT torus form do not bound A at fixed labels d,M. Use the rank-truncated n-boson construction in PROOF section 5. It has A=n v, R<=beta M/2, and Theta''=-b/(4L^3) below a crossover proportional to n. It is not a local-chain counterexample.

[Proved] Actual XX invalidates the proposed global Matsubara cone and real tilt. Locally detected gaps are asinh(pi/beta)+asinh((2m-1)pi/beta) at p=2pi m/beta. Their logarithmic growth makes the real tilt divergent at every nonzero a. C109 instead implies a locally detected low-frequency cone and a logarithmic high-frequency LR rate. See PROOF sections 3 and 4.

[Proved] Positive even-ring transfer spectra are constructed for zero-field XXZ Delta>=0, antiferromagnetic bilinear Heisenberg at any spin, and all SU(2) symmetric integrable bonds. The latter have positive rank-l reflected coefficients (k+1)/(2l(l+1)); a Racah finite-difference proof is included. Positive realignment of Trotter gates gives a positive two-site transfer, and spectral compactness recovers exact even-ring traces at L>=4. An exact XX product includes L=2. The on-site RP criterion fails for Delta<0 and same-representation SU(N) swap/Sym^k with N>2.

[Numerical] Delta=-0.9 gives R_6>R_4 at beta=21.928721377 and increasing c_A on L=4 through 12. SU(3) gives D(4)=-1.406345439. Its vacuum lengths pass the Theta test. All supplied Sym^2 tables pass both signs after adding the proved E_N=0 point, which is absent from the JSON. N=8 gives c_A(8)/M=7.125295634 and A/M=6.597344573.

[Work record] RESULT.md is the summary, PROOF.md contains derivations and qualifications, CHECKS.md records scope and numerical failures. checks.py, audit_checks.py and spin_certificates.py reproduce the local work. checks.json, audit.json and spin_certificates.txt are outputs, not orientation sources. run_typhon.slurm is prepared but unsubmitted. No outside file was modified. The project-level STATE.md is outside the authorized write folder, so this file serves as the permitted handoff.

[Assumed or unresolved] No fixed-d nearest-neighbor counterexample to the proposed A bound. No classification of arbitrary positive representations restricted to N>2 vacuum subsequences. No general reconstruction of time-circle rotations and local insertions from the compactness-built limiting transfer spectrum. These unresolved points do not affect the negative sufficiency theorem.
