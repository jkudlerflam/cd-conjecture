# C27 progress

2026-10-28. G7 only. No preexisting progress file.

[proved, source accounting] Read project STATE.md, route-A background and assembly, STATUS.md, C20 Sections 1-4, C13 Section 6, C8 joint_fluid.md and length_fluid.md, and C19 grading conventions. Read C21/C22 progress and C8 run-8 progress. C25/C26 and the three named Claude contribution progress files were absent at the first read. No compiler or second-order wrapper work is being duplicated.

[assumed] Work within C8's projected wrapper with identity payload transport, an ordinary tensor product of cell codes, and one common logical parity. The microscopic compiler and higher-order preservation of the projected description are separate gaps.

[proved, derivation in preparation] The common Jordan-Wigner boundary sign is (-1)^(N_total-1), shared by every flavour. The ordinary cell translation, not a graded permutation, determines the payload character. Its many-body momentum is integral even though the individual orbital grid can be half-integral. C8's quotient then gives m_gap+m_payload+ell=0 modulo M.

[heuristic, calculation being audited] Physical-coordinate harmonics have an integer reciprocal branch D and momentum rho(q_payload+2 pi D) modulo 2 pi. Their length-field dressing must change with that branch. At a physical resonance, q_payload+2 pi D=2 pi t(b+nu), so the cell resonance q_payload-2 pi t nu=0 modulo 2 pi is recovered. The free length contribution to the dimension becomes t^2 there. This could resolve C20's gluing caveat without identifying all nonresonant physical dimensions with cell dimensions. A derivation and finite-ring checks are next.

[proved, workspace restriction] All writes stay in this folder. Project STATE.md will not be changed because the user explicitly forbids writes outside it. No network, external computation, or subagents are used.

## Step 1 completed, finite gluing and common boundary sign

[proved] The ordinary cyclic cell shift equals second-quantized twisted fermion translation in the total-parity sector with seam sign tau=-sigma. The reordering and seam signs multiply to (-1)^[N_0(N_total-N_0)+N_0(N_total-1)]=1. Every flavour therefore uses beta=(N_total-1)/2 modulo one, with no separate flavour spin structures. The relative winding congruence is 2 delta c_alpha = delta N_total+delta N_alpha modulo two. For parity-even perturbations it reduces to the usual fermion congruence, together with even total charge change.

[numerical, completed] gluing_checks.py checks every energy within each ordinary cell-translation sector for common-graded free payloads (M,n)=(3,1),(3,2),(4,2),(5,2), dimensions 8 through 1024. Largest error 2.31e-14. Independent physical marker graphs at (b,M,r,n)=(2,3,1,2),(2,3,2,2),(2,4,1,1), dimensions 448,512,144, agree with the full quotient spectrum for both number-conserving and parity-only pairing/conversion payloads. Largest full-spectrum error 2.67e-14. These are finite-matrix checks, not central-charge estimates.

## Step 2 completed, continuum character and resonance comparison

[heuristic, standard bosonization input explicitly required] With x(y)=(b+nu)y-phi_l(y)/pi, an anchored cell operator with cell harmonic m and reciprocal branch s has physical momentum rho(q_payload-2 pi m nu+2 pi s). Its length-field coefficient is 2m+rho(q_payload-2 pi m nu+2 pi s)/pi. Set D=s+b m. The physical vertex data are Q_T=rho(q_payload+2 pi D) modulo 2 pi and length coefficient 2rho(D+q_payload/(2 pi)). This is one integer family, not two independent winding labels.

[proved, algebra conditional on that continuum map] Physical resonance means Q_T=2 pi t for an integer t. It is equivalent to q_payload-2 pi t nu being in 2 pi Z, and the length-field coefficient then becomes 2t. Thus C20's parity-even local lattice and resonance/dimension test apply on resonances. Nonresonant physical dimensions are generally different. A pure marker harmonic has free dimension rho^2 D^2, a useful consistency check.

[assumed] Stability under integration of the actual microscopic gapped modes and the existence of the spatial compiler remain separate. G7 is resolved for the specified projected identity-transport wrapper, with the continuum identification stated as physics input. Final derivation and audit are in preparation.
