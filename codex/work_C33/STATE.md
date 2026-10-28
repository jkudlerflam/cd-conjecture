# C33 local handoff

[Proved] The folded 23-letter frozen chain has positive Hermitian second-order hopping and common-field crossings for every n under the explicit finite D(n), g(n), epsilon(n) schedule in PROOF.md Sections 3-6. This controls both one-particle hopping and the quadratic projection. It does not prove all-n Sidon margins.

[Proved, computer-assisted] At n=23 a loop phase choice within the same alphabet gives 46 simple roots and delta2>0.0140193, delta4>3.1558e-7, including repeated-point pair sums. Parameters are a=4,D=1000,beta=1/det(C), td=beta*2^-136, g=2^-365, normalized field 2^-214. Physical loop phases are (pi,pi/2,pi/2), neutral express phase pi/46, port phase pi/2. PROOF.md Section 9 and root_certificate_n23.json explain the finite-g certificate.

[Proved, computer-assisted] Tree positivity certificates reach n=32. Numerical leading matrices reach n=64. Full-in-g one-particle inverses reach n=12 at D=10 and n=8 at D=1000. Failed precision and timeout attempts are documented in PROGRESS.md and RESULT.md.

[Open] All-n distinct-center and Sidon theorem, interacting infrared conclusion, frozen length-sector audit and mobility. This work does not claim c=23 for the microscopic chain. C32 handles a different mobile positive-core design. No alphabet enlargement or mobility modification was made here.

[Scope] Project STATE.md lies outside the allowed write scope and was not changed. This file is the local handoff. Start with RESULT.md, then PROOF.md. Reproduce the certificate using the five commands in RESULT.md. All new code and artifacts are in this folder.
