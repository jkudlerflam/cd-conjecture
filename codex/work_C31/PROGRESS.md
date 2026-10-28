# C31 progress

[proved, scope] Read project memory and the required route-A background, assembly Section 7, switch-cell README, C26 result, and infrared control note. Read C8/C21/C22/C23 progress to avoid duplicating their work. No previous C31 progress existed. All writes stay in this folder; project STATE.md cannot be updated under the explicit task restriction.

[proved, inherited limitation] Equal-port symmetry makes T Hermitian. Exact total U(1) forbids pairing but not pair transfer. Distinct sea centers alone do not exclude every number-conserving quartic resonance. At zero onsite matrix the Bloch symbol has h(k+pi)=-h(k), which creates additional half-filling resonance constraints.

[heuristic, work plan] Derive the leading conversion-squared carrier resolvent, validate it against the existing finite-cell extraction, and test asymmetry, phases and a small uniform logical chemical potential for n=2,3,4. Seek an all-n nonresonance proof or isolate the missing finite statement. Separate the slowest velocity and visibility factors in the coupling hierarchy. The alternative centered-sea RG calculation is needed only if displacement does not furnish a controlled route.

## Leading carrier reduction completed

[proved] The lambda^2 tau^2 coefficient can be obtained by two-carrier Sylvester inverses. At tau=0 each fixed data occupation and carrier-type pair is a Kronecker sum of two (3n+3)-dimensional one-carrier matrices. Conversions act locally between these blocks. This eliminates scalar subtraction and the large many-body basis.

[numerical, n=2] The reduced calculation agrees with a literal 10,950-state two-cell calculation to 6.0e-14 relative in T and 1.9e-11 relative in M. The large calculation took 3.7 seconds; the reduction took 0.06 seconds.

[proved] At equal H/D energies and w_D=-w_H, staggered particle-hole conjugation gives T real and M=-M^T at mu=0. Its Bloch symmetry is h(pi-k)=-h(k)^T. Thus displaced centers at zero field do not suffice: same-chirality root pairs have sum pi and allow quartic umklapp. A uniform chemical potential breaks this constraint without enlarging the alphabet.

[numerical, n=2,3,4] Leading-order scans find 2n transversal roots and positive number-conserving quartic margins when mu=0.3 times the smallest symmetric hopping eigenvalue. Candidate (arm asymmetry, neutral phase) values are (0.05,0.001), (0.001,0.00001), and (0.00001,0.0000001). The n=4 margin is about 2e-6. These numbers still need symmetry-based numerical cleanup and an independent root check. No all-n displacement proof follows from this scan.

[proved, conversion counting] At tau=0 all logical data are spectators. Every connected traceless effective term has at least two conversions. A normal-ordered quartic needs at least four. Thus the lambda^3 coefficient at order tau^2 is Gaussian, while its non-Gaussian terms start at tau^4. This is separate from C26's rho argument.

## Infrared and independent checks

[proved] The symmetric leading hopping is a positive Gram matrix up to its negative overall sign. The heat-kernel functions have distinct first nonzero Taylor powers, proving full rank at every n. A quantitative small-time Gram bound supplies a strictly positive, though very small, lower bound as a function of n. This proves that the eigenvalue hierarchy is not a finite-rank obstruction.

[numerical, n=2,4] The leading coefficient was independently extracted in literal occupation spaces of dimensions 10,950 and 188,418. A tree gauge reduces each carrier-pair inverse to a positive sparse matrix. This improves the reduced calculation substantially: it agrees with the full extraction at about 1e-15 relative. The initial spectral-resolvent implementation lost some accuracy in n=4 momentum shifts and is retained only as a diagnostic option.

[numerical, n=2] The quartic coefficients of lambda^2 tau^4 and lambda^3 tau^4 were checked against an independent 11,680-state extraction, using all 11 vacuum/one/two-particle code columns. Relative errors were about 1e-13 and 5e-14 with the initial spectral implementation.

[numerical, n=2,3,4; heuristic infrared interpretation] At symmetric arms, phase zero, and mu=0.3*s_min, the leading interflavour opposite-chirality forward couplings are negative, the same-flavour couplings positive, and all pair-transfer dimensions exceed two. Their minimum excess divided by tau^2 is 9.0884, 16.5500, 42.2787. The one-loop Cooper/forward flow remains in the irrelevant basin in all three sizes. Cubic bare quartic vertices have the same favorable signs in these cases. No all-n sign theorem has been obtained.

[proved, limitation] These checks do not prove distinct displaced centers or favorable relative stiffness for arbitrary n. The all-n conclusion must remain conditional. The finite-coupling hierarchy closes at each fixed n once one of those infrared conditions is supplied, but that algebra cannot establish the missing condition.

## Final handoff

[proved, file accounting] RESULT.md is complete. It contains the microscopic symmetry proof, the all-n Gram rank proof and explicit lower bound, the missing all-n conditions, displaced targets and complete marginal-momentum enumeration, actual quartic extraction, one-loop equations and flows, and conditional physical energy scales. Only this folder was modified. Project STATE.md was intentionally left untouched under the explicit task restriction.

[numerical, final validation] The positive tree-gauge inverse agrees with literal occupation-space leading coefficients at about 1e-15 relative. At the actual n=4 displaced target the largest independent root discrepancy is 1.68e-13 radians, much smaller than delta_4=2.16352e-6. The final n=2 quartic comparisons, dimension 11,680, agree to 9.53e-16 and 5.58e-16 relative for lambda^2 tau^4 and lambda^3 tau^4. All Python files parse. The final one-loop runs preserve positive pair-transfer dimensions at n=2,3,4.

[heuristic, verdict] Stable fixed-pattern c=n phases are supported at n=2,3,4 under the stated local-integration and bosonization inputs. The general-n conclusion remains open. Continue with an analytic sign decomposition for the centered forward interactions, or an all-n proof of separated displaced centers. No selector, mobile-wrapper construction, full finite-coupling thermodynamic computation or central-charge measurement was performed.
