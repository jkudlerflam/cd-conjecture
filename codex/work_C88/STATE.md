# C88 handoff

[Completed] Read RESULT.md and PROOF.md. The d=13 bond has eleven C84 material letters with the pin removed and two C83-style vacancy colors. All work stays in work_C88. No project-level memory was edited under the user's folder restriction.

[Proved] The automaton has o=0 on L and markers, o=1 on R,D,X. Incoming sets are {0,1} on R_0,E, {1} on R_1, and {0} otherwise. V_q accepts and outputs q. Colorings are unique for nonempty occupancy and N_matbad<=N_bad. Pure-vacancy rings have two constant legal colors.

[Proved] The physical-ring operator F counts squeezed L_1-marker adjacencies and commutes with every literal term. Full finite-K material restriction to F=0 equals C85's. C85's material certificate therefore applies unchanged; finite-K mobile corrections differ from its three-color lift. Frozen pin-free ground-state selection is proved by extending F3. Mobile selection is the explicit assumed S_pin, not a consequence of conservation.

[Proved] D2=-(t^2 C+J^2 A+2J^2 lambda^2 Z)/(JK). C has coefficient one on R_0,E contacts with vacancies and zero elsewhere. It is a bit-dependent contact field. Mixed hopping/material term is zero. A and Z are inherited vacuum-tail and wrong-end port counters. F4 uses kappa=JK and fixed-order volume-independent local remainder. Mobile K conditions are PROOF (30), implicit in C3 and the infrared tolerance.

[Proved, regenerated] n=13 active bits, p=14 data sites, b=15, 29 carrier vertices, 841 product vertices, 26 simple roots, 351 pair sums. Delta2>0.0156742281, delta4>3.2513114e-6. D=1000, Delta=4-2cos(pi/16), phi_A=pi/28, phi_B=0, psi_A=pi, psi_B=pi/2, t_d=beta*2^-66, beta=1/d_15 with recurrence d_(j+1)=1004d_j-d_(j-1), g=2^-272, lambda=2^-608, K>=2^2560, field=lambda^2 alpha^4 beta^2 g^2 2^-120. The root routine's -228 is an intermediate threshold, not the final g. Selector midpoint gives mu_sel=1.8891880765550189 and density gap>2.8116274e-5. The ten certificate sources and exact root, budget and selector outputs match C85.

[Numerical] All eight scripts in run_all.py pass. Legal correction enumeration reaches 19267 states. Nonzero-bit/port/two-vacancy leakage is K^-2 in dimension 1045; additional two-cell/two-vacancy port closure has dimension 3078. All illegal states reachable from the seeds are included. F3 uses complete 576-state sectors; F4 uses 60 states. Pin conservation exhausts 28561 words and 50024 directed moves. A legal component of 3528 states supplies single-letter creation witnesses for every letter, ruling out further deletion by the same argument.

[Assumed; heuristic conclusion] c=14>12 at d=13 remains conditional on massive matching, finite-J vacancy matching, mobile dressed-defect selection including pins, multichannel stability and conformal identification. No central-charge measurement or unconditional mobile pin-selection theorem. A different alphabet redesign below d=13 remains open.
