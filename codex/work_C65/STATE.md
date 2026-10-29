# C65 handoff

2026-10-28. Explicit regulator obstruction and smooth trace-compatible repair. Lemma TP remains unproved.

[proved] PROOF Sections 2 and 3 show that C62's sharp fermion windows fail the common absolute moment norm before regulator removal. At fixed physical cutoff Lambda and T>0, the static chiral density bubble has coefficient -n_Lambda(1-n_Lambda)/(2 pi T) multiplying |p|. A smooth infrared subtraction cannot cancel it. It generates an actual connected g^2 same-chirality quartic, and antisymmetrization leaves a cusp. Finite first spatial moment would imply C1 momentum regularity, so the fifth-moment norm and its thermal replacement fail. This is not a refutation of C62's trace-limit order or of TP after regulator removal.

[proved] PROOF Section 5 repairs the regulator while keeping the same trace. On the excitation Fock space use j f=r f direct_sum sqrt(1-r^2)f and Phi(O)=Gamma(j)^*(O tensor I)Gamma(j). Normal-ordered monomials acquire one r per leg. For the cross-chirality g2 interaction this stays quartic. Positivity and Phi(H0)<=H0 preserve |V| in form sense <=b H0+C_L. Compactness and min-max give eigenvalue convergence; the uniform free Gibbs majorant gives trace convergence at fixed L,beta. The full scalar is C62's oscillator trace after removal.

[proved] Smooth cutoffs remove the static cusp, with uniform spatial moment control of the edge correction. Smooth infrared shells have uniform Gram and fixed-moment bounds, compatible with the collaborator's one-step p=3 localization/winding results. No full ultraviolet or all-depth connected expansion is proved.

[assumed next work] Prove PROOF Section 7 for this smooth trace-defined reference: cutoff-independent analytic ultraviolet preparation, common critical domains, full-kernel reference difference and mixed thermal difference, and scalar derivative/Hessian bounds. C62's double-scale sum and the collaborator's conditional fixed-mark allocation are already available. BM's smooth reference avoids the sharp boundary issue but still needs its own full scalar identification if used instead.

[numerical] checks.py passes, including independent finite CAR attenuation on system dimension 16 and dilation dimension 256, residual 2.23e-16. No RG simulation, subagents, shell network commands or writes outside this folder. Project STATE was left unchanged under the explicit write restriction. New arguments have not received independent or human review.
