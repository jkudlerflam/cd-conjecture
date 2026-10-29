# C77 state

Completed 2026-10-28. Task is the strongest supported restriction on the
smallest local dimension of a counterexample to c<=d-1.

Read RESULT.md for the conclusion and PROOF.md for derivations. The
ordinary complete q-component liquid gives c<=q=d-1, with spin-charge
separation 1+(q-1). The hole is not an extra density. A crucial scope
correction is that literal projected hopping plus density interactions
preserves the material color word and supplies no spin exchange.
Species-independent couplings give a spinless charge spectrum with
q^N frozen degeneracy at fixed particle number on an open chain. The
permutation-chain saturation c=q requires occupied-color exchange.
No unrestricted maximum for all species-dependent density-only models
was proved or assumed.

The number of commuting onsite charges is <=d-1 but their CFT current
algebra only gives c>=r unless the stress tensor is exhausted. Bare
Fermi-root counting is valid for quasi-free chains, not arbitrary
interacting collective sectors. Entropy plus Lieb-Robinson yields
c T_0/v_LR <=3 log d/[pi(1-epsilon)], requiring a specified common
conformal temperature T_0. No uniform lower bound on this scale is known.
The task's pi/6 entropy coefficient uses a chiral-branch sum; the usual
nonchiral sector contributes pi c T/(3v) when both speeds equal v.

The smallest dimension left open by these arguments is d=2. C63/C67/C68
do not prove a global qubit exclusion. The proposed search family is the
parity-symmetric qubit XX/YY/ZZ/DM/field bond from C67, looking for a
neutral Ising criticality inside a finite-speed Luttinger liquid.
It is speculative. No point with c>1 or thermodynamic evidence is new.
An explicit anchor is C68's (gamma,u,D,b)=(0,0.4,0.5,0), known c=1.
The candidate and mixed-coupling obstacles are specified in PROOF.md.

checks.py uses Python and NumPy only, small local matrices. At q=2,L=6,
N=3 the color/charge spectral comparison agrees to 1.20e-14. Candidate
bond algebra and a formal entropy countertest also pass. checks.json
records outputs. No cluster jobs or large computations were performed.
No files outside work_C77 were modified, per explicit task restriction.
