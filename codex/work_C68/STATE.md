# C68 local memory

2026-10-28. Task complete as a qualified assessment, not a proof of the requested global premise.

[heuristic] RESULT.md explains why the full gapless set in the T-only, reflection-only, and symmetry-free residuals need not have codimension one. A weak twisted-XXZ construction gives a physics-level open c=1 basin in all three strata. Untuned allowed commensurate relevant operators do yield local gapping arguments, but cannot establish global genericity of a gap.

[assumed] The construction uses the complete XXZ bosonization and weak oscillatory-perturbation RG. Anchor S=diag(1,1,0.4), D=(0,0,0.5), B=0. K=0.8110785334459236, v=5.444225545344561, alpha=0.46364760900080615. Local vertices exp(i m theta+2i n phi) have dimension m^2/(4K)+n^2 K and momentum m(pi+alpha)+n pi. All five conjugacy classes with dimension at most two are nonresonant. First neutral umklapp has dimension 3.2443141337836945. No rigorous perturbation radius or new thermodynamic simulation is supplied.

[proved] Three explicit small-eta perturbations in RESULT.md lie in the requested strata. Gapped anchors and bounds are inherited from C34. An Ising point has one relevant tuning only when its magnetic field is forbidden and its thermal matching coefficient is transverse. Otherwise there are two, before any mixed-sector operators are counted. Kramers LSM on T requires broken symmetry in a gapped phase or anomaly matching at a symmetric critical point. It does not upper-bound c.

[heuristic] c>1 remains symmetry-permitted with fine tuning. The main candidate mechanism is an independent neutral Ising sector in a floating boson, with all invariant relevant products and marginal mixed terms analyzed. No supplied point is a supported c>1 candidate. Point 1282 is a calibrated c=1 starting point in N. RESULT.md gives its bond and the continuation J=J0+s e_z e_z^T, B=B0+t e_z, without predicting a critical location.

[numerical] checks.py and checks.json contain only deterministic symmetry, bond-rotation, normalization, and continuum vertex checks, plus one extracted parameter record. All assertions passed. The c=1 data quoted for 220 and 1282 are inherited from the local follow-up README and were not rerun.

[proved, file scope] All work and Python execution stayed here. No network command or large computation was used. The project-wide STATE.md could not be edited within the authorized write scope. This file provides the handoff instead. PROGRESS.md records completion and RESULT.md contains the full argument, sources, limitations, and numerical proposal.
