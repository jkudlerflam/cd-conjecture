# C89 handoff

[Status] Read RESULT.md, then PROOF.md. S_pin remains open for the full finite-K d=13 mobile chain. Do not present this folder as a proof of mobile sector selection.

[Proved] Finite logical zero-port vacuum-clock pin exclusion follows by replacing every pinned terminal clock vector by a terminal-plus-penultimate superposition. The energy gain is J times the sum of terminal contact expectations. Perron-Frobenius makes these strictly positive at finite size. No uniform thermodynamic contact lower bound was obtained.

[Proved] In the legal projected chain the one-marker sector obeys E>=E_gap-4J. A selected identical-cell vacuum product with the free bath gives E_good<=E_gap+J a_L N e15. At half occupancy, compatible N>=2865 excludes that entire giant-cell sector, including auxiliary states. This is not an island comparison or a finite-K result.

[Proved] The legal marker-free words are all-L or all-R. Their projected material transitions vanish and their nonnegative field has minimum zero. Hence E_marker-free=E_gap, strictly above the selected trial on compatible N>=15 with positive contact. This does not include finite-K dressing.

[Reference constants] a_half=.29735763271532445, e15=-.004825498950096, pin bare cost=1.96157056080646, leading mobile pin cost=.583287978365481 J. d14=d16=.000449860388753, giving .000133769420252 J per wrong cell at half occupancy. Use the inherited C88 selector certificate for rigorous constants. Log2(lambda^2 Eunit)=-2074.577410080, log2(field)=-2194.577410080.

[Missing input] A one-sided common-bulk residual inequality R_sigma-R_good>=-J epsilon(J/t) W, epsilon->0, for quantum sources and arbitrary lengths. W=sum_wrong lengths+number of pinned good-length cells. The UV first-moment normal form is insufficient because its low-frequency term still needs a ground-state comparison. Finite-K local strengths also need defect matching before they can be charged to W rather than total volume. The scalar-contact continuum argument does not establish this operator-valued inequality.

[Scale bookkeeping] Finite-K local corrections divided by J contain 2/[(J/t)^2 K]+3/K+6lambda^2/K plus C3(nu/t)^3/[(J/t)^3 K^2]. The material lower bound K>=2^2560 cannot be held fixed while sending J/t arbitrarily small. C88 leaves J/t, the mobile K constant, and the infrared tolerance implicit.

[Numerical] check_defects.py gives 24 exact vacuum-clock and 24 scalar-contact comparisons, all positive. Clock dimensions <=6750, scalar <=6435, maximum eigenpair residual 8.89e-12. b15 quantum cases are N15,V2 and N30,V1. Half-density quantum cases use b3,N6,V6. Half-density scalar cases N6,8 insert b15 reference costs but do not model b15 quantum geometry. All use lambda=f=0,K=infinity,t=1, with J=.04,.02,.01 or .02,.01,.005.

[Numerical] second_order.py evaluates the full Gram matrix of h_j acting on clock ground states and contracts it with the finite bath susceptibility. Maximum linear-solve residual 2.49e-13. It explains the observed O(J^2/t) corrections. Coefficients and exact energies are in second_order_results.json. No infinite-volume remainder bound is inferred.

[Files] PROOF.md contains derivations and the conditional budget. RESULT.md is the short outcome. PROGRESS.md records work status. check_defects.py, summarize_checks.py, second_order.py reproduce the three JSON outputs. ed_stdout.txt is the retained first-pass console output. No external files were modified, including the project STATE.md, under the user's explicit folder restriction.
