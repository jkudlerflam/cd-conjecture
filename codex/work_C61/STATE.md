# C61 local state

[status] Partial result. Full C60 (17) remains unresolved; no legal nonzero violating tree found. Do not cite this run as proving the source-free endpoint majorant.

[proved] PROOF Sections 3--5 give the cut telescoping identity, product-rule normal form and the fixed-order bound for derivatives of a weighted propagator. The normal form is partial^b R_m(P_a F)=sum_d binom(b,d)(partial^d P_a) R_max(0,m-b-a+d)(partial^(b-d) F), with its multi-index interpretation. Its algebraic cut demand J is at least m and its unallocated coordinate count B is at most two when b,m<=2.

[assumed] PROOF Section 2 states corrected Lemma N. Free source descendants require their own full Gaussian subtraction. The interaction bound concerns D(u)-D(0), not D(u). Source and reference-comparison hypotheses remain separate from the action majorant.

[open next step] Build a single exact GN/BBF normal form realizing the cut counts while preserving bounded loads, coefficient mass and permitted anchor incidence. Keep derivatives and Taylor relocations on X^b g(X) together. Do not expand Taylor integrals and then take absolute values before using remainder identities; PROOF (7)--(8) shows the spurious third moment this can create. Also retain the engineering-derivative convention in C60 (17).

[checked] checks.py passes finite exact-arithmetic algebra tests. checks.json records their counts. No full routing search or interacting simulation was performed.

[scope] Task rules prohibit writes outside work_C61, so this file substitutes for a parent project-memory edit. No subagents or independent review.
