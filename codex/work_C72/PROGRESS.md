# C72 progress

[assumed, task] Construct a simultaneous moment and survival normal form for the source-free GN/BBF tree, with bounded loads and retained paths, or prove an obstruction for every such reorganization. C60 (17) is the scale target. All writes stay in work_C72.

[proved, imported] C61 proves the local homogeneous-remainder normal form and the cut identity. C64 proves weighted-line transposition with path concatenation. C70 includes the full exterior-field derivative and splits quartic remainders according to port survival. None establishes the depth-independent simultaneous representation.

[record] Read the S and S2 backgrounds and all files in C61, C64 and C70. Read the endpoint-majorant construction and C60's conventions. No subagents or network commands are used. Local STATE.md will provide the handoff because parent memory edits are prohibited.

[assumed, investigation] Test partial Taylor projectors adapted to survival sets under incoming momentum derivatives. The immediate issue is cancellation between the derivatives of explicit surviving factors and their integral kernels. A scalar identity alone will not be counted as the required spatial induction.

[proved, written] PROOF.md now gives the moment-first flag formula with homogeneous factors, bounded local counters and a measure pushforward interpretation. It exhibits cancellations between different survival components that remove third and fourth derivatives of the original kernel.

[proved, restricted obstruction] Keeping a particular surviving factor can require a third derivative even after exact algebra. Section 6 rules out a uniformly bounded formula using only derivatives through order two for its divided-difference coefficient. Recombining or dropping that factor avoids the obstruction. This is not a no-go theorem for the GN problem.

[proved, restricted depth induction] A chain of partial Taylor projectors on nested subsets of one fixed kernel reduces to at most two projectors. With b<=2 incoming derivatives its interpolation has order at most four independently of chain length. This theorem excludes intervening contractions and does not assign spatial scale gains.

[assumed, unresolved] I have not proved closure of the moment-first representation under GN merging. C60 (17), corrected Lemma N, and the alternative universal obstruction remain unestablished. RESULT.md and STATE.md record this limitation explicitly.

[numerical, completed] python3 -B checks.py passes 42 flag identities, 42 integral-remainder identities, 330 moment-first identities, 6,222 homogeneous-factor normal forms, 110 partial commutations, 2,997 local counter cases, 18 cross-component cancellations, 100 checks of the obstruction sequence, 480 nested-flag reductions, 3,040 nested-flag moment commutations and four retained-pair interpolations. The finite flag chains extend through depth 40. All arithmetic is rational/integer. These checks do not verify the missing GN closure.

[record] Primary BM Sections 3.3-3.4 and 3.9-3.11 were checked using the browsing tool. No source text is used as a substitute for the missing induction. All deliverables and checks are local to work_C72.
