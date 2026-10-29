# C79 result

[proved, general-tree moment estimate in the stated scheme] The coefficient-moment problem has a general-tree bound using surplus scale gains. The argument permits several quartic children at a vertex and retarded endpoints of arbitrary degree. It does not require a fixed upper bound on every edge load.

[proved, mechanism] A moment born at scale u and routed to an actual tree edge at scale t>=u supplies an extra factor 2^(-(t-u)), beyond its birth-to-sink localization gain. At most two moments are born at each scale on a given edge's ancestor chain. A load b on that edge therefore has surplus at least

    q(b)=floor((b-1)^2/4), b>=1, q(0)=0.

Half this surplus pays (b!)^2. The other half sums every routing choice. The resulting total is C^(sum_i p_i+n), independent of depth. [PROOF.md Sections 2-6](PROOF.md) give the exact field/kernel induction, shifted-coordinate routing, cut record, and full moment estimate.

[proved, cutoff qualification] The argument uses the all-order bound

    || |X|^b partial^a g_s ||_1
       <= C^(b+1)(b!)^2 2^[s(a-b-1)], a<=4.

Section 8 derives it for a fixed real Gevrey-order-two cutoff and the analytic lattice dispersion. Exponential terminal-tree endpoint measures satisfy the needed endpoint bounds. Conversion legs must have the corresponding all-order control; merely knowing their first four moments is insufficient for this proof.

[proved, scope of (17)] This proves the source-free localization estimate including hidden moment constants for general GN trees. The localization-only form is (20). With explicit engineering-derivative endpoints, (27) retains their transported engineering dimensions once, and (28) gives the same required renormalized power count. Determinants remain Gram determinants, and coordinate weights are on actual tree-line arguments. A stronger localization-only inequality that drops those engineering factors is not asserted.

[assumed, task not fully closed] C61's N-comparison and N-source are not proved, so full corrected Lemma N is not concluded. Its source difference remains D_O(u)-D_O(0), with the free descendant separate. No first branching pattern that violates the endpoint bound was found. These limitations are stated in PROOF.md Section 9.

[numerical] checks.py passes finite exact-arithmetic controls of the geometric induction, load surplus, factorial absorption through load 300, route summation, local homogeneity, and branching cut identities. checks.json records the counts. These checks do not establish an interacting infrared limit.

[record] PROGRESS.md and STATE.md contain the handoff. All writes are confined to work_C79. No subagents, shell network commands, or independent review were used.
