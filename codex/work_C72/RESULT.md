# C72 result

[assumed, unresolved target] C72 is not closed. I did not prove C60 (17), complete corrected Lemma N, or construct a legal tree with unbounded loads under every reorganization.

[proved, moment-first formula] PROOF.md Sections 2 through 5 give an explicit interpolation formula that combines incoming momentum derivatives with a survival-adapted Taylor decomposition before absolute values. For remainder order m<=2 and incoming moment order b<=2, the formula differentiates the current input kernel through order at most max(m,b)<=2. With inherited homogeneous factors, it has the counters

    r=max(0,m-b-a+d),
    B=b-d+r<=2,
    A=a-d+r<=max(a,2),
    J=a+b-d+r>=m.

Its local coefficient mass is bounded, and its coordinate pushforwards preserve a permitted child anchor. B counts derivatives of the kernel supplied to this step. It is not yet a bound on the original endpoint hidden inside an integrated subtree.

[proved, cancellation across survival components] For one contracted variable c and one surviving variable s, the differentiated quadratic split gives

    R^s_2 K_c + s K_sc(c,0) = R^s_1 K_c,

and with two incoming c derivatives the complete split recombines to K_cc(c,s). Expanding the separated Taylor kernels first would introduce third and fourth derivatives. These cancellations require combining different survival components, not merely simplifying each component separately.

[proved, restricted depth induction] PROOF.md Section 7 treats arbitrary-depth chains of partial Taylor projectors on nested coordinate subsets of one fixed kernel. The chain reduces inductively to at most two projectors. With b<=2 incoming derivatives, it has an interpolation using derivatives of that kernel through order at most four, uniformly in chain length. This does not include intervening Gaussian contractions or the assignment of gains to distinct scale cuts.

[proved, restricted obstruction] PROOF.md Section 6 shows that keeping a selected surviving factor explicit can prevent a bounded two-derivative representation. For

    Q_s K=[K(c,s)-K(c,0)]/s,
    K_N=N^(-2) exp(iN(c+s)),

all derivatives of K_N through order two have magnitude at most one, but partial_c^2 Q_s K_N(0,0)=-iN. Dropping the explicit factor gives the bounded representation K_cc(c,s)-K_cc(c,0). This rules out a specified fixed-factor normalization strategy. It does not refute the exponential endpoint norm or any reorganization allowing different carrier labels.

[assumed, remaining step] The unsolved induction must reunite the required survival components after contractions while retaining their distinct ancestor operations, spatial paths and anchor constraints. The local counter B<=2 and the fixed-kernel flag induction do not establish that closure. Section 8 states the missing assertions. Reference comparison and the Gaussian-subtracted source estimate also remain separate requirements for full corrected Lemma N.

[numerical] checks.py contains passing exact rational/integer checks of the displayed formulas, local counters, cross-component cancellations and restricted flag induction. They are not a GN routing test or an interacting RG calculation. Counts are in checks.json.

[record] PROOF.md, PROGRESS.md, STATE.md, checks.py and checks.json are confined to work_C72. No other files were modified. No subagents, large computations or independent review were used.
