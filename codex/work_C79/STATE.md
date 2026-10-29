# C79 local state

[proved, new route] PROOF.md gives a general-tree source-free moment bound with a fixed Gevrey cutoff. It permits unbounded individual edge loads and pays their constants with unused source-to-birth scale gains. This avoids attempting to extend C72's fixed-moment coefficient invariant.

[proved, invariant] Each term retains original endpoint measures, actual BBF covariances and determinants, all previously created displacement factors, and externally differentiated fields at convex combinations of terminal coordinates. C74's homogeneity step keeps at most two localization derivatives per original field. Every new displacement is a linear combination of actual descendant tree arguments and endpoint-tree edges, with absolute coefficient at most one per edge. See PROOF Sections 2-3.

[proved, bound] If b moments choose one edge, their birth scales lie on its unique ancestor chain, at most two per integer scale. Their surplus is at least q(b)=floor((b-1)^2/4), with q(0)=0. This pays (b!)^2 using half the surplus. With rho=2^(-1/2), S_v=sum_(e below v)rho^(t_e-h_v) satisfies sum_v S_v<=N/(1-rho), so the remaining route sum product_v S_v^r_v is exponential in N. The full bound is PROOF (20), with no loss from the required localization cut gain.

[proved, input realization] A fixed real Gevrey-order-two cutoff gives the arbitrary-moment covariance estimate C^(b+1)(b!)^2 E_s^(a-b-1), a<=4. Exponential endpoint measures need only their original terminal-tree norm. Optional fixed conversion legs need a Gevrey moment bound, which holds for an analytic band projection with the chosen Gevrey cutoff. Merely assuming four moments for arbitrary fixed conversion legs is not enough for this method.

[proved, engineering bookkeeping] PROOF Section 9 also includes explicit engineering-derivative endpoints. With e_w engineering derivatives surviving a cut, localization alone covers max(0,m_w-e_w). The transported engineering ratio contributes e_w. Formula (27) retains both contributions once; equivalently the canonical dimension p/2-2+e_w plus the localization gain is at least p/2-2+m_w. Thus the source-free dimensional bound does not need a change of quadratic-absorption scheme. The stronger localization-only display obtained by dropping those engineering factors is not asserted.

[assumed, next required work] The reference comparison and the Gaussian-subtracted source estimates of corrected Lemma N remain separate. Full Lemma N and the task's complete proof-or-obstruction alternative remain incomplete. An arbitrary prescribed conversion leg with only four controlled moments is also outside this Gevrey method; use the actual Gevrey construction or establish another bound.

[numerical] checks.py and checks.json give passing exact finite controls, including branched convex geometry and factorial loads through 300. No full RG computation, independent review, or formal verification.

[record] Parent STATE.md is outside the task's allowed write scope. This local state is the handoff. No subagents or shell network commands were used.
