# C50 local state

Last updated 2026-10-28.

[proved, task status] The requested constructive multichannel RG theorem is not proved. Do not mark S(b) complete. The finished output is RESULT.md, with full partial derivations in PROOF.md. All-scale marginal comparison, counterterm construction and source asymptotics remain open.

[read, important source corrections] MP 2108.13708v3 Theorem 8.5 already treats multichannel reference beta cancellation. PSS 2509.08665v1 Assumption 2.1 allows arbitrary nonzero velocities, without a velocity-separation condition. Its equations (4.48)-(4.53) add and subtract counterterms at fixed Hamiltonian. Both generate arbitrary-degree nonlocal kernels internally. SOURCE_AUDIT.md records the actual comparison. Do not claim that many flavours, unequal speeds, or the mere presence of retarded kernels are absent from those works.

[proved, new finite-scale theorem] For covariance integrated bound K_s and interpolated determinant bound D^r, choose H=8 exp(2) max(1,rho,sqrt(D)), K=max(1,K_s). Input coefficient norm epsilon at H gives connected output norm <=epsilon/(1-K epsilon) at rho if K epsilon<1. Proof uses BBF interpolation, explicit endpoint sums, weighted Cayley. This does not iterate to zero energy without additional RG estimates.

[proved conditional on C48] Its completed enlarged-rooted-model theorem can be reweighted to radius B via external substitution psi=(B/4)chi. Recompute its input support exponent a>=L+(n+1)log(8B^2)+1. Choose B>=2H to pay for degree-dependent moments when converting exponential diameters to polynomial tree length. Rescale time by W. Constants then stay uniform for W/Delta<=1. C48's theorem and physical-sector premises are imported, not independently verified.

[proved, other lemmas] Retarded Taylor localization with two moments; velocity-sensitive continuum covariance bounds; O(epsilon) six-to-four endpoint contraction; conditional summable-flow inequality (7.1)-(7.4); continuous counterterm inversion using Brouwer and optional Lipschitz uniqueness; conditional source normalization and shell summation. No rederivation of the inherited C39/C42 momentum certificate or C45 Gaussian stiffness bound was needed.

[assumed, next concrete step] Prove (7.1) in the explicit kernel class after localization and dimensional rescaling, including degree sums, finite-volume bottom scales and source insertions. The finite-scale lemma loses field radius and cannot substitute for this. A uniform Matsubara ultraviolet initialization is also needed if starting from a Hamiltonian or an uncut action.

[numerical] checks.py passes exact combinatorial/Grassmann controls and synthetic recurrence/root controls. checks.json records them. No microscopic central charge, thermodynamic extrapolation or actual logical RG trajectory was computed.

[proved, scope] All writes are in work_C50. No subagents. Primary papers were read using the session-required browsing tool; local commands used Python 3 only after reading the background. Project STATE.md is outside the permitted write scope and remains unchanged. C48/C49/C51/C52 progress was reread near completion to avoid duplicating their tasks.
