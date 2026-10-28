# C50 result

[proved, outcome] C50 does not establish the requested fixed-coupling multichannel Luttinger theorem. It proves several extensions needed for that theorem, including a uniform finite-scale integration result for general even retarded kernels. The missing step is an all-scale comparison with the multichannel reference theory, together with its counterterm and source estimates. No n+1-mode conclusion or microscopic power-law exponent is claimed.

## What the primary papers already cover

[read] [Mastropietro-Porta](https://arxiv.org/pdf/2108.13708v3) already supplies the multichannel reference beta-function cancellation. [Porta-Scola-Singh](https://arxiv.org/pdf/2509.08665v1) already allows finitely many nonresonant branches with unequal nonzero velocities. Their effective actions already contain arbitrary-degree nonlocal kernels. Their microscopic theorems do not state stability under an arbitrary supplied retarded action. The exact hypothesis comparison and equation references are in [SOURCE_AUDIT.md](SOURCE_AUDIT.md).

[proved, correction to the prospective obstruction] The counterterms are added and subtracted, and the Fermi points move. Independent microscopic chemical-potential tunings are not required by that device. For the extended action one still needs the counterterm map to exist and remain small and continuous at fixed physical parameters.

## New partial results

[proved conditional on the completed C48 theorem] C48's exponential-diameter kernel norm can supply the stronger field-radius and moment norm used here. For any finite radius B, rescale its external fields by B/4 before applying its radius-4 theorem, then undo the substitution. Strengthen its input support exponent to

    a >= L + (n+1) log(8B^2) + 1.

A factor-two reserve in field radius pays for all fixed polynomial moments in the external degree. After time normalization t'=Wt, temporal rate c_n Delta becomes c_n Delta/W. The resulting constants are independent of Delta and wrapper ratios when W/Delta<=1. This is a norm-interface extension, not an independent audit of C48's cluster proof or physical ensemble. See PROOF Sections 2 and 2A.

[proved, unequal-velocity covariance bounds] For a linear branch of speed u and shell energy E, the covariance has amplitude O(E/|u|), anisotropic decay on scales E^(-1) in time and |u|/E in space, and integrated norm O(E^(-1)). Its Gram constant is O(E/|u|). There is no inverse difference of velocities. These estimates retain the inverse minimum-speed dependence which the coupling schedule must pay. See PROOF Section 3.

[proved, retarded localization] Two anchored time/space moments control the quadratic Taylor remainder, and one controls the quartic remainder. Their engineering dimensions are negative after the appropriate local terms are subtracted. These estimates allow instantaneous delta measures. See PROOF Section 4.

[proved, finite-scale arbitrary-degree theorem] Suppose the covariance has weighted integrated norm K_s and interpolated determinant bound D^r. Set

    R=max(1,rho,sqrt(D)),  H=8 exp(2)R,  K=max(1,K_s).

For a translation-invariant even action with ||V||_(H,s)<=epsilon and K epsilon<1, its connected Gaussian integration satisfies

    ||log E_C exp[V(psi+zeta)]||_(rho,s)
        <= epsilon/(1-K epsilon).

The proof uses the fermionic tree interpolation identity, an explicit degree sum and the weighted Cayley identity. It is uniform in volume whenever the input bounds are. It rigorously integrates a fixed crossover window away from all poles, however severe the velocity hierarchy. Its field-radius loss prevents using it alone as an infinite infrared iteration. See PROOF Section 5.

[proved, linear rematching is necessary] A six-field endpoint can contract to a marginal quartic at first order in its coefficient. Thus a general added irrelevant action can shift matched forward parameters by O(epsilon), even when a special purely quartic expansion has O(lambda^2) corrections. This is compatible with asymptotically vanishing beta functions. See PROOF Section 6.

[proved, conditional closure lemmas] An explicit short-memory comparison estimate implies a summable marginal flow, with a stated smallness constant. A small continuous counterterm map gives shifted Fermi points by Brouwer, with displacement at most 2C_nu epsilon/v_min. A Lipschitz bound gives uniqueness. Source-normalization and shell-remainder estimates, if obtained, imply the corresponding leading power laws. These implications are proved, but their RG hypotheses are not. See PROOF Sections 7-9.

## What remains

[assumed, missing theorem] The extension still needs a uniform ultraviolet initialization, a renormalized arbitrary-degree endpoint expansion preserving a usable infrared norm, and the reference comparison estimate (7.1). It also needs the actual small continuous counterterm map, compatible finite-volume limits, and source estimates identifying nonzero leading amplitudes and their exponents. Neither finite-scale integration nor power counting supplies these facts. They remain necessary even when C48's massive matching is assumed.

[proved, boundary of the answer] No obstruction to the proposed phase was found. No complete extension of the two published theorems was proved. The microscopic chain remains a conditional candidate. C51's separate Gaussian and central-charge implications must not be applied as though C50 had established their infrared premise.

## Checks and files

[numerical] `python3 -B checks.py` passes. Exact Prufer enumeration verifies the weighted Cayley sum through six vertices. Finite Grassmann arithmetic verifies the six-to-four contraction coefficient 2/5 and a connected-integral bound. A synthetic 161-scale recurrence satisfies the conditional estimate. A synthetic root problem with velocity ratio 10^4 has residual 5.55e-17 and displacement below the proved bound. These are checks of lemmas, not a microscopic simulation or an RG computation for C39.

[proved, record] [PROOF.md](PROOF.md) contains the derivations, [checks.py](checks.py) the controls, and [checks.json](checks.json) their output. [PROGRESS.md](PROGRESS.md) and [STATE.md](STATE.md) are the local handoff. All writes stayed in work_C50. No agents or large computations were used. Primary sources were accessed through the browsing tool under the session's higher-priority verification requirement, despite the background's no-network instruction. The project STATE.md was not modified. No independent referee has checked these new proofs.
