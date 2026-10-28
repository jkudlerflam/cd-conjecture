# C61 result

[unresolved] C60 (17) is neither proved nor refuted for the proposed full routing. No legal, nonvanishing single-channel counterexample tree was found. The requested proof-or-counterexample alternative remains open.

[assumed, corrected statement delivered] [PROOF.md Section 2](PROOF.md) states Lemma N with separate source-free action, reference-comparison and source hypotheses. The source difference is decomposed as

    D_O(u)=D_O(0)+[D_O(u)-D_O(0)].

Only the bracket is assigned an interaction prefactor. D_O(0) retains the full free lattice-to-reference difference, including C60's nonzero pair descendant 3/[4 pi^2 x^2(x^2-1)] at p=pi/3, x=3m. This repairs the Gaussian normalization of the proposed source statement, without claiming its nonlinear estimate.

[proved] For an oriented zero/derivative incidence record in a reduced scale tree,

    sum_f a_f s(f)-sum_l b_l s(l)
       =-sum_v F_v ell_v,

where F_v counts paired paths crossing a cluster edge. Cut coverage F_v>=z_v and a coefficient-mass bound C^(sum p_i+n) imply (17). Consumed zero/derivative pairs must be concatenated rather than erased. Section 3 gives the proof and the necessary nonvanishing qualification for a counterexample based on a deficient cut.

[proved] Section 4 gives a local Taylor normal form including all product-rule terms. With remainder order m, incoming moment count b, exposed derivative degree a and d consumed derivatives, its residual Taylor order and cut count are

    r=max(0,m-b-a+d),
    J=a+b-d+r>=m,
    B=b-d+r<=2 when b,m<=2.

This is the local algebra needed by the proposed routing. Realizing these counts simultaneously in the actual spatial tree remains unproved.

[proved, local weighted-line estimate] Once exact operator identities justify a weighted-propagator representation, retain X^b g_s(X) under its derivative and relocation operators. All its Leibniz terms have scale cost 2^[s(a-b)], relative to the unweighted line norm. Shifting the propagator alone leaves an internal offset that the anchor-Jacobian argument does not control. Section 5.1 proves this rule; Section 6 lists the remaining global construction. It is not claimed as a complete corrected routing.

[numerical, finite algebra controls] `python3 -B checks.py` passes 1,755 polynomial normal-form identities, 294 local count cases, 100 branched cut identities and 31 Taylor-integral cancellations, using exact rational/integer arithmetic. These do not verify the GN routing or an interacting RG theorem.

[record] PROGRESS.md and local STATE.md contain the handoff. All writes remain in work_C61. No subagents, large computation, or independent review were used. No Luttinger, thermal, or central-charge conclusion is claimed.
