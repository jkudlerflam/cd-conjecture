# C58 result

[proved conditional on the imported finite construction and C48] Estimate (c) closes in its thermodynamic form for the original wrapped C39 chain, d0=4464. At each fixed n and sufficiently small fixed positive couplings, the original rooted pressure and the pressure of the exact matched thermal logical action obey, uniformly in a real length-source interval,

    |p_original(T,s)-p_matched(T,s)|<=C T exp(-delta/T),
    delta>0.

The same activated accuracy holds after conversion to fixed physical length, both with marker number summed and in canonical sectors M_L/L tending to an interior selected density. The physical cyclic identification and periodic boundary condition are included. The order of limits is L to infinity first, then T to zero. Full derivations and the hypothesis ledger are in [PROOF.md](PROOF.md).

## Gaps closed here

[proved] The collaborator's full-space matching applies to a positively completed bond, not directly to the original conditional pins. The completion changes no selector-legal dynamics. Each illegal ordered cell has an invariant mismatch flag, costs at least U+|e|ell, and can remove pins only on itself and its successor. Cutting the weak interactions touching those cells and comparing with a compression to the complete code gives

    1<=Z_original/Z_selector-legal
       <=[1+A_def exp(-delta_def/T)]^M,
    delta_def>=U/4.

Since Z_selector-legal<=Z_completed<=Z_original, this transfers the completed matching theorem back to the original bond. No pin norm is incorrectly counted as a small perturbation, and unbounded multi-token readout is retained as a positive invariant onsite term.

[proved] The existing pinning result sums marker numbers. For a fixed-M physical sector a new isometry gives, with ||h||<=J,

    (L/M)Z_pin<=Z_phys<=L exp(4 beta J L/M)Z_pin.

Adding the bounded rooted seam changes the free energy per physical site by at most

    (B+4J/rho_L)/L+T log L/L,       rho_L=M/L.

This vanishes at positive limiting density. It handles finite-wrapper length-dependent interactions, without imposing an incorrect independent periodic fermion trace.

[proved] A block-cut and positive-convolution argument supplies canonical ensemble equivalence for the all-length rooted spaces. For rho=M/L in the controlled density interval,

    f_original(T,rho)=e+sup_s[-rho p_original(T,s)-s].

Uniform pressure errors therefore transfer by the elementary supremum inequality. C44's reference length convexity and the finite norm-error budget give a strict finite-difference saddle bracket for a nonempty interval of the scheduled densities. This confinement does not require an infrared susceptibility theorem. The proof does not extract fixed-length coefficients from a bound on their generating function.

## Imported interfaces checked

[proved application] The completed all-string cell gap and heat trace, positive token-readout split, C45 local SW identity, C49 physical-support bound, charge-preserving separate Jordan-Wigner representation, and C48 vacuum-annihilating coefficient norm satisfy the needed hypotheses after strengthening the support exponent and choosing lambda last. The unbounded length charge stays in the auxiliary semigroup. The collaborator's tilted embedding and holomorphic matching estimates give a uniform source window and fixed-order endpoint derivatives. The full scalar port energy is included when bracketing the physical pressure zero inside this window.

[proved, thermal normalization] The matched action retains the auxiliary scalar logarithm, Gaussian determinants and the exponentially small finite-temperature massive endpoint. Its temperature-independent endpoint differs from the full second-order action by O(lambda^3), which is kept as an interaction. That fixed error is never discarded as o(T^2). The collaborator's gapless winding theorem is a single-scale fixed-history result; the full thermal-history estimate remains part of (b).

## Final conditional theorem

[assumed, only remaining infrared inputs] Assume (a) the all-scale matched-action RG, counterterm and length-source estimates, yielding positive dressed chiral speeds, a differentiable ground pressure and a strictly varying selected mean length; and (b) the subtracted scalar thermal estimate uniformly on the same small source interval, with the actual thermal matched action and all scalar normalizations. These are stated explicitly in PROOF Section 9. A source-free interpretation of (a) alone would not suffice for the fixed-density conclusion.

[proved conditional on (a),(b)] The original fixed-length periodic chain has

    f(T)=e_phys-(pi/12) sum_(a=1)^(n+1)
                 [1/u_(a,R)+1/u_(a,L)] T^2+o(T^2).

Here u_(a,R/L)=ell_0 v_(a,R/L) are dressed physical speeds, with ell_0 the ground-state mean cell length. For the marker-summed ensemble ell_0=-p_0'(s_*) at p_0(s_*)=0 and e_phys=e-s_*. For fixed marker density rho, ell_0=1/rho and

    e_phys(rho)=e+sup_s[-rho p_0(s)-s].

The corresponding coefficient is evaluated at -p_0'(s_rho)=1/rho. The matched and physical free energies differ by an activated error, so they have the same ground energy and T^2 coefficient. C51 gives c_R=c_L=n+1 for the complete compact Gaussian theory. Thermodynamic transfer by itself does not construct the microscopic stress tensor or identify the full physical scaling theory.

[proved, scope] No o(1/L) spectral theorem is obtained. Bounded seam errors and the C49 finite prethermal window are insufficient for it. The free-energy alternative is the completed result. There is no unconditional c-d counterexample or claim that (a) or (b) has been proved.

[numerical] All new local checks pass: 75 quantum fixed-marker trace comparisons, compression residual 3.29e-16, 6400 literal selector moves, and nine noncommuting defect comparisons. Maximum matrix dimension is 64. See [checks.py](checks.py) and [checks.json](checks.json). No Form W many-body simulation or central-charge measurement was performed.

[proved, record] All files are confined to work_C58. [STATE.md](STATE.md) is the local handoff; the project STATE was not modified under the explicit write restriction. No subagents or network commands were used. These new proofs have not been independently or human reviewed.
