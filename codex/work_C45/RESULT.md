# C45 result

[heuristic under the explicit standard input S] Yes at the project's theoretical-physics standard. C39/C42 supply the all-n finite band input, C44 has now completed the finite wrapper audit, and the local reduction below closes the infrared scale hierarchy. For every fixed n there are fixed nonzero coefficients of the repeated d0=4464 bond for which the infrared theory has n payload compact bosons and one length boson, hence c=n+1 in the project's multi-velocity convention. Comparing these fixed-n infrared limits gives unbounded c at fixed d0.

[assumed, precise remaining input] The one continuum input is ordinary Wilsonian decoupling and Gaussian RG stability for a nonresonant critical subsystem coupled weakly to massive local degrees of freedom. Its quantitative statement, including smooth filling shifts, frequency moments, the collective irrelevant norm and uniformity as W/Delta decreases, is S in [PROOF.md](PROOF.md), Section 5. This is a controlled physics inference using that input, not a rigorous microscopic CFT theorem. No microscopic central charge was measured.

[proved, qualification about C40] C40's exact connected-action hypothesis is not proved unconditionally here. Its stronger all-order construction is unnecessary for this physics conclusion. An exact finite-order local change of basis supplies the required microscopic smallness, and reduces the remaining low-frequency question to S. Finite-cluster coefficients alone would not do this.

## What closes the microscopic remainder

[proved] In ordered cell coordinates the full retained manifold has local dimension 2^(n+1), consisting of the two lengths and the n-flavor code. Write

    H=A+lambda V+lambda^2 F,       A_j>=Delta(1-P_j).

V is the compensated port. F includes marker motion, the readout and the scalar chemical adjustments. All physical compensator and transport norms are retained. Delta is the minimum complete wrapper exclusion scale, including its small charge tunneling splitting.

[proved] A range-three generator S_loc=lambda S1+lambda^2 S2 gives the exact identity

    exp(S_loc) H exp(-S_loc)=D+R,
    ||R||_a <= C_SW lambda^3,
    P D P=lambda^2[PFP-PVQA^(-1)QVP].

The last expression contains the complete finite-g second-order bond. D is locally block diagonal. The spatial interaction norm and C_SW are independent of ring size. PROOF Section 3 derives the estimate by summing the BCH series of this fixed generator, without summing an all-order SW series or isolating an extensive code band.

[proved, explicit constants] With ||Phi||_a=sup_j sum_(X contains j) exp(a|X|)||Phi_X||, put

    v0=max(1,||V||_a), f0=||F||_a,
    G=f0+48v0^2/Delta,
    C_SW=2^16(v0 G/Delta+v0^3/Delta^2).

The displayed bound holds for lambda<=min(v0/G,Delta/(192v0)). One can take a=16(n+1) to pay for conversion to complete retained coefficient norms and spatial moments. On malformed internal words the potentially unbounded per-cell readout is kept as a positive invariant-sector contribution, as explained in PROOF Section 3, rather than incorrectly counted as bounded.

[proved] A product completely positive extension of the code Hamiltonian gives

    D >= E_code+(Delta-2||D-A||_loc) sum_j(1-P_j).

At the chosen lambda this leaves at least Delta/2 per auxiliary defect in the reference problem. This proof allows defect motion and does not require D to commute with each individual defect projector. Physical length constraints are handled by the positive length cost and the tangent chemical potential, with a smaller reserved defect margin. The marker-free sector remains nonnegative, while the selected phase retains negative physical energy density. See PROOF Sections 2 and 4.

[proved conditional on S] After massive matching, subtracting all propagation reducible through retained fields, the connected low-frequency remainder obeys

    ||delta S_R||_LE <= A_dec,n C_SW lambda^3,
    ||delta S_R(omega)-delta S_R(0)||_LE
          <= A_dec,n C_SW lambda^3 |omega|/Delta.

The norm sums connected spatial kernels, all external legs, two time moments and filling derivatives. It refers to integrating the massive complement, before the subsequent gapless band matching. This is the concrete low-frequency statement needed from S; it is not a bound on a bare off-code correlator at all energies.

## The coupling schedule

[proved, inherited normalization] The correct payload unit is

    W=lambda^2 J0/s,
    J0=2g^2 q^(2n+2)D^(-2n-4)/S_q^2,
    S_q=sum_(j=0)^(n+1)q^(2j).

The propagation and visibility factors cannot be dropped. In particular, equalization must control the full absolute bond error below J0 times an infrared tolerance.

[proved conditional on S, ordered choices] PROOF Section 8 gives the complete inequalities. The order is as follows.

1. Obtain C42's rational witness at n, including the momentum margin, physical minimum speed and simple-root patches. Choose an interior length density avoiding the finitely many mixed resonances. Select epsilon_n below the root-shift, stiffness and pre-detuning RG tolerances. The matching constants may depend arbitrarily badly on n and the velocity ratios.

2. Fix finite D using C42's explicit Cauchy bound with the smaller tolerance. Choose nonzero g so that both the finite-g quadratic error and

       3(1+eta)(2+64n^2) g^2 D^(2n+4)

   fit the static epsilon_n budget. The displayed quantity bounds the full normalized non-Gaussian second-order error. The reference momentum shape is held fixed while g decreases.

3. Choose the wrapper vacancy ratio w_c/U small enough, then s/delta_fast small enough, that all four complete bonds satisfy

       max_(p,q)||s K2_(p,q)-K_bare||
                           <= J0 epsilon_n/(16 C_b),
       C_b=16(n+1)^2 exp(48(n+1)).

   C44 supplies explicit source norms and the nonasymptotic charge bound

       delta_ch,p >= w_c [ (w_c/U)/(k+p-1) ]^(2kp)
                              /[binomial(k+p,p) kp],
       k=n+2, p=1,2.

   Thus this step does not rely on an unquantified tunneling asymptotic. Recompute the complete positive Delta after fixing the wrapper.

4. Set the physical marker hopping to t_m=W/beta_n, with C44's positive identity-transport overlap beta_n. Set the readout coefficient to lambda^2 J0 mu S_q/(s q^2). Retune the marker chemical potential, including the full common scalar port energy, to select the chosen length density. Include all these physical norms in F.

5. Finally choose one positive lambda with the local SW restrictions and

       lambda <= epsilon_n (J0/s)/(16 A_dec,n C_SW).

   Also require the retained bandwidth and tangent energy below the complete auxiliary gap, and the per-physical-site perturbation below the marker-free energy-density margin. PROOF (35) gives these extra positive bounds. Then

       A_dec,n ||R||_a/W <= epsilon_n/16.

   No parameter depends on L.

[proved] This schedule pays the full loss s/J0 in the remainder comparison. It needs no C26 rho^2 suppression. The fourth-order reflection-cell result cannot be imported into the two-bus wrapper, and was not used.

[proved] Scalar equalization error gives a genuine length-density coupling. If each scaled length-pair error is at most epsilon, its coefficient is at most 4epsilon and the corresponding traceless mixed term at most 8epsilon. Those errors are included above. One-port/one-marker terms vanish by vacuum-versus-bus conservation. The surviving higher marker, mixed and readout terms are included in the local R bound. Reducing lambda controls them after all earlier choices are frozen.

## Why control persists to 1/L

[proved within the compact-boson dictionary] Exact payload U(1), marker number and identity transport leave only finitely many dangerous dimension-at-most-two vertices. C42 detunes the payload ones. An interior density can detune the mixed length-bilinear ones with a positive finite margin. All remaining non-Gaussian primaries have free dimension at least three. A relative stiffness deformation at most 1/4 keeps their dimensions above 3sqrt(3/5)>2, even with arbitrary velocity ratios.

[heuristic under S] Keep every retained band while matching through the hierarchy. With e_det the normalized lowest linearization/detuning scale and ell_*=log(e_top/e_det), choose the collective initial coupling so that A_beta X0 ell_*<1/2. Below W e_det only Gaussian marginal terms and irrelevant interactions remain. Their feedback gives finite stiffness and filling shifts, while the irrelevant norm decays as a positive power of E/(W e_det). The slowest mode therefore remains gapless at arbitrarily low energies, including the finite-size scale.

[proved, interpretation] A fixed extensive error norm cannot be smaller than 1/L for all L. That is the wrong criterion. Gaussian corrections renormalize velocities by a small fixed relative amount. Irrelevant corrections to the dressed finite-size spectrum vanish relative to its 1/L spacings. The thermodynamic/infrared limit is taken separately at each fixed n, before n is increased.

## What fails, and what was checked

[proved] The tempting shortcuts fail: a global-band SW comparison, using U as the wrapper gap, comparing equalization only with g^2, dropping scalar length interactions, using lambda to fix the second-order conversion error, importing C26's special suppression, or summing C40's fixed-order coefficient bounds. The argument above uses none of them. No failure of the actual scale schedule was identified.

[proved, completed neighboring input] C44 finished during this calculation. Its final RESULT.md and PROOFS.md establish the repeated 4464-letter bond, all-length/illegal-sector bounds, complete second-order equalization and identity-transport gluing for every n. C45 read that completed proof and incorporated its explicit constants. C44's diagnostic parameters do not meet the lacunary infrared budget; their numerics are not evidence that those particular couplings realize c=n+1.

[numerical, new controls] `python3 -B checks.py` checks the new local SW algebra on complete two- and three-cell matrices of dimensions 16 and 64, with code dimensions 4 and 8. The second-order identity error is below 6.6e-17. At lambda=.02,.01,.005, the remainder divided by lambda^3 approaches .341 and .621, respectively. The defect-threshold inequality passes to 3.3e-15 roundoff. These are independent algebra controls, not microscopic wrapper or central-charge computations. Exact rational budget checks also pass.

[proved, limits of the conclusion] No exact all-order massive integration, common Lorentz velocity, every-ring ground-state uniqueness, full microscopic many-cell ED or measured central charge is claimed. S remains an assumed standard physics input. Subject to that explicitly stated input, the infrared step now closes and the all-n fixed-d0 conclusion is established at the theoretical-physics standard.

[proved, files and scope] [PROOF.md](PROOF.md) contains the derivations, [checks.py](checks.py) and [checks.json](checks.json) the small controls, and [PROGRESS.md](PROGRESS.md) and [STATE.md](STATE.md) the handoff. No network, subagents, cluster computation or edits outside work_C45 were used. The project STATE was left unchanged because it is outside the permitted write scope.
