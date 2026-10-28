# C46 result

[proved, tagging convention] Each tag covers its displayed formulas and continuation prose up to the next tag or heading. Conditional results name their imported inputs.

[proved, verdict] The infrared step for the full wrapped d0=4464 chain is not closed rigorously. This work proves a uniform thermal comparison for C45's locally block-diagonal reference and the common-bond cyclic quotient, but not the exact massive-sector action, its constructive RG flow, or c=n+1 for the microscopic Hamiltonian. The stronger raw Grassmann statement cannot be inferred from a cell gap alone. Full proofs are in [PROOF.md](PROOF.md).

## New theorems

[proved, thermal reference decoupling] Let D=sum_j A_j+sum_X B_X, with A_j>=Delta Q_j and each B_X block diagonal between its local all-code space P_X and its complement. Put J=sup_j sum_(X containing j)||B_X|| and delta=Delta-2J>0. For finite identical cells,

    1 <= Z_D/Z_P <= [1+(d_Q/d_P)e^(-beta delta)]^M,
    0 <= f_P-f_D <= beta^(-1)log[1+(d_Q/d_P)e^(-beta delta)].

This holds for every M and beta>0, permits moving defects and changes of defect count, and assumes no gap inside the logical chain. Section 2 proves it using a product extension of H_P=PDP and a partial-trace inequality. Infinite cell spaces are allowed by replacing d_Q exp(-beta delta) with the excited one-cell heat trace.

[proved conditional on the imported cell/reduction bounds] C44's linear wrong-length cost gives an explicit heat-trace constant C_th in PROOF (9). C45's constants give J<=Delta/48 and ||R||_a<=C_SW lambda^3, so in its rooted ensemble,

    |f_H-f_P| <= C_SW lambda^3
          + beta^(-1)log[1+C_th e^(-beta Delta/2)]

for beta>=log(2d0)/(t_sel eta_b-U). The thermal defect contribution is smaller than every inverse power of beta. The fixed leakage error is not. Multiplying it by beta^2 diverges at fixed lambda, even after subtracting ground energies. This estimate does not determine a central charge.

[proved conditional on C44's exact reference spectrum] At physical length L=bM+r, the common-bond quotient has

    b e^(-4pi beta |J_len|) Z_product
       <= Z_physical,reference
       <= (b+1)e^(4pi beta |J_len|) Z_product.

Thus cyclic gluing leaves its thermodynamic pressure unchanged after conversion from cells to physical sites. This holds for any common translation-invariant payload bond. It does not treat the finite-wrapper length-dependent interactions. See theorem Q in Section 7.

[proved, finite-order time control] For the actual C44 pair source, let

    B_src=3(1+eta)V0(1+M_h^4).

Its vacuum excursion K(t) satisfies, for kappa<Delta,

    integral_0^infinity t^m e^(kappa t)||K(t)||dt
       <= lambda^2 B_src m!/(Delta-kappa)^(m+1).

The spatial local norm adds at most 2exp(2mu) for bond weight exp(mu|X|). This includes the complete finite-g second-order source, its scalar and mixed terms, and every time moment at that order. It is an operator-valued transfer-kernel estimate, not an all-order external-field action bound.

## Why the action theorem is still missing

[proved, finite-temperature obstruction] Astra B's example H_c=q(Delta+u n_1n_2) has a complete logical code and a positive gap, but its reduced logical thermal state is non-Gaussian already at zero coupling. Section 3 gives the full proof. A uniform all-beta theorem whose non-Gaussian action vanishes with the coupling must include thermal terms or additional factorization hypotheses.

[proved, stronger raw-expansion obstruction] In the same cell, adding lambda sigma_x gives a vacuum auxiliary excursion

    lambda^2 e^(-Delta t) exp[-u integral_0^t N(s)ds].

The absolute integrated coefficient with m external quartic insertions is lambda^2 u^m/Delta^(m+1). Its sum over field degree diverges for u>=Delta, even at zero temperature. A positive external-leg weight only worsens it. Meanwhile the operator effective Hamiltonian is the bounded polynomial

    -lambda^2/Delta
       +lambda^2 u/[Delta(Delta+u)] n_1n_2.

Thus the gap controls operator propagation but not the raw summed Grassmann norm. A resummation of logical products and contact contractions is necessary for that route. This counterexample does not refute a suitably resummed action for the wrapped chain. It identifies a step that a polymer or determinant argument must actually prove.

## RG hypothesis audit

[proved] The actual microscopic local rules preserve payload U(1), and the rooted representation preserves total length. Section 5 proves both selection rules. The missing-U(1) issue in the older literature notes is repaired for C39/C44.

[proved conditional on C39/C42's finite-band input] A strengthened finite search adds the test that no payload pair sums to zero. Its termination follows from C39's same-sign sea centers. Choosing the length density outside a finite forbidden set then gives a full pair-sum margin for all 2(n+1) Fermi points, explicitly

    delta_all >= min(delta_pay,delta_zero,pi/[8(F_n+1)]),
    F_n=3 binomial((2n+2)(2n+3)/2,2).

The extra zero-pair test is unnecessary if the RG uses separate payload/length conservation, but it allows a stronger unrestricted pair-sum hypothesis. No alphabet enlargement is needed.

[proved] C42 certifies nonzero physical speeds, not all distinct absolute speeds. Already its n=1 scalar reference has equal speed magnitudes at its two roots. The local literature summaries allow arbitrary velocities; no requirement of all distinct absolute speeds is inferred from them. Distinct Fermi points and transverse crossings are separate conditions.

[assumed, still needed for an application] The missing inputs are an exact norm-small retarded action after resummation, a proved RG initialization for its arbitrary-degree vertices, control of dressed Fermi points at the fixed microscopic parameters, and the physical ensemble transfer beyond the common-bond reference. Section 6 separates these from the finite checks. Only the supplied local reviews of 2108.13708 and 2509.08665 were available, so their full primary hypotheses have not been independently verified. No network was used and no unpublished extension is treated as a cited theorem.

## Thermodynamics and central charge

[proved, quadratic reference] For the C42 quadratic payload plus the free length band, theorem F gives

    f_G(beta)=e_G-(pi/12) sum_(i=1)^(2n+2)|u_i|^(-1) beta^(-2)
                      +o(beta^(-2)),

where u_i are physical signed chiral speeds in the chosen length units. The proof uses the exact free partition function and the simple-root density of states. Section 7 also proves the canonical free-length transfer, controls the common-bond quotient by theorem Q, and minimizes over marker count for the selected Gaussian reference. Physical speeds are cell speeds multiplied by b+nu.

[proved, Gaussian fixed point] For n+1 nonchiral Gaussian modes the coefficient is

    -(pi/12) sum_(a=1)^(n+1)
         [v_(a,R)^(-1)+v_(a,L)^(-1)].

The proposed -(pi/6)sum_a v_a^(-1) form needs equal right/left speeds in each pair. A common Lorentz velocity is not assumed. Gaussian current quantization gives c_R=c_L=n+1, hence c=n+1 in the conjecture's convention, if this is proved to be the microscopic scaling limit.

[proved, limitation] Normalized correlation functions discard source-independent vacuum factors. They therefore do not by themselves supply the thermal coefficient. A pressure comparison with error o(beta^(-2)), or a complete current/stress scaling-limit theorem, is still required. Under an independently established unitary CFT limit, n+1 positive-level chiral currents would give c>=n+1 by the fully proved Sugawara bound in Section 8. The microscopic charges alone do not establish those emergent currents.

## Remaining conjecture and next step

[assumed, remaining conjecture] At fixed nonzero couplings in the C45 schedule, the full wrapped chain has n+1 Gaussian modes and c=n+1. This remains a conjecture here, with the finite construction and physics-standard argument intact. No counterexample to that phase was found.

[heuristic, shortest path] Prove resummed low-frequency massive matching with two time/space moments and source control, including winding and excluded-sector bounds. Then check and extend the constructive RG to that action at fixed bare parameters. Retain the scalar vacuum terms and prove an o(T^2) thermal remainder, or construct the complete chiral current/stress limit. Theorem T has removed the thermal defect trace of the block-diagonal reference from this list; it has not removed leakage.

[numerical, verification] Nine small partition-function comparisons pass, including defect-changing interactions. Independent controls reproduce the raw coefficient growth and the chiral-speed thermal coefficient. They use local dimensions at most 3 on four cells and a two-band quadrature. They do not simulate the 4464-letter chain or measure its central charge. Run `python3 -B checks.py`; outputs are in [checks.json](checks.json).

[proved, files and scope] [PROOF.md](PROOF.md) contains all new proofs and the assumption ledger. [PROGRESS.md](PROGRESS.md) and [STATE.md](STATE.md) are the local handoff. All writes stayed in work_C46. The project STATE was left unchanged under the user's restriction.
