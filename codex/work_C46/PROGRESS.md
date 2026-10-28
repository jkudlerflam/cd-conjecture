# C46 progress

[proved, work record] Read project STATE, STATUS, the infrared literature and leakage reports, C45 RESULT and reduction proof, C40 kernel definition, C42 finite band proof, and C44 wrapper bounds. No preexisting C46 progress file was present. All writes stay in work_C46. The parent project STATE will remain untouched under the task restriction.

[assumed, inherited inputs] C39/C42 finite-band certificates and C44 finite-wrapper estimates are inputs, with their proof locations to be recorded. Their all-order infrared interpretation is not assumed as a conclusion.

[proved, located obstruction] The proposed all-temperature, coupling-small action theorem already has the explicit thermal counterexample in astra_B Section D. A cell gap does not give a Gaussian logical thermal state at zero coupling. Any replacement must specify thermal winding terms and which logical observables it matches.

[heuristic, current route] Investigate a uniform partition-function comparison for the locally block-diagonal reference, an explicit failure of absolute Grassmann coefficient norms for generic logical-dependent excited cells, and the precise extra RG and thermodynamic inputs. The goal is a rigorous partial result with no unsupported claim of microscopic central charge.

[assumed, source limit] Only local literature summaries have been located for 2108.13708 and 2509.08665. No network will be used. Their exact hypotheses and conclusions must not be asserted as independently verified primary-source facts.

[proved, new result to be written out] A trace comparison follows from the product completely positive extension used in C45. For finite cell dimension d_P+d_Q and locally block-diagonal D=A+W, with J=sup_j sum_{X containing j}||W_X|| and delta=Delta-2J>0,

    1 <= Z_D/Z_P <= (1+(d_Q/d_P)exp(-beta delta))^M.

The proof uses scalar convexity on the partial trace of the logical Hamiltonian. It allows moving and changing defect patterns inside local Q blocks. For an infinite direct sum, d_Q exp(-beta delta) is replaced by the one-cell excited heat trace.

[proved, new obstruction to be written out] In the one-cell model q(Delta+u n_1 n_2)+lambda sigma_x, a raw vacuum auxiliary excursion produces the factor exp[-Delta t-u integral_0^t N(s)ds]. The absolute integrated coefficient of its m-th external quartic insertion is lambda^2 u^m/Delta^(m+1). Its field-degree sum diverges for u>=Delta, although the cell is gapped and has exact logical U(1). This is a failure of the raw external-field expansion, not a proof that every resummed low-frequency action is impossible.

[proved, normalization correction] For a quadratic reference with 2(n+1) simple Fermi roots, the thermal coefficient is pi/12 times the sum of inverse absolute chiral speeds. Equal speeds within each left/right pair are needed to rewrite this as pi/6 times a sum over n+1 mode speeds. Neither the root certificate nor exact U(1) supplies that equality.


## Completed work

[proved, theorem T] PROOF Section 2 proves 1<=Z_D/Z_P<=product_j(1+t_j/d_j), where t_j=Tr_Q exp[-beta(A_j-2J_j)]. The result permits defect hopping and changing defect number. C44's linear length penalty gives an explicit C_th and beta_0 for the rooted direct sum. Combining the theorem with C45's imported exact local reduction gives the stated pressure error C_SW lambda^3 plus an exponentially small thermal term. It does not control the T^2 coefficient at fixed lambda.

[proved, action obstruction] PROOF Section 3 proves both the finite-temperature non-Gaussian counterexample and a zero-temperature divergence of the raw absolute external-field expansion at second order. The bounded operator effective Hamiltonian is computed in the same example, showing that this is an obstruction to the raw expansion rather than a no-go for resummed effective actions. Source insertions and thermal winding terms are explicitly distinguished.

[proved, finite matching] PROOF Section 4 gives the complete second-order operator-valued time-moment bound using C44's source constant B_src=3(1+eta)V0(1+M_h^4). Section 5 proves the exact charges and, conditional on the imported C39 root construction, a terminating optional zero-pair test and an explicit full-channel density margin. Distinct absolute chiral speeds are not falsely inferred from transversality.

[proved, reference thermodynamics] PROOF Section 7 derives the free thermal coefficient from simple roots and from Gaussian oscillators. It proves canonical free-length ensemble transfer, bounds the common-bond cyclic quotient by factors b exp(-4pi beta |J_len|) and (b+1)exp(4pi beta |J_len|), and transfers the selected Gaussian reference to variable marker count and physical length. It retains separate left/right speeds. Section 8 proves the Gaussian central charge and the conditional unitary-current lower bound.

[assumed, unresolved] D* is the resummed massive-sector action with uniform moments, source matching and dressed sector control. R* is the constructive RG application at fixed microscopic parameters, including retarded higher-degree vertices and dressed roots. T* is scalar vacuum-functional control to o(T^2), or an appropriate complete stress/current scaling limit. None has been proved for the actual 4464-letter chain here. The exact primary hypotheses of the two requested RG papers remain unverified because only local reviews were located and network access was forbidden.

[numerical, final controls] python3 -B checks.py passes. Nine small pressure comparisons include nonconservation of defect number; the operator lower bound is consistent with zero to 2.40e-16 roundoff. The thermal counterexample and two-band chiral-speed coefficient controls pass. These are toy controls, not microscopic chain or central-charge computations.

[proved, handoff] RESULT.md and PROOF.md are complete. STATE.md is the local memory, used because the parent project STATE is outside the permitted write scope. checks.py and checks.json make the numerical controls reproducible. No network, agents, large computation or edits outside this folder were used.
