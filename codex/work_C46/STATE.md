# C46 local handoff

[proved, completed] RESULT.md is the concise report and PROOF.md contains all new derivations. The rigorous infrared closure for the wrapped d0=4464 chain was not obtained. The exact status is stronger reference-level control with three remaining microscopic assumptions, D*, R* and T*, stated in PROOF Section 8.

[proved, new theorem] For a locally block-diagonal reference D=A+sum_X B_X, with J_j=sum_(X containing j)||B_X||, code rank d_j and A_j>=Delta_j Q_j,

    1 <= Z_D/Z_P <= product_j[1+Tr_(Q_j)e^(-beta(A_j-2J_j))/d_j].

The proof uses a normalized partial-trace extension and scalar convexity. It needs no logical gap or individual defect conservation. Finite cells give exponential thermal suppression for every beta; the rooted all-length space has an explicit low-temperature heat-trace bound from C44's length cost. Combining it with C45's imported C_SW lambda^3 remainder controls pressure but not its T^2 coefficient at fixed lambda.

[proved, obstruction] The raw auxiliary influence expansion can fail the fully summed Grassmann coefficient norm already at lambda^2 and zero temperature. For H=q(Delta+u n_1n_2)+lambda sigma_x, its m-th quartic insertion has absolute integrated coefficient lambda^2 u^m/Delta^(m+1). The sum diverges for u>=Delta. The operator effective Hamiltonian remains the bounded polynomial -lambda^2/Delta+lambda^2 u n_1n_2/[Delta(Delta+u)]. Any proposed polymer proof must resum these operator/contact effects. This is not a no-go for another effective action or for the actual C39/C44 phase.

[proved, finite results] Actual local rules preserve payload U(1), and rooted cells preserve total length. An optional zero-pair test plus a finite density exclusion strengthens C42 to full pair-sum nonresonance of all 2(n+1) branches, conditional on its imported all-n construction. The pair-source time moments have explicit constant B_src=3(1+eta)V0(1+M_h^4). Neither statement is a dressed self-energy estimate.

[proved, thermodynamic reference] The free coefficient is -(pi/12) times the sum of inverse absolute chiral speeds. The -(pi/6) sum over n+1 speeds assumes equal left/right speeds within each pair. PROOF proves canonical free-length transfer, the common-bond quotient partition bounds and equilibrium marker-density transfer for the selected Gaussian physical reference. It does not identify this reference with the finite-coupling microscopic chain. Gaussian current quantization gives c=n+1 if that scaling limit is established; positive-level chiral currents give a conditional unitary-CFT lower bound.

[assumed, remaining work] Prove resummed low-frequency massive matching D*, verify and extend the exact cited RG hypotheses to that action at fixed bare parameters R*, and control the scalar thermal functional or complete stress/current limit T*. Only local literature summaries were available, so the full primary theorem hypotheses were not independently verified. Do not report a microscopic c=n+1 theorem from this folder.

[numerical, checks] checks.py uses numpy, scipy and the standard library. Nine pressure inequalities pass on 2-4 cells with local dimension 3. A thermal Taylor-radius diagnostic and a two-band free-energy quadrature pass. checks.json records the output. No wrapper simulation or central-charge measurement.

[proved, scope] PROGRESS.md records the work. All writes stayed in work_C46, with no network or delegated agents. The project STATE was not modified because the user explicitly forbade external writes.
