# Primary-source infrared audit and a finite-cell reduction

Status: independently reviewed finite-cell and arbitrary-cell massive-matching theorems, a primary-source RG audit, and an obstruction to a naive dedicated-mode extension. The conjecture is not solved here. No new microscopic central charge is asserted.

The strongest new result is [thermal operator-valued matching](INFINITE_CELL_PROGRESS.md). It permits a countably infinite excited cell space: only the finite logical code is expanded into fermion monomials, while the auxiliary operators keep their operator norm. A coarse-time polymer expansion supplies a cutoff-independent analytic radius, spatial/time moments, and a correction proportional to the leakage residual. With an explicit auxiliary heat-trace assumption, it controls the **full thermal influence functional**, including sources and scalar density; thermal kernels converge exponentially to the periodized zero-temperature endpoint. No inverse logical bandwidth enters. [The independent audit](../alternatives/INFINITE_CELL_AUDIT.md) accepts the abstract argument.

The earlier [finite-cell determinant theorem](FINITE_CELL_MATCHING.md) gives a separate route with the simple bound 4r. The operator-valued theorem avoids its dedicated-auxiliary-mode coefficient explosion and is the relevant extension for arbitrary cell lengths.

Its precise inputs are a finite code, a positive auxiliary gap, a finite excited heat trace at some time, exponentially summable bounded coefficient interactions, and a baseline that annihilates the auxiliary vacuum from both sides. Constants depend on the code rank, gap, heat-trace time and support-weight reserve, but not on volume, auxiliary cutoff or the logical bandwidth. C46 supplies a candidate heat bound on its specified legal rooted space. The complete microscopic invariant-space and physical-ensemble application is still conditional; general unbounded length sources are not included in the bounded-source theorem.

Provenance: Codex subagent `infrared_theorem`, 2026-10-28, working in the direct-attack campaign. Derivations below were made in this session. The finite-cell theorem and its arbitrary-cell extension received separate Codex-agent reviews, with corrections incorporated; see [IR_AUDIT.md](../alternatives/IR_AUDIT.md) and [INFINITE_CELL_AUDIT.md](../alternatives/INFINITE_CELL_AUDIT.md). The parent independently checked the published determinant theorem. This is AI review, not human review. Sources were opened in full text. The parent coordinates commits and the shared status.

## What the constructive RG papers actually require

[read] Porta, Scola and Singh, [2509.08665, Assumption 2.1 and (2.18)](https://arxiv.org/html/2509.08665v1), require distinct nondegenerate Fermi crossings and a pair-difference nonresonance condition. They do **not** require distinct absolute velocities; their two-branch example has velocities of opposite sign and equal magnitude. Their microscopic model has finite-range hopping and finite-range density-density interactions. Theorem 2.6 establishes dynamical density/current response. Applying it to arbitrary-degree retarded kernels, or extracting the thermal central-charge coefficient, requires additional work.

[read] Mastropietro and Porta, [2108.13708, Assumption 2.1 and Sections 6, 10](https://arxiv.org/html/2108.13708v3), use the corresponding distinct-root and elastic-scattering conditions for edge modes and construct a comparison with a multichannel Luttinger reference. This is a useful technical framework; it is not a theorem about an arbitrary encoded spin chain coupled to a massive many-body complement.

[proved] The pair-difference condition is precisely uniqueness of unordered pair sums, including repeated entries: rewrite

    p_a-p_b = p_c-p_d  as  p_a+p_d = p_c+p_b.

Equality of the unordered multisets {a,d} and {c,b} gives exactly the two allowed alternatives. Thus C42's repeated-index pair-sum certificate and C46's strengthened length-band search target the correct condition. The extra requirement of all distinct speed magnitudes in the C46 problem statement was unnecessary. This removes a possible distraction, not the infrared matching problem.

## A concrete finite-cell coefficient bound

[proved] Let an operator O act on the Fock space of N complex fermion modes. Expand it in normally ordered monomials and use the coefficient norm

    |O|_h = sum_I |c_I| h^(degree I),      h >= 1.

There is an expansion satisfying

    |O|_h <= (8 h^2)^N ||O||.                         (1)

Proof. There are 4^N occupation-basis matrix units, and every matrix entry has magnitude at most ||O||. A matrix unit can be written as a product of creation operators, the vacuum projector, and annihilation operators. Expanding the vacuum projector produces at most 2^N terms. After the CAR relations are used, each surviving normal monomial has degree at most 2N. Summing the resulting coefficients proves (1). Cancellations can only improve this upper bound. The numerical constant is deliberately loose.

[proved] For a finite parity-graded cell K whose code is F_n, the spectral embedding of [MASSIVE_AUXILIARY.md](../../chat_followup/MASSIVE_AUXILIARY.md) therefore has an explicit coefficient conversion. If m auxiliary modes encode its excited eigenvectors and F=n+m, an interaction supported on r cells obeys

    |J O_X J^dagger|_h <= exp[r F log(8 h^2)] ||O_X||.   (2)

The same bound holds for a difference between the embedded interaction and a chosen logical extension, with the sum of their operator norms on the right. If the original quasi-local interaction is controlled with exponential support weight a > F log(8 h^2), (2) supplies a summable coefficient norm. The support weight must be chosen for the *embedded* number of modes. C45's illustrative a=16(n+1) is not automatically enough for this particular embedding.

[proved] This extra condition does not create a circular scale schedule for a **fixed finite cell**. One can choose a first, compute its finite interaction-norm constants, and choose the final port coupling smaller. The cost can be arbitrarily large with n. It remains finite at fixed n and fixed finite cell space.

## What this says about C46's raw-expansion obstruction

[proved] C46's example has a positive cell Hamiltonian

    A=q(Delta+u n_1 n_2),       u >= Delta > 0.

The invalid raw expansion treats u n_1 n_2 as arbitrarily many external insertions during an auxiliary excursion. It creates coefficients proportional to u^k/Delta^(k+1). A small bound on C45's intercell perturbation J does not by itself control u, since u belongs to A.

[proved] Nevertheless, for a finite cell this is a representation obstruction, not an obstruction to every controlled massive expansion. Diagonalize A first and apply the spectral embedding. All of A is then a quadratic Hamiltonian on dedicated massive auxiliary modes, with energies at least Delta/2. The factor u is part of those energies and is never expanded. Intercell couplings acquire the finite coefficient cost (2). Decreasing their coefficient then makes that norm arbitrarily small without requiring u<Delta. This is the exact finite-dimensional resummation that the raw series omitted.

[read] Pedra and Salmhofer, [Theorem 4.5](https://www.wias-berlin.de/people/koenig/www/FG/PedraSalm.pdf#page=18), control general fermionic interaction kernels with external fields by a pinned-coordinate, field-weighted norm. They require covariance decay and determinant bounds and a small interaction in a stronger field weight. Their theorem is not restricted to quartic vertices.

[proved] The auxiliary covariance in the spectral embedding is diagonal in cell and species. For an auxiliary energy e>=g=Delta/2, its thermal time integral is at most 1/g. More generally, with circular time distance d_beta and 0<=kappa<g,

    integral_0^beta exp(kappa d_beta(t,0)) |C_e(t)| dt
        <= 1/(e-kappa) <= 1/(g-kappa).                 (3)

For 0<t<beta the absolute covariance is e^(-et)/(1+e^(-beta e)); use d_beta(t,0)<=t and extend the integral to infinity. Row and column integrals obey the same bound by time translation and antiperiodicity. These bounds do not depend on the number of cells or beta. Species diagonality means no sum over all auxiliary species appears in a fixed row.

[proved in the accompanying finite-cell theorem] Equations (1)-(3) are combined with the determinant bound, a weighted connected-tree estimate, a normal-ordered equal-time prescription, and explicit scalar/source bookkeeping in [FINITE_CELL_MATCHING.md](FINITE_CELL_MATCHING.md). A Cauchy estimate gives a correction proportional to the residual r, rather than the square of the larger auxiliary coupling. This applies to finite cells with an exactly auxiliary-vacuum-reducing baseline. The finite-temperature baseline and the infinite-cell extension are explicitly outside that conclusion.

## A precise obstruction to taking the all-length limit this way

[proved] A naive cutoff-uniform extension of this *coefficient norm* fails. The spectral embedding uses one or two dedicated auxiliary modes per excited eigenvector. With m modes, the projector onto the auxiliary vacuum is

    P_vac = product_(s=1)^m (1-a_s^dagger a_s).

Its unique normal-monomial expansion has

    |P_vac|_h = (1+h^2)^m.                              (4)

Every subset of the m modes contributes one different monomial, with coefficient +1 or -1; summing proves (4). This diverges as m tends to infinity for every fixed h>0. In particular an operator-norm-one vacuum projector is not uniformly bounded in this field-weighted coefficient norm.

[proved application boundary] C45's ordered-cell Hilbert space includes arbitrary lengths. A length cutoff makes the cell finite, but changes the Hamiltonian: marker motion can carry a cell across the cutoff. Applying (2) after each cutoff gives no uniform estimate, and (4) shows why operator norm and a growing excitation energy alone do not repair that route. The new operator-valued theorem supplies a different resummation: it never expands the auxiliary vacuum projector. Thus this obstruction still rules out the naive dedicated-mode limit, but no longer rules out the available abstract all-length matching route.

This is an obstruction to an often tempting proof shortcut. It neither disproves massive matching for C44/C45 nor produces a competitor state of lower energy.

## Which route currently looks shortest

[assessment] The strongest next mathematical target is now the retarded-action multichannel RG and physical ensemble transfer. The arbitrary-cell matching theorem removes a concrete analytical obstruction for the mobile-wrapper route, provided the complete invariant-space and heat hypotheses are instantiated. The frozen finite-cell version remains a useful control with fewer ensemble complications; its physical scope is treated separately by the parent.

[assessment] If the goal is a physics-level counterexample at the smallest stated local dimension, the d=24 vacancy candidate remains more economical, but its spatially weighted quantum-selector and dressed-defect comparison is still a substantive extra problem. C44/C45 is useful as a structurally different backup: its selected cells have a genuine local auxiliary threshold, at the price of an all-length Hilbert space and extraordinarily severe scales.

[assumed, unresolved] Neither route has been proved here to realize the claimed Gaussian infrared limit at fixed nonzero couplings. No present result licenses reporting the c-d conjecture as refuted.
