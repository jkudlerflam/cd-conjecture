# C40. Infrared transfer for the hierarchical chain

[heuristic, verdict] The hierarchy of velocities does not obstruct the fixed-n infrared limit. Given a nonresonant quadratic reference, controlled local massive-sector matching, and sufficiently small interactions measured against the slowest band, the conventional RG argument gives n compact bosons for frozen cells and n+1 for mobile cells. Every coupling can be fixed and nonzero before taking the thermodynamic limit.

[proved, limitation] This is not an unconditional all-n result for the supplied bonds. The missing all-n finite-dimensional statement is a simple-root, nonresonant reference for the actual spatial core. The remaining dynamical input is a bound on the exact irreducible low-frequency action after massive-sector elimination. Finite-order connected locality is proved below, with an explicit coefficient bound. It does not prove that last input. The requested allowance for a precisely located obstruction is therefore used.

[proved, file map] [PROOFS.md](PROOFS.md) contains the derivations, assumptions and order of choices. [checks.py](checks.py) and [checks.json](checks.json) contain small exact combinatorial and rational checks. [PROGRESS.md](PROGRESS.md) records the work and neighboring results. No microscopic central charge was measured.

## 1. What the argument assumes

[assumed, F] The full selected-cell code and positive excluded-sector bounds hold, including comparison with wrong-length, illegal and marker-free sectors. The complete physical port and marker norms, including compensators and visibility factors, are included. C32 supplies the positive mobile core and its inherited wrapper bounds. Some unrestricted-sector bounds for the original frozen 23-letter bond remain unfinished in its README.

[assumed, B] At the chosen fixed n, the leading normalized quadratic symbol has 2n simple Fermi roots, nonzero minimum speed, and positive number-conserving quadratic and quartic momentum margins. The mobile length density has the additional mixed margins. A finite-g quadratic certificate must also leave room for the non-Gaussian interactions, not just for further quadratic perturbations.

[assumed, M2] Integrating out the massive cell variables, with retained propagation and operator mixing subtracted locally, gives a symmetry-preserving, spatially summable low-frequency action. Its complete second-order coefficient is the calculated bond, and its remaining local vertices are o(lambda^2) at fixed core parameters. Frequency and filling derivatives are controlled in the same norm. PROOFS section 3 gives a concrete sufficient kernel norm and quantitative version. The excluded sectors remain energetically disfavored.

[assumed, R] Conventional weak-coupling lattice matching and compact-boson RG apply continuously in that norm, including the collective irrelevant-vertex norm, the dressed fillings and the full density-current stiffness matrix. This is the theoretical-physics input, not a rigorous lattice-to-CFT theorem.

[assumed, convention] The mobile local operator lattice and cyclic gluing are those of C32/C8's identity payload transport. A collection of compact bosons with different velocities is assigned the sum of its mode central charges. A common Lorentz velocity and exact ground-state uniqueness on every finite ring are not inferred.

## 2. Massive-sector elimination, what is controlled

[proved, inherited] In ordered cell coordinates the mobile model is a sum of independent on-cell Hamiltonians plus bounded nearest-cell perturbations. The on-cell space includes all lengths. The selected projectors are fixed, even though markers move in physical coordinates. The C8/C36 reduction therefore applies to C32's positive core once its stated wrapper hypotheses hold.

[proved, new finite-order bound] Let A_j>=Delta(1-P_j), let J bound the local perturbation norm before multiplying by its expansion parameter, and use an exponentially weighted interaction norm. The connected order-r coefficient has support in at most r+1 cells and obeys the conservative bound

    ||Phi_r||_mu <= C_r(mu) J^r/Delta^(r-1),
    C_r(mu)=(r+1)^2 2^(r+1) exp(mu(r+1)) [16(r+1)]^r.

On-site terms can be assigned to bonds or absorbed by enlarging the constant. The proof uses a Riesz projection on an interval of at most r+1 cells at radius Delta/[16(r+1)J], its analytic canonical direct rotation, and inclusion-exclusion over bond subsets. Disconnected coefficients vanish because the canonical code identification factorizes. No extensive interacting band is assumed isolated. See PROOFS section 2.

[proved] Each irreducible unperturbed fast interval has energy at least Delta. Reduced-resolvent integral formulas give finite temporal moments at every fixed order. A word with r-1 such intervals acquires at most the factor

    (r-1)(r)...(r+m-2)/Delta^m

for its m-th total-duration moment. Retained intervals are kept in the slow theory. The canonical fourth-order normalization term is included, as in C26. Symmetries survive every step.

[proved] A connected r-insertion term has physical support at most (r+1)(b+1)+2r sites. Moving markers do not destroy fixed-order locality. Conversion from operator norms to local fermion coefficients costs finite, n-dependent factors, stated in PROOFS.

[proved, remaining obstruction] These constants do not give a summable all-order series. The restricted active-defect estimate in C36 does not include histories that create another defect pair away from existing defects and join it later, or histories that propagate a changed retained state after a local defect disappears. Their cancellation against retained propagation is the missing connected-action estimate. A global bare-Q inverse, an ordinary connected bare-source correlator, or fixed-slow-time dynamics cannot replace it. C8's explicit controls disprove those substitutions.

[assumed, sufficient replacement] With q the full local perturbation strength, an adequate exact low-frequency bound is

    ||R_>=3||_(LE,mu) <= A_n q^3/Delta_n^2,
    ||K2(omega)-K2(0)||_(LE,mu)
                         <= A_n q^2 |omega|/Delta_n^2.

It need not come from a convergent perturbation series. At fixed g and wrapper parameters, q=O(lambda), while the payload scale W is proportional to lambda^2. Hence the first error divided by W vanishes with lambda; at omega=O(W), the second does also. This is the correct second-order transfer criterion. Counting the whole second-order term as an error would leave an order-one mismatch.

## 3. The RG through separated bandwidths

[proved, normalization] Write the leading logical Hamiltonian as W Hbar, with W=lambda^2 g^2 for the frozen core, or W=lambda^2 g^2/s for the mobile core. All attenuation and small Gram eigenvalues remain in Hbar. Let bar_v_* be its smallest speed, including the length mode when present.

[heuristic] First lower the cutoff through the massive cell gaps using M2. Keep every retained band with its full lattice dispersion. The fast bands enter their linear regimes first. A slower band remains in the action and enters its own linear regime at its much smaller curvature scale. A small bandwidth is not a mass gap, and integrating a fast flavor out globally is also invalid because it has its own Fermi points.

[proved, finite reference data] Simple roots admit disjoint monotone patches with a positive isolation gap and a positive energy minimum outside the patches. These give a common E_lin>0 below which every band has linear patches. If delta is the smallest reserved payload or mixed mismatch, a conservative final matching scale is

    E_det=min(E_lin, W bar_v_* delta/32).

Below E_det, the momenta within four retained patches cannot compensate delta. All dangerous oscillatory quartics are then absent from the low-energy momentum support.

[heuristic, hierarchy cost] Matching from the largest bandwidth E_top to E_lin has finite constants at fixed n. These can contain powers of bandwidth ratios, not just logarithms. Include them in A_match,n. The subsequent marginal running can be bounded conservatively over the entire interval

    ell_* = log(E_top/E_det).

A criterion using only log(1/delta) can omit a very large bandwidth-ratio contribution. No uniform-in-n estimate of A_match,n is asserted.

[proved conditional on the RG majorant] Let x0 be a collective dimensionless norm of all marginal deviations, including forward couplings, and Y0 the irrelevant norm after matching. For a positive irrelevance margin sigma, set

    X0=x0+C(Y0+Y0^2)/sigma.

The standard weak-coupling majorant gives

    x(ell)<=X0/[1-A X0 ell].

Thus A X0 ell_*<1/2, together with a reserved Gaussian and root-shift budget, prevents strong coupling before every slower band has formed and every dangerous mismatch has cut off. Including forward terms in x is necessary because their stiffness changes alter the dimensions of oscillatory quartics. Flavor sums belong in the collective norm; an entrywise bound alone is insufficient.

[heuristic] The already-linear modes stay gapless throughout this interval because their roots cannot annihilate within the C1 perturbation margin and their non-Gaussian couplings cannot reach strong coupling. Below E_det only Gaussian marginal terms and irrelevant interactions remain. This supplies the requested stability mechanism during the successive band crossovers.

## 4. Generated operators and their fate

[proved] Exact total U(1) forbids anomalous pairing at every order, including third order. It does not forbid interflavor pair transfer with two creations and two annihilations. That process requires a quartic momentum margin or a separately established irrelevant basin.

[proved, source-specific orders] C31's central switch has Gaussian terms at lambda^2 tau^2 and lambda^3 tau^2, and non-Gaussian terms beginning at lambda^2 tau^4 and lambda^3 tau^4. Two conversions cannot produce a normal-ordered quartic. Its cubic Gaussian term changes velocities and Fermi points. C31's favorable centered pair-transfer flow at n=2,3,4 is numerical and does not supply an all-n result.

[proved, source-specific distinction] For the actual 23-letter folded bus and C32's core, alternate-cell excitation parity removes odd port orders in connected bulk coefficients. The leading non-Gaussian second-order remainder is O(lambda^2 g^4). The mixed u lambda term vanishes; pure marker excursions can generate length interactions. A rho^2 suppression from C26's reflection family is not assumed for either spatial core.

[proved conditional on the bosonization dictionary] For chiral charges z_i with sum z_i=0,

    D0(z)=sum_i z_i^2/2,
    Q(z)=sum_i z_i p_i mod 2pi.

The only nonzero primaries with D0<=2 are bilinear differences and four-distinct-branch two-in/two-out vectors. They are excluded by delta_2 and delta_4. Equal incoming and outgoing multisets yield densities and current products. Apparent exchange with that same multiset is a density interaction after anticommutation. Different outgoing pairs are oscillatory and stop at E_det.

[proved conditional] Mobile length harmonics add h^2 to D0 and -2pi h nu to Q. Marker conservation forbids length dual-field vertices. The only additional low-dimension exclusions are h=+/-1 alone and h=+/-1 times a payload bilinear. An interior nu avoiding finitely many values gives a positive mixed margin. Irrationality alone is neither necessary nor sufficient for this finite test.

[proved within Gaussian theory] Let H0 be the positive free stiffness, including all speed ratios, and H the dressed stiffness. The relative norm

    eta=||H0^(-1/2)(H-H0)H0^(-1/2)||

controls the scaling metric, with lower bound sqrt((1-eta)/(1+eta)) times its free value. If eta<=1/4, every remaining nonzero primary of free dimension at least three has dimension at least 3 sqrt(3/5)>2. This bound includes density-current and length-payload mixing. It does not require comparable velocities.

[heuristic] Oscillatory quartics give finite pre-cutoff renormalizations, forward/exchange density terms change Gaussian couplings, and the remaining higher vertices decay. Their complete feedback must fit the reserved eta and filling budgets. Under those bounds, the generated non-Gaussian terms do not remove a mode. This conclusion depends on the margins and smallness, not on U(1) alone.

## 5. A consistent order of choices

[proved conditional on the stated inputs] First fix n and a certified leading nonresonant band shape, including its loop and normalized field. Then reduce g until the lambda^2 g^4 remainder divided by W meets the complete matching/RG tolerance. Keep the leading shape fixed. Tying the loop to g^2 without a separate margin estimate need not improve the error-to-detuning ratio.

[proved conditional] For mobile cells next choose w/U and then s/delta_fast to make all four second-order length-pair bonds agree within the same tolerance. Include their scalar parts. Dividing C32's scaled equalization error by g^2 gives its contribution in W units. Recompute the actual positive excluded gap after these choices. Choose the length density and set u beta_n=W j_l.

[proved conditional] Finally take one sufficiently small lambda>0. The higher-port errors divided by W vanish even under a crude cubic bound with an inverse g^2 loss. The sharper lambda^4, lambda^2 u and u^2 bounds also vanish after u=O(lambda^2). A fixed exact code has no residual splitting that grows relative to W as lambda decreases. Retune the marker chemical potential for the length-independent scalar energy. Keep the remaining scalar length interaction in the budget.

[proved, order of limits] Hold every parameter fixed, take the thermodynamic limit, then examine energies below E_det. Only after defining that limit for each n may the family be compared as n increases. No coupling is sent to zero with ring length.

## 6. Current finite witnesses and acceptance

[proved, source attribution] C37 has now completed certificates for the original frozen 23-letter family at every n=1,...,23, with simple roots, all pair-sum margins and a finite-conversion schedule. Its largest carrier inverse has dimension 2209 and pair-product degree 1081. C33 independently provides a sharper n=23 finite-g certificate with delta_2>0.0140193 and delta_4>3.1558e-7. These are inherited finite results, not recomputed in C40. They do not certify C32's different positive mobile core or the thermodynamic action.

[proved, coordination] C38 establishes that bandwidth separation alone does not imply the required nonresonance and gives a conditional endpoint target. C39 is investigating a modified fixed-alphabet core; its last progress read does not claim a completed all-n construction. C40 does not duplicate their band or compiler work.

[heuristic, answer to the acceptance question] Yes for the conditional physics inference at a fixed n satisfying F, B, M2 and R. The hierarchy can be absorbed into finite matching constants and stricter choices of nonzero couplings. It does not force an additional gap or invalidate bosonization. No for an unconditional statement that the supplied full frozen and mobile families have c=n and c=n+1 for every n. The all-n band premise and the exact connected low-frequency remainder remain unestablished, and the frozen global-sector audit remains separate. A careful theoretical physicist can use M2 as an explicitly stated Wilsonian assumption, but should not describe it as proved by the available finite-order calculations.

[heuristic, infrared interpretation] Under those inputs the oscillator count is n or n+1, with c equal to that count. Frozen cells retain broken-translation sectors. Mobile gluing changes zero modes and boundary twists. Unequal velocities give a velocity-weighted Casimir energy, so fitting that coefficient without the velocities does not measure the unweighted c.

[proved, checks and limits] Running `python3 -B checks.py` passes exact charge-vector enumeration at n=1 through 6, a rational slow-mode stiffness example, the Gaussian threshold arithmetic and an RG-majorant consistency check. No floating-point numerical result is used for a new finite-dimensional sign claim. No network, external file modification, subagent, cluster computation, microscopic many-body simulation or independent referee was used. The project STATE is outside the authorized write scope; local STATE.md supplies the handoff.
