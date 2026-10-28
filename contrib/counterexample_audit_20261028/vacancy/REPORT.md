# Independent vacancy, locality, and CFT-scope audit

**Status: exact kinematic identities re-derived; finite checks passed; conditional infrared conclusion.** No previously recorded endorsement is used as evidence. I found no hidden low-dimension neutral vertex or charge-momentum cancellation that refutes the proposed vacancy completion. The resulting statement is conditional on local finite-coupling matching, preservation of the material-defect gap, and a quantitative Gaussian neighborhood. A positive Gaussian matrix by itself is insufficient. These qualifications distinguish a missing hypothesis from a refutation of the model.

Audit by a GPT-6-family Codex subagent on 2026-10-28; no human review. Assigned scope: the vacancy map, finite-ring seam, charge stiffness, microscopic locality and vertex lattice, arbitrary-ring gluing, and velocity/CFT scope. The internal finite-dimensional certificate and its inertia proof are inputs here, audited by other workers. Read `CONTRIBUTING.md`, the complete d=24 synthesis, its charge/infrared appendices, and the literal 23-letter bond and code description.

## 1. Exact ring map, including what the gap variables omit

Fix occupied number N, with 0<N<L. Ordered physical positions and labels give the exact vector-space identification

    H_(L,N) = C^(choose(L,N)) tensor (C^23)^tensor N.

A vacancy exchange away from the cut preserves the ordered label word. A particle crossing L−1 to 0 rotates the word right by R. Therefore, with the −t amplitude included in B,

    H_t = T_bulk tensor I + B tensor R + B† tensor R†.

The corresponding fermion cut coefficient is (−1)^(N−1)R. This parity factor is necessary even though the microscopic particles are hard-core bosonic labels. An occupied interaction is exactly the oriented material bond h_(j,j+1) multiplied by the adjacency indicator A_j. In particular the seam bond has its original orientation; no reflection symmetry is assumed.

The gap representation needs more than a list of N integers. A complete rooted coordinate is

    (x_root; g_1,...,g_N; sigma_1,...,sigma_N),
    sum_j g_j=L−N,

with the root equivalence

    (x;g;sigma) ~ (x+1+g_1; rotate_left(g); rotate_left(sigma)).

It has N equivalent roots. Its dimension is consequently

    L choose(L−1,N−1) 23^N / N = choose(L,N) 23^N.

Dropping both x_root and the quotient would give the wrong full-ring Hilbert space. At fixed physical momentum the root coordinate becomes phases, together with the cyclic consistency condition. At momentum zero and with a distinguished material tag the simple gap Hamiltonian is literally obtained. For the actual untagged chain, the same bulk Hamiltonian has the required seam and quotient. This is a precision correction to any unqualified assertion of an exact tensor product of gaps and labels; it is not a bulk obstruction.

Moving material particle j by one physical site transfers one vacancy between its two neighboring gaps, with amplitude −t independent of occupancy. Thus the bulk expressions are exactly

    H_gap = −t sum_j (B_j† B_(j+1)+h.c.),
    B|g> = |g−1> for g>0, B|0>=0,
    V = J sum_j 1[g_j=0] h_(j,j+1).

No square-root factors are present. The auxiliary infinite alphabet is a coordinate description, not the physical site dimension.

Within an R=e^(i theta) sector the charge orbitals are

    k_n = [2 pi n−theta−pi(N−1)]/L.

For a consecutive sea, Wick contraction gives

    alpha_(L,N) = N/L − sin²(pi N/L)/(NL sin²(pi/L)),
    alpha(rho) = rho − sin²(pi rho)/(pi² rho).

Projection onto all minimizing charge seas gives alpha_(L,N) times the original material Hamiltonian in each cyclic sector. Different retained seas have different physical momenta; uniform V cannot mix them. This includes the two half-flux seas. Individual rooted adjacency probabilities need not equal alpha. The projection is a first-order identity, not a uniform elimination of the gapless charge fluid.

**Independent check:** `check_independent.py` constructs the physical matrix directly from local exchanges and a generic complex oriented two-site bond, then separately constructs the squeezed matrix. Five cases, L=4,...,7, dimensions 24–280, give zero kinetic and interaction entry error. The full charge-sea projection error is at most 7.33e−15, including half-flux degeneracy. No campaign implementation is imported.

## 2. Charge stiffness independently recovered

Write r=(L−N)/N and rho=1/(1+r). The free-fermion physical energy per physical site is −(2t/pi)sin(pi rho). Dividing by rho converts to energy per material gap:

    e_gap(r) = −(2t/pi)(1+r) sin[pi/(1+r)],
    e_gap''(r) = 2 pi t rho³ sin(pi rho).

Material position j≈rho x, so the velocity becomes

    v_gap = rho v_physical = 2t rho sin(pi rho).

Using the compressibility normalization e''=pi v/K gives

    K_gap=rho^(−2),
    rho=1/2: v_gap=t, K_gap=4,
    alpha=1/2−2/pi²=0.2973576327... .

Equivalently the displacement field satisfies phi_gap=−phi_particle/rho to linear order; the physical free-particle K=1 therefore becomes K_gap=1/rho². There is no conflict between these two stiffnesses. The local neutral phase-model harmonics e^(2im phi_gap) have dimensions 4m² at half density.

The primary q-boson paper's q=0 operator amplitudes, velocity and dressed charge agree: multiply its hopping by 2t and translate its K convention by pi. See [Bortz–Sergeev, §§2.1, 4.1, equations 42–46](https://arxiv.org/html/cond-mat/0603093). This source was opened directly; the argument above does not assume a value copied from the campaign.

**Independent stiffness check:** diagonalizing the integer-gap hopping matrix, and measuring both particle-number curvature and adiabatic flux curvature, gives K=4.03087, 4.01713, 4.01089, 4.00754, 4.00553 for 3–7 material sites at r=1. The corresponding velocities approach t from 1.03918t to 1.00705t. Ground energies agree with −2t/sin(pi/(2N)) within 8e−15. This tests the nontrivial factor of four, not merely the thermodynamic derivative formula.

## 3. Microscopic local algebra and the fractional harmonics

The code charge is the sum of its binary data occupations. A finite-support, charge-conserving code matrix unit has equal numbers of raising and lowering operators. Under the Jordan–Wigner map, all strings outside that support cancel. What remains is a finite even fermion polynomial, possibly with finite products of number/parity factors inside the support. On the ring, the usual Jordan–Wigner boundary condition depends on the total internal charge parity; a neutral Hamiltonian term stays in that parity sector. Consequently internal half-integer disorder/twist vertices are not extra local Hamiltonian terms merely because the underlying variables are spins. Odd fermion operators or a single spin raising field carry the exact internal U(1) charge and are forbidden as perturbations.

This reasoning requires the low-energy matching to be quasi-local and symmetry preserving. It is not a proof that an arbitrary nonlocal change of variables preserves locality. Here a finite physical interval contains at most as many occupied particles as physical sites, so its neutral operator acts on a bounded material interval. Projection onto fixed finite cells preserves that property. A microscopic neutral quasi-local interaction can therefore be treated through the finite even-polynomial algebra. Semi-infinite Jordan–Wigner disorder strings would require additional justification and are not generated by that local algebra.

Marker harmonics are different: they really are local physical observables with fractional charge-field exponents. At zero-defect period b=24, write

    N(x)=rho x−phi_p(x)/pi, y(x)=N(x)/b.

A payload field with cell momentum q, including its cell comb, has the branches Q=q+2pi ell, ell integer. Substitution gives

    e^(i Q y(x)) = e^(i rho Q x/b) e^(−i Q phi_p/(pi b)),
    k_phys = rho Q/b mod 2pi,
    Delta_ref = Delta_pay + [Q/(2pi b)]².

The marker harmonic itself is q=0, ell=1: its dimension is 1/b² and its correlation exponent is 2/b². It cannot be discarded by occupied-number conservation. Its nonzero physical momentum nevertheless forbids it as a uniform perturbation. Exact occupied-number conservation excludes a net dual-charge vertex; exact internal U(1) excludes a net payload charge. The independent marker symmetry and the positive mismatch-defect gap are needed to restrict the theory to the selected locked material sector.

Uniformity at rho=1/2 requires Q=4pi b h, h integer. Therefore

    Delta_ref,res = Delta_pay + 4h².

For h=0 the selection rule is exactly q=0 modulo 2pi. For h nonzero the charge contribution alone is at least four. This includes fractional marker branches and ordinary charge umklapp in one lattice; adding a separate unconstrained charge exponent would double-count degrees of freedom.

At the free payload reference a chiral charge vector l has dimension (1/2)sum_a l_a². Neutral nontrivial operators below dimension three are quadratic or quartic fermion operators, including their lower-dimensional contractions and descendants. Equal incoming/outgoing branch multisets give densities and current products, hence Gaussian terms. Repeated same-chirality fields require derivatives or have dimension at least three; checking pairs with replacement is a conservative enlargement. All remaining neutral polynomials have dimension at least three. There is no unexamined arbitrarily high-order operator below this cutoff.

**Independent enumeration:** `vertices.py` reads only the saved rational root intervals. It redoes exact circular root and pair-product comparisons, then explicitly evaluates every ordered transfer and every quartic marker branch with reference dimension below three. It covers 46 roots, 2,070 ordered nonforward quadratic transfers, 1,081 unordered incoming pairs with replacement, 1,167,480 ordered nonforward quartic transfers, and 56,039,040 dangerous quartic harmonic branches. No resonance occurs. Its rational lower bounds are

    delta_2 = 0.014019386452133372...,
    delta_4 = 3.155876454096123e−7...,
    delta_phys >= min(delta_2,delta_4)/48
               = 6.574742612700256e−9... .

The root intervals' inertia certificate is inherited, not re-proved by this script. The exact pair-product comparison supplies the rigorous bound; the all-transfer floating-point enumeration is a separate completeness check.

## 4. Boundary twists do not provide a cheap momentum-cancellation channel

At half filling an arbitrary charge excitation has equal particle and hole counts. For d(k), the distance to {−pi/2,+pi/2},

    |cos k| = sin d(k) >= 2d(k)/pi.

Replace every particle and hole momentum by its nearest Fermi edge. The signed edge sum is an integer multiple of pi. The triangle inequality therefore yields, for any number of particle-hole excitations,

    E_exc >= (4t/pi) dist(P,pi Z).

A material-twist change Delta theta shifts the charge-sea physical momentum by Delta theta/2 modulo pi. To conserve physical momentum requires

    E_target−E_source >= (2t/pi)dist(Delta theta,2pi Z)
                         −2t tan[pi/(2L)].

Minimizing over the b branches Delta theta=(q+2pi ell)/b leaves the positive limit 2t dist(q,2pi Z)/(pi b). Thus the O(1/L) splitting between unconstrained seas is not an available nonforward process at conserved physical momentum. Soft charge forward modes must still be retained. I re-derived the bounded-rounding extension: if N−L/2 stays bounded, comparison with the set of all negative-energy orbitals introduces only O(1/L) corrections; its signed Fermi-edge contributions still cancel modulo pi.

**Independent exhaustive test:** all occupations and all internal twists at L=6,8,10,12,14 satisfy the excitation inequality, with worst negative numerical slack 2.35e−15. The L=14 test includes 24,024 state/twist combinations. The inequality itself, not this finite test, controls arbitrarily many excitations.

## 5. Precise conditional infrared statement

For this fixed n=23, the following is sufficient within conventional local weak-coupling bosonization/RG:

1. The selected code has a positive dressed gap to auxiliary excitations and wrong-length/illegal material defects. Its neutral local algebra is the one described above.
2. Local matching produces a charge boson plus 23 payload right/left pairs, with all velocities nonzero. Their reference root and low-order momentum margins persist at the matching scale.
3. In canonical chiral coordinates, let Sigma=diag(I_24,−I_24), let V_0 be the positive diagonal reference velocity matrix, and let V be the matched Gaussian Hamiltonian. Require (1−eta)V_0<=V<=(1+eta)V_0 with eta<5/13, rather than positivity alone.
4. The remaining quasi-local non-Gaussian coefficients are sufficiently small for the finite matching interval above the oscillatory cutoffs and for perturbative RG below it. The matching preserves the exact charges and the vertex lattice. No new massless defect sector appears.

Then every resonant non-Gaussian perturbation has dimension greater than two; every other potentially relevant or marginal non-Gaussian term has a fixed nonzero momentum cutoff. The surviving infrared theory has 24 right and 24 left Gaussian oscillator families, compact zero-mode gluing, and Virasoro c_L=c_R=24, conventionally reported as c=24, not 48.

For completeness the dimension metric is

    G = V^(−1/2) |V^(1/2) Sigma V^(1/2)| V^(−1/2)
      = (V^(−1)) # (Sigma V Sigma),

where # is the positive operator geometric mean. Its monotonicity bounds G below by sqrt[(1−eta)/(1+eta)] times the identity. Hence every reference dimension at least three stays above two if eta<5/13. This also covers asymmetric right/left velocities and mixed density-current couplings. See the coordinating audit's `../GAUSSIAN_BOUND.md` for its independent derivation. No equality of cone slopes is needed, and no comparison with the fastest velocity alone suffices.

The physical velocity corresponding to a cell-coordinate velocity u is (b/rho)u; the gap-chain velocity corresponds to v_physical=v_gap/rho. These finite conversion factors must accompany any physical finite-size fit. A local small correction compared with t can still be huge compared with the slowest payload speed; the relative bound above is substantive.

The statement does not assert that the four hypotheses follow from the finite-band certificate. In particular the exact J=0 projection and free gap stiffness do not prove hypothesis 1 or the operator-valued, finite-J part of hypotheses 2–4. Locality of the gap representation makes the required matching a conventional coupled-chain problem, but does not itself derive its coefficients or bound its defect energies. This is the remaining conditional step in this audit, not an actual refutation.

## 6. Arbitrary rings, compactness, and the original conjecture

Assuming the dressed defect gap, N can be rounded to a multiple of b with |N−L/2|<=b/2. Positive charge compressibility makes this cost O(tb²/L), whereas a wrong material cell has a fixed positive cost O(J). This works on prime and odd physical rings and does not require a fixed vacuum zero-mode sector for every L.

The b material-origin branches satisfy K=(P+2pi r)/b and become charge fluxes. They are not b independent frozen physical crystals. Jordan–Wigner parity and Fermi-sea shell effects add further zero-mode conditions. They restrict or shear the joint charge/payload lattice; they do not project out a local oscillator. Fermi momenta may give irrational coefficients in that shear. This is not decompactification: the full joint charge lattice remains discrete and full rank even when its projection onto one coordinate is dense. Its fractional origin sectors supply the complementary charge windings. The marker harmonic has positive exponent 2/b², so ordinary physical positional long-range order is absent at the reference.

The original conjecture asks about a periodic one-site-invariant Hermitian nearest-neighbor chain and the sum of the central charges of its unitary low-energy CFTs. It expressly permits different propagation velocities. It imposes no all-L uniqueness requirement. See [Latorre–Sierra, §.2, equation (3) and the following velocity discussion](https://arxiv.org/html/2403.17242v2), opened directly. Unequal tilted chiral cones are handled by the positive Gaussian form above; they change the Hamiltonian's velocity coefficients, not the oscillator Virasoro anomaly. The finite-ring gluing must still be specified when comparing exact spectra, but these features do not exclude a 24-mode Gaussian phase from the conjecture's stated scope.

## Reproduction and disposition

Run from the repository root, using Python with NumPy and SciPy:

    python3 contrib/counterexample_audit_20261028/vacancy/check_independent.py
    python3 contrib/counterexample_audit_20261028/vacancy/vertices.py

Outputs are `independent_results.json` and `vertices_results.json` in this folder. Recorded runtime is below one second per script in the provided runtime. The scripts import no campaign test implementation and write only this audit folder.

**Disposition:** no vacancy-specific refutation found; exact ring structure and charge normalization verified independently; low-order vertex exclusion verified conditional on the supplied certified roots. Correct the full gap-coordinate statement by retaining its root quotient. Treat c=24 as following from the explicit conditional RG statement above, not as an independently measured microscopic central charge or a theorem proved by the finite-size checks.
