# C117 result. Bootstrap bounds on A

2026-10-29.

[Open] None of the three routes yields a general crossover-free inequality A <= C M f(d). No universal C or f is established. The work produces a fixed-support obstruction, rigorously certified small SDP bounds, and a microscopic example showing that even the q^4 energy spectral moment receives massive contributions. The unresolved inputs can be stated as explicit inequalities.

## Ground-state positivity and spectral moments

[Proved] For an optimized bond with norm M, write

    m_1(q) = <[h_-q,[H,h_q]]>/(2L)
           = b_1(cos q-1)+b_2(cos(2q)-1).

Locality and stationarity give |b_1|<=8M^3 and |b_2|<=4M^3. Therefore

    |[q^4] m_1(q)| <= 3 M^3.

This is an actual lattice inequality with constant 3, but it does not bound A. Long-distance nested energy commutators vanish identically once their supports separate.

[Proved] The optimized transverse-field Ising bond

    h_i = -J X_i X_(i+1) - (g/2)(Z_i+Z_(i+1)), 0<J<g,
    M = sqrt(g^2+J^2),

has A=0 but

    [q^4] m_1(q) = gJ^2/4 + O(J^4/g) > 0

for small J/g. Thus the coefficient at the same q^4 order as the CFT contribution cannot be identified with sum c_a v_a^3/(24pi). The calculation uses a gapped model in an optimized gauge, not a density coboundary. Even a bound on sum c_a v_a^3 would not by itself bound sum c_a v_a when arbitrarily many slow sectors are allowed. Derivations and a small-ring check are in [PROOF.md](PROOF.md), Section 6.

[Assumed infrared matching; proved conditional implication] A shrinking positive spectral window does exclude fixed-gap contributions. Define

    F(q) = integral_(0,Lambda|q|] Sbar_h(q,omega) d omega/omega,
    Lambda > every chiral velocity.

Stress-tensor matching gives liminf F(q)/q^2 >= A/(24pi). The missing microscopic lemma is

    limsup_(q->0) F(q)/q^2 <= K M f(d).

It would prove A <= 24pi K M f(d), without a fixed matching scale. Nonnegative-moment upper bounds and finite local norm bounds do not supply this inverse-moment bound. Operator matching itself is an additional assumption beyond the ring energy expansion.

## Ring SDPs and certificates

[Proved] Every translation-invariant local program whose complete constraint supports fit in R sites becomes independent of circumference for L>2R. A shorter n-site ring ground state is feasible whenever n>2R, including all such local ground-state constraints. Thus its optimum obeys

    ell_R <= E_0(n)/n.

For A>0, choosing a sufficiently large fixed compatible n gives e_infinity-ell_R>0. The putative upper bound on A from E_0(L)>=L ell_R consequently grows as L^2. With compatible shorter lengths n=O(R), even R=o(L) cannot give energy error o(1/L). This obstruction does not cover global Gram matrices, winding operators, or explicitly circumference-dependent constraints.

[Proved] Fifteen rational local lower-bound certificates and ten rational full-ring energy intervals passed independent integer checks. The local experiments use k=2,3 on Heisenberg, XXZ at Delta=1/2, SU(3) swap, and seeded real random bonds at d=2,3. The full-ring reference SDPs use L=4,6 for d=2 and L=3,4 for d=3. Exact three-site optima, unchanged by all one-site ground-state positivity constraints, are

| Bond | Exact three-site lower energy per site |
| --- | ---: |
| Heisenberg S.S | -1/2 |
| XXZ, Delta=1/2 | -(1+sqrt(33))/16 |
| SU(3) swap | -1 |

[Proved] Heisenberg's dual is the local SOS

    H_L+(L/2)I = sum_i [(h_i+h_(i+1))/2+I/2],

with each three-site summand positive. Its energy density is fixed by the four-site ring and does not track the large-ring Casimir term. Full-ring certificates can resolve the finite energy to widths below 3.83e-6, but this is the full k=L problem and supplies no uniform asymptotic theorem. General ground-localizing certificates are positive in ground-state expectation and are not automatically operator SOS certificates for C114.

[Numerical limitation] Random d=3 k=3 searches hit their iteration cap. The verified rational lower bounds remain valid, but their optimality is unproved. No significant improvement from the added ground constraints is established in this small sample. [CHECKS.md](CHECKS.md) records endpoints, residuals, solver limitations, and reproduction commands.

[Assumed ring-certificate lemma; proved conditional implication] A concrete sufficient target is a certified lower ell_L<=E_0(L) and upper u_2L>=E_0(2L) satisfying

    limsup L(u_2L-2ell_L) <= K M f(d).

The Casimir expansion then gives A <= (4K/pi) M f(d). This is the unresolved certificate-rate lemma. Two lower energy bounds cannot be subtracted to obtain it.

## Modular and conformal bootstrap

[Proved limitation of the stated inputs] A k-site algebra bounds the rank of a specified scaling map by d^(2k). It does not establish a surjective map onto the fast-sector fields. One stress tensor counts one factor, independently of its central charge. Modular invariance contains no lattice support size or M. The standard modular gap bound is in the wrong direction to convert a gap lower bound into a central-charge upper bound. See the primary references and qualifications in [PROOF.md](PROOF.md), Section 8.

[Assumed operator lemma; proved conditional implication] Let V=4eM and N(u)=sum_(v_a>=u)c_a. If a field space of dimension at least eta N(u) is the image of operators on

    k(u) <= k_0(d)+beta(d) log(V/u),
    gamma(d)=2 beta(d) log d<1,

then

    A <= 4e M d^(2k_0) / [eta(1-gamma)].

This is the C106 operator reduction with the improved C109 velocity constant. Neither the required central-charge count nor the quantitative support bound is derived from modular bootstrap here. The hypothesis permits unboundedly many slow sectors.

[Proved work record] PROGRESS.md, PROOF.md, CHECKS.md, the Python solvers and verifiers, and local STATE.md are complete. Only work_C117 was modified. No packages were installed. Typhon returned "Connection timed out during banner exchange", so larger-block and larger-ring experiments remain unrun. A SLURM driver is provided in run_typhon.slurm; no heavy computation was moved to the local machine or a login node.
