# C133 result

[heuristic, adversarial verdict] No fixed-d fast tower was constructed. No general no-go theorem was proved. The concrete controlled encoder is obstructed, while hierarchical emergence, weak charge visibility and genuinely neutral interacting encodings remain open. Failure to supply their microscopic bonds is not evidence that fixed-d locality forbids them.

[proved, correction to the target framing] An unbounded number of order-M velocities is sufficient, but not necessary, to violate the bound. The sector data c_k=1, v_k=M/k give A/M=H_N, although the number above any fixed positive velocity threshold stays uniformly bounded. A theorem about fast sectors therefore needs an integrable bound as the threshold tends to zero. Also, bounding infrared A alone does not automatically bound the full finite-size supremum C(d).

## 1. Quantitative diagnosis of the fixed-d construction

[proved, leading-reference theorem] For n encoded payload channels, cell lengths n+3 and n+4, core scale s, port coupling lambda, bus mass D and conversion g, define delta=a-2w>0. The inherited hopping matrix is

    T=-W[Gamma(w)+i eta Gamma(w v_bus)],
    W=(2 lambda^2 g^2/s) V_n D^(-2n-4).

Two new bounds are

    V_n <= (n+2)^(-2),
    tr Gamma(w) <= 1/(4 delta^3).

The second follows by summing the nonnegative bus filters before squaring them, using the heat-kernel semigroup identity. With the physical length conversion included,

    A_payload/M <= [(n+4)/(n+2)^2]
                   lambda^2 g^2(1+eta)/(s M delta^3 D^(2n+4)).

An isolated microscopic express hop certifies M>=s. In the controlled range lambda sqrt(1+eta)<=kappa s,

    A_payload/M <= kappa^2 g^2(n+4)
                   /[delta^3(n+2)^2 D^(2n+4)].

Each payload velocity satisfies the same upper bound times M. PROOF.md also includes the length mode. These statements concern the leading Gaussian reference. Their transfer to the full scheduled chain requires S and relative control of dressed velocities.

[proved diagnosis] The losses come from silent first-order ports, two conversion events, propagation through 2n+4 high vertices, and endpoint visibility. The stability schedule adds further suppression. The particular extreme exponents are not universal. Encoding capacity alone gives only n<=b log_2 d and does not prove slow propagation. The cells themselves do not have bounded dimension. A further attempt to remove the bias and lower the bus masses runs into the unbiased path gap of order s/n^2: when that excitation is retained and the port remains perturbative, the effective physical speeds are at most order M/n. PROOF.md section 2.5 states the spectral hypothesis and its proof.

[proved conditional on Gaussian matching and the stiffness dictionary] There is a separate obstruction to accelerating this charged weak-band phase. Its physical on-site charge index is S_Q=101664/31. If matching changes the current quadratic form and total reference speed by relative epsilon<1, then

    A_payload <= 4pi S_Q (1+epsilon)/(1-epsilon) M.

Thus an acceleration that remains uniformly close to this charged Gaussian reference cannot give unbounded A/M, independently of the particular bus losses. The one length boson adds at most 4eM. A small-K counterexample must leave this uniformly weak matching regime.

## 2. Status of the attempted towers

| Candidate | Result and remaining obstruction |
| :-- | :-- |
| Hierarchy xi_k=b^k with v_k of order M | [heuristic, unresolved] Thermodynamic entropy capacity permits it. PROOF.md constructs bounded entropy profiles with an arbitrarily large infrared slope. No single fixed-d repeated bond realizing the hierarchy was obtained. Free, direct-stack, first-order Gaussian and controlled-bus realizations are separately obstructed. |
| Many weakly charged bosons | [proved within the continuum dictionary] K_k=kappa/k^2 keeps the diagonal-charge stiffness bounded while A=Nv0 grows. Suitable incommensurate momenta can avoid the finitely many relevant winding vertices at each N. [heuristic, unresolved] No fixed-d microscopic realization was obtained. Charge rescaling alone cannot produce the effect. |
| Neutral or coset sectors | [proved restricted obstruction] One minimal model or Z_p parafermion cannot have unbounded c. Full symmetry-matched WZW parents are covered by C104, including their same-velocity cosets. [heuristic, unresolved] Multiple neutral interacting factors or a coset without such a matched parent remain unbounded by the available argument. |
| Literature constructions | [proved/source audit] The explicit black-hole spin chain is non-Hermitian and has Casimir c_eff=2 despite large Virasoro c. [heuristic] Staircase plateaus are successive regimes, not additive sectors. Direct high-rank WZW implementations increase d. No in-scope counterexample was identified. |

[proved, restricted locality lemma] For an exact first-order product-block compression whose retained Hamiltonian is quadratic and number conserving,

    rank T <= floor((d^2-1)/2),
    number of complex gapless sectors <= rank T.

The first inequality uses the operator-Schmidt rank of one physical cut. The second uses the degree of det(F+zT+z^(-1)T^dag). It permits arbitrary block length and arbitrary onsite mixing. A second-order resolvent can already have arbitrarily large rank, so the lemma does not extend by counting port insertions. Its proof and limitation are in PROOF.md section 3.

## 3. Numerical tests and their limits

[numerical] Computed 98 unrestricted ring ground states for 21 parameter choices in XXZ, Ashkin-Teller, spin-orbital, spin-1 bilinear-biquadratic and FZ parafermion families, d=2 through 6. The largest space had dimension 65536. The full table includes Q_L, 4Q_L/pi, fitted A/M where enough sizes were available, and pi H_(d-1). M is gauge optimal in these reflection-symmetric families. See NUMERICS.md.

[numerical] The XX and Ising exact-energy controls agree within 1.8e-14. Sparse and matrix-free implementations agree within 5.4e-15 on their controls. The largest eigenvector residual in the scan is 4.1e-10. Gram-integral and endpoint-visibility checks pass. These residuals are convergence checks, not interval-certified ground-energy bounds.

[numerical, interpretation] No unbounded growth signal was found. Short-ring values can exceed the proposed asymptotic envelope. For example, ULS gives A_eff/M=5.326 at L=3, compared with the envelope 4.712 and its known infrared A/M=4pi/3. Heisenberg's fitted value 3.157 is slightly above its exact pi. XXZ at Delta=1.05 still fits near-critical behavior despite being asymptotically gapped. These are finite-size artifacts, not violations inferred from this scan.

[numerical] The C106 two-band family was also tested on vacuum-compatible rings using both physical fermion parity sectors. In one scan v/M is held at 1/2 as the crossover length increases. Its exact A/M stays 1, and the estimates at L=64 ell run from 1.000243 to 1.000080 through ell=128. The largest doubled ring has 16384 sites; this is a small single-particle sum, not a large many-body computation. See CHECKS.md.

[numerical, limited mobile test] Nine further unrestricted d=6 mobile-toy energies at L=3,4,6 agree with the complete legal-sector minima. Gauge M is bracketed. These short rings do not supply a certified vacuum-compatible charge sequence, so MOBILE.md reports raw Q intervals and does not assign an infrared A. This toy has only one payload and one length mode, with no unbounded-channel tuning parameter.

[proved, scope of the numerical work] The abstract hierarchy and small-K proposals have no microscopic bond to diagonalize. The n=2 projected route-A local space is already 8-dimensional, and its microscopic d0=4464 chain is outside the requested small-d scan. Their infrared coefficients were not claimed measured. No DMRG, Typhon run or full microscopic route-A simulation was performed.

## 4. Missing theorem and handoff

[assumed, sufficient missing lemma] A crossover-independent bound on the weighted count

    N_c(u)=sum_(v_a>=u)c_a <= B(d)(M/u)^p,  0<p<1,

would give

    A/M <= B(d)(4e)^(1-p)/(1-p).

An arbitrary explicit integrable function of u/M would also suffice. A 1/u bound would not. The missing step is a microscopic estimate on all dressed boundary-channel weights, including neutral and weakly charged sectors, through arbitrarily long crossovers. It is not supplied by entropy capacity, stiffness alone, a finite operator count or the first-order rank lemma.

[proved, lesson audit] PROOF.md section 8 explains how the result respects lessons 1 through 5. No transfer-positivity-only proof, local-response proxy, bounded crossover, raw state-counting argument, renamed reconnection assumption or eventual convexity was used.

[proved, reproduction] Run from this folder with NumPy and SciPy:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B scan.py
    python3 -B summarize.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B mobile_probe.py

[proved, pending computation] large.sbatch is a SLURM-ready five-family array for larger rings using the matrix-free implementation where needed. It was not submitted because this sandbox cannot reach Typhon. It assumes an existing Python environment with NumPy and SciPy. Larger ED may clarify the finite-size drift of these controls, but cannot test a proposed tower until a microscopic bond is specified.

[proved, file record] All changes are confined to work_C133. STATE.md is the local handoff because the project-level STATE.md is outside the authorized write scope. No repository-wide files, commits or pushes were changed.
