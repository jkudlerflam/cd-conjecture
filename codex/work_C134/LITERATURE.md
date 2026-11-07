# C134 literature audit

2026-11-07. This separates finite-volume identities, asymptotic calculations, and uniform inequalities. The new compact theorem is in PROOF.md and is not attributed to the sources below.

## Rigorous thermodynamic and finite-size results

[Proved, literature statement] [Kozlowski, On condensation properties of Bethe roots associated with the XXZ chain, arXiv:1508.05741v2](https://arxiv.org/html/1508.05741v2), assumes even length. Proposition 3.4 supplies the finite-boundary counting expansion; Corollary 4.3 explicitly supplies the energy correction for density below 1/2. Theorem 3.5 at half filling, -1<Delta<=1, gives convergence with an accuracy-dependent threshold, not a rate. Its proof first cuts off rapidities, takes subsequential limits, and then removes the cutoff. Theorem 4.1 gives convergence of bounded Lipschitz root sums. The remark after Lemma 4.1 identifies missing uniform energy control near half filling. The leading estimates are density-uniform on the stated sub-half-filled sets at fixed Delta, not explicitly Delta-uniform.

[Proved, deduction] Full support makes the linear dressing problem solvable by Fourier transform, but removes the lower bound on density at a finite endpoint. One cannot take q to infinity in a remainder containing inverse powers of rho(q). PROOF.md replaces that step on Delta<=0 with a finite counting coordinate, a cumulative discrepancy estimate, and a cutoff at rho(A) proportional to 1/L. It does not derive the universal Casimir coefficient.

[Proved, literature statement] [Duminil-Copin, Kozlowski, Krachun, Manolescu and Tikhonovskaia, On the six-vertex model's free energy, arXiv:2012.11675v2](https://arxiv.org/html/2012.11675), Theorem 2, requires

    n <= N/2-C min{zeta^(-2), log(N)^2}

and has O(1/N) error in per-site free energy. Theorem 22, with a!=b and |Delta_6V|<1, requires n<=N/2-C/zeta and improves the error to O(1/[zeta(N-2n)N]). Both exclude n=N/2. Section 8 Claim 1 displays a quadrature denominator N-2n-C_0. Their convention Delta_6V=-cos zeta is opposite to the antiferromagnetic-transverse quantum convention used here.

[Proved, deduction] These estimates do not provide an O(1/N^2) per-site error at half filling. Even such an estimate for a transfer eigenvalue would require a uniform bound on its spectral-parameter derivative before taking the Hamiltonian limit. Commutation and entrywise positivity alone supply neither derivative control nor the positive spatial operator assumed in C125.

## Nonlinear integral equations and classical finite-size calculations

[Proved, literature statement] [Destri and de Vega, Unified approach to thermodynamic Bethe ansatz and finite size corrections for lattice models and field theories, hep-th/9407117v2, Nucl. Phys. B 438 (1995) 413](https://arxiv.org/html/hep-th/9407117), Sections 2-3, derive finite-size NLIEs and energy/transfer-eigenvalue quadratures. Section 3 calls the NLIE exact at finite N. In Section 7.1 the bulk region is discarded to orders 1 and 1/N, rapidities are shifted by (gamma/pi)log N, and a limiting equation is substituted into the energy integral. The dilogarithm lemma then evaluates that limiting integral.

[Assumed, unresolved estimate] The step needed for this task is a quantitative bound on the discarded bulk and the shifted tails, or a dominated-convergence estimate strong enough to retain an O(1/N) energy remainder uniformly in gamma. Section 7.1 does not supply such an explicit norm bound. This is a limitation of the inference to the present theorem, not a claim that the finite-N contour identity or dilogarithm algebra is false. PROOF.md bypasses this limit exchange on Delta<=0.

[Proved, source statement with access limitation] [de Vega and Woynarovich, Method for calculating finite size corrections in Bethe ansatz systems, Nucl. Phys. B 251 (1985) 439](https://doi.org/10.1016/0550-3213(85)90271-8), describes an exact contour representation and a saddle-point treatment yielding leading finite-size corrections. The publisher abstract was accessible through search; the full paper was not accessible. I therefore do not attribute a specific missing lemma or a rigorous uniform remainder theorem to its unread body.

[Proved, source statement with access limitation] [Hamer, Quispel and Batchelor, Conformal anomaly and surface energy for Potts and Ashkin-Teller quantum chains, J. Phys. A 20 (1987) 5677](https://doi.org/10.1088/0305-4470/20/16/040), uses XXZ equivalences and finite-size methods to derive leading corrections. The abstract and indexed paper passages were accessible, but the complete proof was not. This is not used as a uniform finite-L inequality. The distinct Hamer paper, J. Phys. A 19 (1986) 3335, treats XXZ sectors; the author list and date should not be conflated.

[Proved, source statement with access limitation] [Kluemper, Batchelor and Pearce, Central charges of the 6- and 19-vertex models with twisted boundary conditions, J. Phys. A 24 (1991) 3111](https://doi.org/10.1088/0305-4470/24/13/025), derives NLIEs and extracts central charges using Rogers dilogarithms. The abstract and indexed opening of the PDF were checked. Full-PDF retrieval timed out or was unavailable. No line-by-line rigor verdict on that paper is claimed here.

[Proved, author's account] [Kozlowski's RAQIS 2016 slides](https://lapth.cnrs.fr/conferences/RAQIS/RAQIS16/slides/Kozlowski.pdf), in the indexed slide on the counting function, list the traditional working inputs as bounded roots and a positive counting derivative, or a priori root growth and local derivative control. The slide explicitly notes that the NLIE form depends on these inputs. Full-PDF retrieval returned a service error, so this is an indexed-slide check, not a complete slide-deck review.

## Later rigorous work and special exact energies

[Proved, literature statement] [Goehmann, Goomanee, Kozlowski and Suzuki, Thermodynamics of the spin-1/2 Heisenberg-Ising chain at high temperatures: a rigorous approach, arXiv:1811.12020](https://arxiv.org/html/1811.12020), Theorem 1.1, proves an NLIE free-energy representation for T>T_0. The paper establishes dominance and the relevant limit interchange in that regime. It is a rigorous result about high-temperature thermodynamics, not a zero-temperature finite-ring Casimir estimate.

[Proved, literature input] Exact supersymmetric energies supply the isolated point Delta=+1/2 in our convention. The sign conversion, twist, cited equations, and deduction of C=2pi^2/9 are recorded in PROOF.md Section 11. The cited papers' Delta=-1/2 has negative transverse exchange, so using that label unchanged in the present Hamiltonian would put the result at the wrong anisotropy.

[Assumed, literature completeness unresolved] Searches for later rigorous half-filled finite-size energy estimates, including Kozlowski's later work and six-vertex papers, did not identify a theorem with the requested uniform explicit bound on the entire massless interval. This is a search outcome, not proof that none exists. The inaccessible classical papers above limit the scope of the audit.

## Consequences for the proposed routes

[Proved] Route (a) supplies exact representations but still needs controlled cancellation and tails. Route (b) can be completed for the one-sided estimate on [-1+delta,0] by the real-axis argument in PROOF.md. Route (c), using only the stated six-vertex errors or C125 monotonicity, loses the needed factor of L. Route (d) would follow from a bound uniform up to Delta=1 by finite-dimensional continuity, without an additional monotonicity assumption. Such a bound has not been established in this run.
