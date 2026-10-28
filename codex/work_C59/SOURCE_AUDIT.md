# C59 primary-source audit

2026-10-28. Sources were accessed with the browsing tool. No shell downloads or network commands were used. No author was contacted.

## 1. Search for the announced thermal result

[read, announcement only] The official [127th Rutgers short-talk program, page 10](https://math.rutgers.edu/images/127_SMC_short_talks.pdf#page=10) lists Bruno Renzi's talk with Alessandro Giuliani, Vieri Mastropietro and Venanzio Alessio Santini. Its title is "Free energy universal finite temperature correction in non-solvable quantum chains". The program dates this session to December 17, 2024. The abstract announces an RG proof for weakly interacting quantum chains. It supplies no Hamiltonian, theorem hypotheses, ensemble, dressed-velocity definition or error bound.

[read, search result] Searches for the title, the author combinations, and variants using thermal, finite temperature, quantum chains, central charge and free energy did not locate a corresponding 2024-2026 arXiv manuscript. ArXiv-domain searches and university-domain searches were included. This is a report of the search, not a proof that no manuscript exists.

[read, unrelated 2026 paper] [Fabbri-Mastropietro-Renzi, arXiv:2602.03705](https://arxiv.org/abs/2602.03705), "Non-perturbative renormalization for lattice massive QED2: the ultraviolet problem", was located. The [full-text introduction](https://arxiv.org/html/2602.03705v1) concerns a massive gauge theory with Wilson lattice regularization and ultraviolet convergence. It is not the announced quantum-chain thermal paper. No theorem from it is invoked for XXZ or the matched retarded action. In particular its abstract's integration of massive bosons is not a theorem for every supplied retarded endpoint class.

[read, search limitations] Giuliani's publication list appeared in indexed search results, but opening its main page timed out. An alternate university-host URL and the attempted Mastropietro homepage were inaccessible. Direct arXiv author-search pages returned cache misses. The SISSA seminar archive was opened and its Renzi entry concerns interacting dimers. These failures limit the negative search conclusion. Search snippets were not treated as complete theorem statements.

[unverified scope] There is no located thermal manuscript whose theorem can be stated or applied exactly. Coverage of off-half-filled XXZ, both signs of small anisotropy, fixed field versus fixed density, or arbitrary retarded actions cannot be certified from the announcement. C57's announcement lead remains unresolved.

## 2. The Ising comparison

[read] [Giuliani-Mastropietro, arXiv:1211.3030](https://arxiv.org/pdf/1211.3030), (3.61), (3.68)-(3.73), (3.84)-(3.91), was checked directly. Equation (3.90) reads, in j<=0 notation,

    C |lambda| 2^((2+theta/2)j) (2^(-j)/ell)^(3-epsilon).

Its free-Majorana reference and decaying local couplings are substantive model restrictions. PROOF Section 2 explains the comparison boundary.

[proved application boundary] The joint estimate in PROOF Section 7 is the missing replacement when the reference quartic survives. PROOF Section 1 handles XXZ spatial boundary conditions independently by a bounded end-bond comparison. It does not transfer the Ising scalar theorem.

## 3. Zero-temperature XXZ matching

[read] [Benfatto-Mastropietro, arXiv:0907.2837](https://arxiv.org/pdf/0907.2837), (81)-(87), (92)-(95), supplies the cited coupling convergence and field-normalization matching. The quartic limit is nonzero in general; the anomalous wavefunction exponent is retained. Theorem 1.1 concerns exponents and susceptibility. No finite-temperature vacuum estimate is extracted from it here.

[proved scope] PROOF Section 3 separates these localized zero-temperature data from a full norm bound on the effective lattice/reference defect. Section 4 gives a separate smooth-symbol argument for Gaussian scalar factors. Neither substitutes for the uniform marked GN estimate.

## 4. Repository inputs

[read, imported results] C57's complete proof and source audit, C58's result, C50-C56 summaries and C55 Lemma N were consulted. The collaborator's thermal criterion, winding estimate, thermal localization, moment reserve, quadratic-memory example, ultraviolet audit, single-channel reduction, infrared-source review and multivelocity review were used to identify what is already established. Their single-step estimates and reference oscillator formula are cited as inputs rather than presented as new all-scale theorems.

[proved provenance boundary] C59 proves the explicit lemmas in its PROOF, not the nonlinear RG input in C55 or the thermal comparison it asks for in TP. The new derivations have not had independent or human review.
