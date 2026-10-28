# C51 progress

[proved, scope] The background, project memory, lead claim, expert note, C45/C46 results and the relevant proofs, C39/C42/C44 handoffs, literature review and C47 progress have been read. C48, C49, C50 and C52 have no PROGRESS.md at this checkpoint. C50 supplies a requested hypothetical input, not a completed theorem available to inspect. All writes stay in work_C51. No network or delegated agents are used.

[proved, inherited facts] C46 already derives the Gaussian thermal coefficient, the common-bond quotient pressure comparison, and the Sugawara central term. These are cited as prior results. The two-bus chain has no imposed left-right symmetry. The safe thermal coefficient sums all 2(n+1) inverse chiral speeds.

[assumed, theorem boundary] Power-law fermion correlations and a vanishing beta function alone are weaker than a complete Gaussian scaling limit or a thermal vacuum-functional estimate. The report will state separate sufficient hypotheses and will not claim that C50 has proved them.

[heuristic, work in progress] Develop a general positive multichannel Gaussian Hamiltonian with density-current mixing, prove its chiral diagonalization, track the compact charge lattice and moving irrational-density sectors, and establish the oscillator multiplicities and central terms. Give a quantitative scale-summation theorem for transferring the thermal coefficient once temperature-sensitive scalar RG estimates are available. Identify the missing microscopic estimates explicitly.

[proved, completed continuum proof] PROOF.md now derives c_R=c_L=n+1 for a complete compact Gaussian scaling limit, using a canonical diagonalization valid for arbitrary positive velocity ratios and density-current mixing. The current central-term proof works directly in nonzero charge sectors. It distinguishes current Fock multiplicities from irreducible Virasoro multiplicities.

[proved, moving sectors and thermal criterion] A uniform shifted-lattice bound proves that irrational-density zero modes do not contribute to the thermodynamic pressure, although they obstruct an uncorrected Casimir fit. The finite-size result is conditional on spectral exhaustion and o(1/L) matching. A scalar RG summation theorem proves an o(T^2) remainder from a temperature-sensitive bound with p>2 and any correction sequence tending to zero. The power-law rate and its logarithmic borderline are explicit.

[proved, contemporaneous cross-read] C48/C49/C50/C52 progress files appeared during this run and have now been read. C50 PROOF.md explicitly leaves the all-scale RG and source theorem open. Its finite-scale scalar norm does not yet provide the temperature-sensitive estimate required here. The statement that no C50 output was present applies only to the initial checkpoint, not the final report.

[numerical, checks] The canonical transformation, finite-energy Virasoro commutators at three charges, partition counts, moving-sector subtraction, shifted theta bounds and all scalar summation regimes pass. Maximum canonical error is 3.56e-14 and maximum Virasoro error is 6.22e-15. An initial numpy warning came from evaluating an unused overflowed branch of np.where; the exponent was selected before exponentiation and the rerun is clean. These controls do not test the microscopic chain.

[proved, complete] RESULT.md, PROOF.md and the local STATE.md are complete. The conditional continuum c=n+1 statement is proved. Microscopic pressure, full tower matching and complete Gaussian scaling identification remain explicitly conditional. The final file audit verifies Python syntax, JSON validity, the required deliverables and the requested prose restrictions. No parent project file was modified.
