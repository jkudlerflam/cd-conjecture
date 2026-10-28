# C57 primary-source audit

[read] Sources below were opened through the browsing tool on 2026-10-28. Equation numbers refer to the linked versions. Search results alone are not treated as theorem statements. No source was downloaded with a shell command.

## 1. A directly relevant announcement

[read, announcement only] Bruno Renzi, with Alessandro Giuliani, Vieri Mastropietro and Venanzio Alessio Santini, is listed for short talk B12 at the 127th Rutgers Statistical Mechanics Conference, 17 December 2024. Its title is "Free energy universal finite temperature correction in non-solvable quantum chains". The [official abstract, page 10](https://www.math.rutgers.edu/images/127_SMC_short_talks.pdf#page=10) announces a fermionic RG proof of the leading thermal correction for a class of weakly interacting quantum chains.

[unverified scope] The abstract states neither the Hamiltonian nor the central-charge value, parameter range, ensemble, remainder, or shell estimate. In particular it does not identify XXZ as a covered model. Searches by title, authors, and thermal/central-charge keywords did not locate a full manuscript. An attempted opening of Giuliani's publication page timed out. This is a concrete prior-work lead, not evidence that an available theorem has been checked. It prevents a responsible claim that nobody has addressed this problem. Obtaining the manuscript is the most useful literature follow-up. No message to any author was sent.

## 2. The spinless Benfatto–Mastropietro result

[read] Benfatto–Mastropietro, [Universality relations in non-solvable quantum spin chains](https://arxiv.org/pdf/0907.2837), J. Stat. Phys. 138 (2010), 1084–1108. Theorem 1.1, equations (14)–(18), identifies exponents and susceptibility. The useful flow estimates are (81), the constant-history beta bound, and (82),

    |lambda_j-lambda_-infinity| <= C_theta lambda^2 gamma^(theta j),
    j<=0, theta<1.

Equations (83)–(87) control field/source normalization ratios. Equation (23) distinguishes the prescribed Fermi point, the physical field, and the sound-velocity counterterm.

[application boundary] These estimates provide a vanishing infrared factor for particular normalized couplings. They do not state a thermal difference of vacuum scalars with power greater than two. Theorem 1.1 is not used here as a thermal theorem. Its first-order K formula (16) and susceptibility identity (18) are used only in the perturbative consistency check in PROOF.md, Section 6.

## 3. The Drude-weight companion

[read] Benfatto–Mastropietro, [Drude weight in non solvable quantum spin chains](https://arxiv.org/pdf/1012.0403), J. Stat. Phys. 143 (2011), 251–260. Theorem 1.1, equations (14)–(16), gives a current response with a Holder remainder and D=uK/pi, u^2=D/kappa. This identifies a physical sound speed independently of thermal pressure. It does not prove a heat-capacity or T^2 free-energy coefficient. The abstract's statement about completing the Luttinger-liquid conjecture must be read with these explicit conclusions.

## 4. Free-energy construction for spinning fermions

[read] Benfatto–Falco–Mastropietro, [Universality of one-dimensional Fermi systems, I](https://arxiv.org/pdf/1303.3681), CMP 330 (2014), 153–215. Section 3.1, equations (3.1)–(3.4), retains the scalar shell contributions and proves the zero-temperature/thermodynamic limit. The majorant used there is |E_j^(beta,L)|+|E_j| <= C gamma^(2j). Section 3.6, equations (3.22)–(3.24), establishes the analytic/Borel control of the ground-state expansion.

[application boundary] A convergent scalar expansion is present, but its displayed bounds do not determine the coefficient after division by T^2. This is a spinful model with its own repulsiveness assumptions, not a direct spinless XXZ application. [Part II](https://arxiv.org/pdf/1303.3684), CMP 330 (2014), 217–282, concerns response/exponent identities and the effective-model comparison. Neither is invoked here to supply estimate P.

## 5. The closest published scalar estimate

[read] Giuliani–Mastropietro, [Universal finite size corrections and the central charge in non solvable Ising models](https://arxiv.org/pdf/1211.3030), CMP 324 (2013), 179–214. With q_j=2^j, circumference ell and their notation for finite/infinite-volume histories, equation (3.90) states

    |e_j^infinity(v_tilde)-e_j^infinity(v_*)|
      <= C |lambda| q_j^(2+theta/2) (1/(ell q_j))^(3-epsilon).

The paragraph after (3.91) gives the same bound for the fixed-history geometric difference. Equation (3.69) controls histories, (3.84) bounds zero-temperature scalars, (3.87) the missing infrared tail, and (3.91) sums the correction to O(|lambda| ell^(-2-theta/2)).

[application boundary] This is the published analog of P, with p=3-epsilon>2 and eta_j proportional to q_j^(theta/2). It applies to their perturbed Ising model. Its local four-fermion interaction vanishes; the limiting theory is a free Majorana field. Spinless XXZ instead has a surviving right-left quartic. Replacing the Ising reference by the interacting Luttinger reference changes a substantive part of the comparison. This paper supplies a proof pattern, not an XXZ theorem used by C57.

## 6. Multichannel papers

[read] Mastropietro–Porta, [Multi-Channel Luttinger Liquids at the Edge of Quantum Hall Systems](https://arxiv.org/html/2108.13708v3), Section 8.2, especially (8.55), constructs an analytic free energy from bounded normalized running data. Theorem 8.5, (8.61), states

    |beta_j^v(s,...,s)| <= C |s| 2^(theta j),
    |beta_j^lambda(s,...,s)| <= C |s|^2 2^(theta j).

Section 10 compares lattice and reference beta functions. These are not bounds on a thermal vacuum difference. Their utility here is the reference cancellation and the known role of normalized field ratios, not a claimed evaluation of the scalar coefficient.

[read] Porta–Scola–Singh, [Large Scale Dynamical Response of Interacting 1d Fermi Systems](https://arxiv.org/html/2509.08665v1), Theorem 2.6 and Section 6, establish response for their stated instantaneous lattice model. The multichannel framework permits unequal nonzero branch velocities. Its theorem is not an arbitrary retarded-endpoint stability result or a thermal-pressure theorem. The local C50 and collaborator source audits give the detailed hypothesis ledger; C57 does not duplicate it.

## 7. Integrability does not automatically remove the rigor question

[read] Faulmann–Gohmann–Kozlowski, [Low-temperature spectrum of the quantum transfer matrix of the XXZ chain in the massless regime](https://arxiv.org/pdf/2305.06679). Their Theorem 1.1 constructs the candidate dominant Bethe eigenvalue and compares it within the analyzed solution class. Theorem 1.4's identification with the physical free energy assumes Conjectures 1.2 and 1.3, concerning the dominant eigenvalue and exchange of Trotter/thermodynamic limits. Section 7.6.2 explains this boundary.

[application boundary] This source cannot be cited as an unconditional thermal theorem for the chain by omitting those hypotheses. This audit does not prove that all other Bethe-ansatz approaches have the same gap. It also does not dispute the standard XXZ thermal prediction. Even an independent integrable proof of the total coefficient would not supply P or stability under a nonintegrable retarded correction.

## 8. Answer to which published estimate is used

[proved, scope of this audit] No verified published estimate is used here to assert P for interacting XXZ. The closest actual scalar estimate is GM (3.90), with the model restriction above. BM (81)–(82), BFM (3.1)–(3.4), and MP (8.61) control different necessary pieces. The Rutgers announcement may contain the desired missing theorem, but its proof and hypotheses were not available in the located material. C57's mathematical results are the conditional comparison and first-order calculation in PROOF.md, not a completed new all-order theorem.
