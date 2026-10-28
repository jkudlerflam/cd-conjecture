# Literature check of the integrate-out-then-bosonize step

## Claim

The step is the one used in the candidate counterexample to the c-d conjecture (Latorre and Sierra, 2403.17242). A ring of M identical cells has a gapped sector, with gap Δ above an exactly degenerate code space. At fixed small coupling λ, one removes the gapped sector by local Schrieffer-Wolff (SW) or degenerate perturbation theory, obtains a quasi-local effective Hamiltonian on the code space, and concludes that the infinite chain flows to a multi-component Luttinger liquid with n flavours, widely separated velocities, and c = n.

No published theorem covers this. Rigorous tools control energy densities, finite times, or gapped targets. Physics checks reach four flavours and two velocity scales. The documented failures involve marginal flows or symmetry-allowed relevant operators, which the design (exact on-site Z₂, detuned Fermi points) is meant to exclude.

Tags. [read] means I read the cited passage. [abstract] means abstract only, and a star means the abstract was seen only in a search-engine summary.

## Status

Heuristic, a literature survey only. Section 4 states the theorem that would upgrade the step.

## How it was established and checked

- Author: Claude (Sonnet 5.5), sub-session of the coordinating session, 2026-10-28.
- AI assistance: none beyond the author.
- Checked by: nobody yet. I read passages of 1701.05182 (Definition 23, Lemmas 37 and 39), 2003.13753 (§4), 2303.07386 (§7) and the first page of Ogata and Shiba, and 2509.20431 through a page extract. Papers [read] in `contrib/ir_control/` are not re-read.

## Argument or method

### 1. Rigorous control

Hamiltonian simulation. Cubitt, Montanaro and Piddock, 1701.05182 [read], define a (Δ, η, ε)-simulation through an operator-norm error ε on the spectrum below a cutoff Δ, and their gadget lemmas need Δ polynomial in ‖H‖, 1/η, 1/ε. Since ‖H‖ is extensive, the simulator strengths grow with system size, and the error is absolute, not per site. Kohler, Piddock, Bausch and Cubitt, 2003.13753 [read], get translation-invariant 1D universality from history states, with strengths polynomial in size, target norm and 1/ε, so it is efficient in spins but not in norm. Aharonov and Zhou, 1804.11084 [abstract], show degree reduction needs growing strengths. Bravyi, DiVincenzo, Loss and Terhal, 0803.2686 [abstract], give n-independent couplings, but the ground-state energy is reproduced only up to an extensive error O(εn). So no version works at fixed couplings in the thermodynamic limit at the 1/M resolution a CFT identification needs.

Locality of SW. Bravyi, DiVincenzo and Loss, 1105.0675 [read, per ir_control], give local SW with an extensive error, and degenerate local ground spaces leave no global gap. Yin and Lucas, 2209.11242 [abstract], prove that for a gapped local H₀ plus a local perturbation of size ε ≪ Δ, dynamics of local observables stays in the low-energy subspace of H₀ up to times exp[(Δ/ε)^a], a < 1/(2d−1), whether or not the perturbation closes the gap. The review by Chen, Lucas and Yin, 2303.07386 [read, §7], explains the proof by an iterated quasi-local SW. It states that a unitary block-diagonalizing completely cannot be quasi-local, because of the orthogonality catastrophe, so the iteration stops at k* ~ (Δ/ε)^a with an exponentially small but nonzero remainder. By my estimate, comparing that remainder with the 1/M level spacing limits spectral statements to M up to about exp[(Δ/λ)^a/3], a prethermal window and not the thermodynamic limit.

Exact block diagonalization. Fröhlich and Pizzo, 1812.02457 [abstract], and Del Vecchio, Fröhlich and Pizzo, 2108.13907 [abstract], use local Lie-Schwinger conjugations to prove the gap persists uniformly in length, so they need a gapped infrared. Datta, Fernández and Fröhlich, math-ph/9809007 [abstract], build volume-independent unitarily equivalent effective Hubbard-type Hamiltonians by convergent expansions, with Pirogov-Sinai theory for the classical part. Their predecessor with Rey-Bellet (Helv. Phys. Acta 69, 1996, no arXiv) [abstract*] fixes the low-energy spectrum only to some finite order beyond one. Ueltschi, cond-mat/9810320 [abstract], proves an analytic Hubbard free energy only for βt²/U small.

Constructive RG. Mastropietro and Porta, 2108.13708, and Porta, Scola and Singh, 2509.08665, stay the closest frameworks. They control correlations of fermionic models, not spectra or c, and I found no paper applying them to a model obtained by integrating out a gapped sector.

### 2. Independent confirmations

Large-U Hubbard is the cleanest case. Ogata and Shiba, PRB 41, 2326 (1990), no arXiv [read, first page], show that at U to infinity the Lieb-Wu wave function is a spinless-fermion Slater determinant times a spin-1/2 Heisenberg wave function, the exact version of the SW statement. The finite-U picture, with v_s of order J = 4t²/U and v_c of order t, follows from the Lieb-Wu solution (cond-mat/0207529 [abstract]) and the finite-size spectrum of Frahm and Korepin (PRB 42, 10553, 1990, no arXiv) [recalled, not retrieved]. That is two flavours, with v_c/v_s growing as U/t.

Assaraf, Azaria, Caffarel and Lecheminant, cond-mat/9903057 [abstract], compare bosonization of the SU(N) Hubbard chain with Monte Carlo for N = 3 and 4, with a Mott transition at nonzero U. A DMRG study, 2412.13320 [abstract*], finds c = N−1 in the Mott phase and c = N in the metal from entanglement entropy. Sikkema, Affleck and White, cond-mat/9702143 [abstract], treat the doped Kondo chain by bosonization at weak coupling and a modified t-J model at strong coupling. Khait et al., 1710.04847 [abstract], find a heavy Luttinger liquid there by DMRG, with RKKY effects crucial, unlike in large-N mean field.

Scope. The largest case is n = 4 with equal spin velocities. I found no check with many flavours, hierarchical velocities, or many-body excited cells at fixed λ.

### 3. Failures and warnings

I found no case where a controlled SW with a robustly critical, non-marginal target failed. These fail through marginal flows or relevant operators the symmetry does not forbid.

- Itoi and Kato, cond-mat/9605105 [abstract*]. The SU(3) point of the spin-1 chain is the Kosterlitz-Thouless boundary between the gapped Haldane phase and an extended critical phase, and the sign of an SU(3)-breaking marginal operator decides which side is reached.
- Itoi, Qin and Affleck, cond-mat/9910109 [abstract*], with Azaria et al., cond-mat/9903047 [abstract*], and Boulat et al., cond-mat/9910218 [abstract*]. For the SU(4) spin-orbital chain, a small deviation opens a gap, exponentially in one RG treatment and as a power law in DMRG, and marginal couplings blur the boundary.
- Prembabu and Verresen, 2509.20431 [read, §IV.4.2 and IV.5]. SW in the coupling shows no relevant operator at orders one and two, and a symmetric relevant operator of dimension 8/5 appears at third order, which would destabilize the c = 8/5 line. DMRG sees no effect at J = 0, so this is unresolved.
- Oshikawa and Affleck, cond-mat/9706085 [abstract]. A symmetry-allowed operator of dimension 1/2 and amplitude h opens a gap proportional to h^{2/3}, so one missed term is enough.
- Allen, Essler and Nersesyan, cond-mat/9907303 [abstract*]. Weakly coupled critical chains in a frustrated ladder flow, through a marginal interchain coupling, to spontaneous dimerization with a gap exponentially small in the coupling. The parent is gapless, so this is an RG warning, not an SW one.

Gaplessness at fixed λ needs every symmetry-allowed relevant operator forbidden at all orders, which a finite-order SW calculation cannot show.

### 4. Sharpest theorem that would justify the step

The unitary route is blocked, since sectors have no global gap and the remainder is only exponentially small. The natural form is Euclidean.

Theorem D (imaginary-time decoupling). Fix the cell data and the target logical chain. There is λ₀ > 0, independent of M and of inverse temperature β, such that for λ < λ₀ the partition function with insertions of code-projected local operators equals a factor Z_exc times an integral over logical Grassmann fields against exp(−S_eff). S_eff has a Gaussian part with 4n simple Fermi points within O(λ) of the target. Its other terms are even, translation-invariant kernels that decay exponentially in space and in imaginary time, with time rate at least cΔ, and weighted norm O((λ/Δ)^θ), uniformly in M and β.

Two further steps finish the argument. First, extend 2108.13708 and 2509.08665 to U(1)-free, general even interactions with sum-type non-resonance, giving multi-velocity Luttinger asymptotics of logical correlations. Second, obtain c, either as c ≥ n from the 2n chiral currents under the c-d premise, or from a spectral theorem for rings with corrections o(1/M).

D needs only e^{−Δs} decay of cell excursions, not a spectral gap between sectors, and it feeds the constructive RG. The special case of one free auxiliary fermion per cell is Astra B §D to §E in `contrib/ir_control/`. I have not checked that D is provable for general cells. The excursion-polymer expansion named there is the missing mathematics, and the sharpest existing statement is the Yin and Lucas prethermal theorem.

## Reproducing

No scripts. PDF text was extracted with pdftotext in the session scratchpad.

## Relation to the rest of the repository

This addresses H6 and H7 of `contrib/ir_control/` and the missing rigorous SW reduction in STATUS step 7. "Where things stand" could add that no literature theorem covers the step and that section 3 reinforces the Z₂ and detuning requirements.
