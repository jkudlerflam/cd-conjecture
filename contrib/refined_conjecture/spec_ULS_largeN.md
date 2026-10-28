# Thermodynamic crossover scales of the SU(N) Uimin–Lai–Sutherland chain

You are working closed-book. Give full derivations, state which steps are controlled and which are heuristic, and include consistency checks. You may write and run small Python scripts (numpy, scipy only, no installs) in the current directory to check claims; report the numbers they produce. Keep runtime under 20 minutes.

## Setting

H_L = sum_{i=1}^{L} P_{i,i+1} on (C^N)^{⊗L}, periodic, where P is the swap on C^N ⊗ C^N (eigenvalues +1 on symmetric, −1 on antisymmetric pairs). The bond norm, defined as half the spread of the eigenvalues of the bond operator, is 1 for every N. At infinite temperature the entropy per site is log N and the energy per bond is 1/N. The model is Bethe-ansatz integrable (Sutherland 1975). Its low-energy theory is the SU(N)_1 WZW model, central charge c = N − 1, with one velocity v (all spinon species share it), plus a marginally irrelevant current–current perturbation.

For N = 2, P = 2 S·S + 1/2, so this is the spin-1/2 Heisenberg chain with exchange 2: v = π and e_0 = 1 − 2 ln 2 per site.

## Definitions

- s(T): thermodynamic entropy per site; s_CFT(T) = π c T/(3 v).
- For a tolerance ε (take ε = 0.1), T*_ε = sup{ T : |s(T')/s_CFT(T') − 1| ≤ ε for all 0 < T' ≤ T }. Log corrections from the marginal operator are allowed to be part of the deviation.
- Onset length ℓ* = v/T*_ε, and entropy fraction σ* = s(T*_ε)/log N.
- The exact constraint s ≤ log N forces T*_ε ≤ 3 v log N/(π (1 − ε)(N − 1)).

## Questions

1. Derive v(N) and the ground-state energy per site e_0(N) from the Bethe ansatz, with their large-N expansions. Check against N = 2.
2. Describe the elementary excitations (spinon species, dispersions, bandwidths). At what momentum and energy scales (as functions of N) does the dispersion stop being linear?
3. Large-N thermodynamics. Determine the scaling with N of T*_ε, ℓ* and σ*. In particular decide between T* ~ a/N^2, T* ~ a log N/N^2, or something else, and give the prefactor if you can. Describe s(T) in the intermediate regime T*(N) ≪ T ≪ 1. One possibly useful observation: on rings with L ≤ N the Schur–Weyl decomposition makes the totally antisymmetric irreducible representation the ground space, with GL(N) multiplicity binomial(N, L), so short segments carry entropy of order log(N/L) per site at zero temperature; the CFT entropy appears only on lengths larger than about N. Make this, or a better picture, quantitative (for example by a controlled large-N limit at fixed T N^2, or by a 1/N expansion of the thermodynamic Bethe ansatz).
4. Finite volume. For which ring lengths L (as a function of N) does the spectrum up to scaling dimension Δ ≈ c/4 + 2 take its CFT form? Heuristics allowed, if stated.
5. Give a table of your best estimates of v, e_0, T*_{0.1}, ℓ*, σ* for N = 2, 3, 4, 5, 8, and the large-N asymptotics.

Be explicit about any formula you recall rather than derive.
