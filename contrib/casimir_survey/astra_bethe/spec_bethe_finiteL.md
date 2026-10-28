# Finite-size Bethe ansatz (ideal strings) for the rank-k symmetric SU(N) chain

Write ONE self-contained Python 3 script (numpy and scipy only), `bethe_strings.py`, that computes the ground-state energy E_0(L) of a periodic ring from the nested Bethe equations in the ideal k-string approximation, validates it against the exact-diagonalization numbers below, and then produces the large-L data described at the end. Return the script in full plus its printed output.

## Model and normalization

Each site carries the irrep Sym^k(C^N) of SU(N). On two sites, Sym^k ⊗ Sym^k = ⊕_{j=0}^{k} V_j with V_j the irrep of Young diagram (2k−j, j). The bond is h = Σ_{j=0}^{k} H_{k−j} Π_j, H_m = Σ_{l=1}^{m} 1/l the harmonic numbers (H_0 = 0), Π_j the projector on V_j. So the fully symmetric component (j = 0) has energy H_k, and (k,k) has energy 0. For k = 1 this is h = (1 + P)/2 with P the swap. For N = 2 it is the Takhtajan–Babujian spin-k/2 chain with v = π/2. H = Σ_i h_{i,i+1} on a ring of L sites, L a multiple of N.

Expected thermodynamics (to be tested, not assumed): e_0 = H_k − (1/N) Σ_{j=1}^{k} [ψ(1 + (j−1)/N) − ψ(j/N)] (ψ digamma), and E_0(L) = e_0 L − π c v/(6L) + o(1/L) with v = π/N, c = k(N²−1)/(N+k).

## Bethe ansatz (determine signs and conventions by matching the ED data)

Write e_n(x) = (x + i n/2)/(x − i n/2). Reference state: fully polarized in colour 1, energy L·H_k. Nested roots λ^{(r)}_a, r = 1, …, N−1, with M_r roots at level r. Standard nested equations:

e_k(λ^{(1)}_a)^L = Π_{b≠a} e_2(λ^{(1)}_a − λ^{(1)}_b) Π_c e_{−1}(λ^{(1)}_a − λ^{(2)}_c),

1 = Π_{b≠a} e_2(λ^{(r)}_a − λ^{(r)}_b) Π_c e_{−1}(λ^{(r)}_a − λ^{(r−1)}_c) Π_c e_{−1}(λ^{(r)}_a − λ^{(r+1)}_c), r ≥ 2,

(up to the overall convention for which side each factor sits on; fix it so that k = 1, N = 2 reproduces the Heisenberg chain). Energy:

E = L H_k − (1/2) Σ_a k/((λ^{(1)}_a)² + k²/4).

(Check: one level-1 root is a single colour defect with E − L H_k = −(1/k)(1 − cos p), e^{ip} = e_k(λ).)

Ground state (singlet, L multiple of N): at every level r, n_r = L(N−r)/N ideal k-strings λ + i(k+1−2j)/2, j = 1..k, so M_r = k n_r. Reduce to real string centres with the usual Takahashi phases θ_n(x) = 2 arctan(2x/n) (θ_0 ≡ 0), Θ_{nm}(x) = θ_{|n−m|} + 2θ_{|n−m|+2} + … + 2θ_{n+m−2} + θ_{n+m} for same-level string–string scattering, the adjacent-level analogue θ_{|n−m|+1} + θ_{|n−m|+3} + … + θ_{n+m−1}, and the source phase of a k-string from e_k^L, together with consecutive symmetric (half-)integer quantum numbers. Derive these reduced equations yourself, write them in the script's docstring, and solve them by Newton's method (dense Jacobian is fine). The energy of a k-string centre is the sum of the member energies.

## Validation targets (exact diagonalization, this normalization)

E_0(L):
- N=2,k=1: L=8: 2.3489065911, 12: 3.6126090826, 16: 4.8577036394
- N=2,k=2: L=4: 1.6492189406, 8: 3.8393150668, 12: 5.8948058277
- N=2,k=3: L=8: 4.9267655562, 12: 7.5551867405
- N=2,k=4: L=8: 5.7817491256, 10: 7.3280616238
- N=3,k=1: L=6: 0.6972243623, 9: 1.2102821205, 12: 1.6877365977
- N=3,k=2: L=6: 1.2020467253, 9: 2.0739788158
- N=3,k=3: L=6: 1.5976165898
- N=4,k=1: L=8: 0.53949513, 12: 0.9445081373
- N=4,k=2: L=8: 0.9558401418

Report the ideal-string BA value next to each target with the difference. Ideal strings are not exact at finite L; quantify how the deviation scales with L (for N=2,k=2 go to the largest L where you can still compare, and for k=1 the equations are exact so agreement must be at round-off). If the reduced equations fail for k ≥ 2, say exactly how (e.g. no real solution, or a deviation that is O(1/L²) rather than exponentially small), and stop there rather than forcing agreement.

## Production (only if validation succeeds)

For (N,k) in {(2,2), (3,2), (4,2), (7,2), (8,2)} compute E_0(L) for L = multiples of N up to at least 50N (more if cheap), then A(L) = −6L[E_0(L) − e_0 L]/π using the e_0 formula above, and fit A(L) = A_∞ + a/L² + b/L⁴ (also try including 1/(L² ln L)-type terms if the k ≥ 2 data show slow drift; SU(N)_k at the integrable point can carry a marginal operator). Report A_∞ with an error estimate from varying the fit window, the predicted A = k(N²−1)π/(N(N+k)), and the ratio A_∞/(2π‖h‖) with ‖h‖ = H_k/2. Also check that the BA e_0 extracted from a free fit agrees with the formula. Output a JSON file with all E_0(L).
