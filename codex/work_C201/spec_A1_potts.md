# Task: spectral inclusion of the periodic Q-state Potts quantum chain in twisted XXZ rings, with a rigorous proof

You are asked for a complete, rigorous proof (every step shown), written for a mathematical-physics referee. Tag every step [rigorous], [rigorous given cited theorem X] (state X precisely with reference), or [not proved]. Do not hand-wave. If a step cannot be proved, say so and isolate the missing statement. Plain declarative prose, no em-dashes.

## Definitions

Fix an integer Q >= 2 and L >= 2. Potts ring: Hilbert space (C^Q)^{⊗L}, sites j in Z_L, basis |s>, s in Z_Q. Define on-site and bond operators
  E_j = Q^{-1/2} sum_{a,b in Z_Q} |a><b|  on site j  (so E_j = sqrt(Q) times the projector on the uniform vector),
  F_{j,j+1} = sqrt(Q) sum_{s in Z_Q} |s s><s s|  on sites j, j+1 (indices mod L).
Set N = 2L and e_{2j-1} = E_j, e_{2j} = F_{j,j+1} for j = 1..L, with indices of e mod N. They satisfy the periodic Temperley-Lieb relations with loop weight delta = sqrt(Q): e_i^2 = delta e_i, e_i e_{i±1} e_i = e_i, [e_i, e_k] = 0 for |i-k| >= 2 cyclically (verify this). Each e_i is Hermitian. Define
  H_P^{(s)}(L) = -s sum_{i=1}^{N} e_i,   s = +1 (ferromagnetic critical Potts chain) or s = -1.
As a translation-invariant nearest-neighbour chain it has bond h = -s [ F_{j,j+1} + (E_j + E_{j+1})/2 ].

XXZ ring with twist: N spins 1/2, Pauli matrices sigma. For real gamma and real phi define
  H_X(N; Delta, phi) = (1/2) sum_{i=1}^{N-1} [ sigma^x_i sigma^x_{i+1} + sigma^y_i sigma^y_{i+1} + Delta (sigma^z_i sigma^z_{i+1} - 1) ]
                       + (1/2) [ 2 e^{i phi} sigma^+_N sigma^-_1 + 2 e^{-i phi} sigma^-_N sigma^+_1 + Delta (sigma^z_N sigma^z_1 - 1) ].
It commutes with S^z = (1/2) sum sigma^z. Write Spec(N; Delta, phi, m) for its spectrum in the sector S^z = m.

## Claims to prove (or correct)

Let delta = sqrt(Q) = 2 cos(gamma) for Q <= 4 (gamma = pi/4, pi/6, 0 for Q = 2, 3, 4).

Claim 1 (inclusion). For Q in {2, 3, 4}, every eigenvalue of H_P^{(+1)}(L) = -sum e_i lies in the union over integers m with |m| <= L and over real phi of Spec(2L; Delta = cos(gamma), phi, m). Every eigenvalue of H_P^{(-1)}(L) = +sum e_i lies in the union over m and real phi of Spec(2L; Delta = -cos(gamma), phi, m) (use the rotation by pi about z on odd sites of the XXZ ring to relate + and - signs and track the twist, N = 2L is even).
It is essential that only REAL twists phi occur (unimodular affine Temperley-Lieb parameter z), because the application uses Hermitian twisted XXZ rings. Prove this point with care. Possible route: the periodic (affine) Temperley-Lieb algebra and its standard (cell) modules W_{j,z} of Graham and Lehrer and of Martin and Saleur; simple subquotients of any finite-dimensional module are quotients of standard modules; the spectrum of sum e_i on W_{j,z} equals the twisted XXZ sector spectrum (isomorphism for generic parameters plus continuity of characteristic polynomials); in a *-representation with Hermitian e_i and a unitary two-step translation, only |z| = 1 occurs. Another possible route is the explicit result of Aufgebauer and Klumper, J. Stat. Mech. (2010) P05018, arXiv:1003.1932, on periodic Temperley-Lieb chains. You may cite published theorems, stating them precisely, but must supply every step connecting them to Claim 1, including the roots-of-unity cases (all of Q = 2, 3, 4 are roots of unity or q = 1) and the precise dictionary between z, the number of through lines 2j, the sector m and the twist phi. If the correct statement needs a different twisted-XXZ normalization from the one above, give the corrected statement.

Claim 2 (ground state). Identify which sector and twist contain the ground state energy of H_P^{(+1)}(L) (expected: m = 0 and phi = ±2 gamma, i.e. non-contractible loops also have weight sqrt(Q)), and prove at least the inequality E_P^{(s)}(L) >= min over m and real phi of min Spec(2L; s cos(gamma), phi, m), plus an upper bound of the form E_P^{(s)}(L) <= min Spec(2L; s cos(gamma), phi_0, m_0) for one explicit (phi_0, m_0) (so that the bulk energy densities agree).

Claim 3 (consequence). Suppose the following is known for some Delta_0 in (-1, 1] (it is a cited theorem in our project; treat it as hypothesis (X)): there are constants J e_bulk(Delta_0) and K such that for every even N >= 2, every sector m and every real twist phi,
   min Spec(N; Delta_0, phi, m) >= N eps_0 - K / N,   where eps_0 = lim_{N->infinity} (1/N) min_{m} min Spec(N; Delta_0, 0, m).
Deduce: the Potts ring ground energy obeys E_P(L) >= L e_0^P - K/(2L) for all L >= 2, where e_0^P is the infinite-chain ground energy density per Potts site of the bond h, and prove e_0^P = 2 eps_0. Here e_0^P is defined as lim E_P(L)/L (prove the limit exists and equals the infimum of the energy density over translation-invariant states of the infinite chain; standard).

Claim 4 (the gauge-optimized half-spread). For a bond h on C^d ⊗ C^d define M(h) = min over Hermitian X on C^d of (lambda_max - lambda_min)/2 of h + X ⊗ 1 - 1 ⊗ X. Compute or bound from below M for the Potts bond h = -s[F + (E ⊗ 1 + 1 ⊗ E)/2] for Q = 3 and Q = 4, with proof. (A valid lower bound is M >= (e_max - e_min)/2 with e_max, e_min the maximal and minimal energy densities over translation-invariant infinite-chain states, each bounded by product states. An exact value is better.) Also give an upper bound.

## Output
A self-contained proof document with the tags above, a list of every cited theorem with precise statement and reference, explicit small checks you performed (for example L = 2, Q = 2 or 3 by hand or by exact computation: list the Potts eigenvalues and the matching (m, phi)), and a summary table of which claims are fully rigorous and which rest on citations.
