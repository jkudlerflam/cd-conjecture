# Literature check for the Casimir bound and for qubit chains with c above 1

## Claim

A literature check plus one computation, for critical, translation-invariant, NN, Hermitian chains on (C^d)^⊗L flowing to unitary CFTs, c the total Virasoro central charge. (i) Nothing in print bounds Σ c_a v_a, c v, or the Casimir coefficient by a local norm such as ‖h‖, and nobody states that c is paid for with velocity or bandwidth. (ii) There is no classification or theorem giving c ≤ 1 for NN qubit chains, only the c–d conjecture itself, and I found no NN qubit chain with c > 1. (iii) The integrable SU(N) chains with the rank-k symmetric representation violate Σ c v ≤ 2π‖h‖, and the central charge and velocity from the Astra run are correct.

## Status

Survey and check, finished 2026-10-27. Items (i) and (ii) are negative literature results. Item (iii) is numerical plus a short Bethe-ansatz argument of my own, and a published velocity for the Andrei–Johannesson chains would upgrade it.

## How it was established and checked

- Author: Claude (Anthropic, Sonnet 5.5), literature agent.
- AI assistance: all searching, reading and numerics. No Codex run.
- Checked by: nobody yet. [read] means I opened the paper text. [abstract] means only the abstract page or a search summary. Pre-arXiv papers were not opened.

## Rank-k symmetric SU(N) chains

The site carries Sym^k of C^N, so d = C(N+k−1, k). The fused bond has eigenvalue 2(H_{k−j} − H_k) on the irrep with j boxes in the second row, H_m the harmonic number. The spread is 2H_k, and averaging over SU(N) shows no single-site shift lowers it, so ‖h‖ = H_k.

Central charge. Führinger et al., 0806.2563 [read], state that these chains are SU(N)_k with c = k(N²−1)/(N+k), cite Takhtajan, Babujian, Andrei and Johannesson for the Hamiltonians, and confirm c by DMRG for (N,k) = (2,2), (2,3), (3,2). Their Eqs. (11) and (12) give H^{[2,3]} and H^{[3,2]}. Appendix A of 1509.01680 [read, through a fetch summary] gives the k = 2 chain for all N as Σ S·S − N/(3N−4)(S·S)² with the same c. I checked that all three have bond eigenvalues proportional to 2(H_{k−j} − H_k), with scales 8, 6/5 and 2N/(3N−4). So c and the spectrum are confirmed. The original papers (Kulish–Reshetikhin 1983, Babujian 1983, Andrei–Johannesson, Phys. Lett. A 104, 370, 1984, Johannesson, Nucl. Phys. B 270, 235, 1986, Affleck, Nucl. Phys. B 305, 582, 1988) are known to me only from citations.

Velocity. I found no arXiv paper that quotes v for general (N,k). For N = 2, k = 2, Vlijm and Caux, 1401.4450 [read in part], use H = (J/4)Σ[S·S − (S·S)²] with spinon dispersion (π/2)|sin q|. That is v = 2π for the eigenvalues (−6,−2,0), hence π in the normalization above, matching 2π/N. For general (N,k) my own argument uses the standard string ground state, with k-strings at every nesting level. Every kernel and the driving term equal the k = 1 form times one common function, (1 − e^{−k|ω|})/(2 sinh(|ω|/2)). Dividing it out leaves k-independent density equations 2cosh(ω/2) ρ̂_a − ρ̂_{a−1} − ρ̂_{a+1} = δ_{a1}, so ρ̂_1 = sinh((N−1)ω/2)/sinh(Nω/2). The spinon energy and momentum follow from ρ_1, whose tails decay as e^{−2π|x|/N}, so v = 2π/N for every k, and e_0 = −2∫_0^∞ dω ρ̂_1(ω)(1 − e^{−kω})/(2 sinh(ω/2)).

Exact diagonalization on rings with kL divisible by N, using this e_0, gives v_eff from E_0(L) − L e_0 = −π c v_eff/(6L). For (2,2) at L = 14 it is 3.201, for (2,3) at L = 10 it is 3.261, and for (3,2) at L = 3, 6, 9 it is 2.718, 2.263, 2.180, against π = 3.142 and 2π/3 = 2.094. A 1/L^{3/2} correction takes the first two to 3.149 and the (3,2) points at L = 6, 9 to 2.08. I skipped (4,2), since gcd(N,k) = 2 and the reachable rings did not show a vacuum ground state.

Consequence. Σ c v/(2π‖h‖) = k(N²−1)/((N+k) N H_k), which is 1 − 1/N at k = 1. Over all k it stays below 0.51, 0.74 and 0.93 for N = 2, 3, 4. It first exceeds 1 at (N,k) = (7,2), d = 28, value 1.016, then (6,3), d = 56, value 1.061, and (5,4), d = 70, value 1.024. At fixed k and large N it tends to k/H_k, which is 4/3 for k = 2 and unbounded in k, so no constant multiple of ‖h‖ bounds Σ c v on this family. [own] A bound with a log d factor is not excluded. The c–d bound itself holds easily.

## Question 1, bounds on c v

No paper I found does this. The c–d paper, 2403.17242 [read], never mentions velocity, bandwidth or a norm. 

Energy current and sum rules. Bernard and Doyon, 1202.0239 [abstract], give the CFT energy current (πc/12)(T_L² − T_R²), and 2507.23439 [abstract] tests it on spin chains. With a lattice bound on the current this limits c T², the repository's item (b), not a published bound. Kapustin and Spodyneiko, 1904.05491 [read], prove equilibrium energy currents vanish for finite-range 1D lattices, so c_L = c_R. Shastry, cond-mat/0508711 [abstract], gives a sum rule for the dynamical thermal conductivity as a local expectation value, with no bound on c or v. Cond-mat/0211633 [abstract] gives a thermal Drude weight proportional to c v T for conformal modes. Charge Drude weights are bounded by the optical f-sum (Kohn, Shastry–Sutherland, recalled), with upper bounds in 2212.08711 [abstract]. I found no extension to non-abelian levels, which the repository's argument 1(d) needs.

Velocity bounds. The review 2206.15126 [abstract] calls identifying the Lieb–Robinson velocity with a maximal quasiparticle velocity a heuristic, and 1904.12464 [abstract] bounds relative group velocities for free-fermion bands. No paper, including the review 2303.07386 [abstract], states v ≤ O(‖h‖) for an emergent CFT velocity as a theorem.

Stress-tensor sum rules. Cappelli, Friedan and Latorre (Nucl. Phys. B 352, 616, 1991) [abstract] and Hartman and Mathys, 2405.10137 [abstract], give sum rules for the change of c along continuum flows. They use Lorentz invariance and have no lattice norm version. Koo–Saleur generators (hep-th/9312156, 1706.01436, 2107.13834 [abstract]) obtain c from a commutator of the energy density with the energy current, but [own] that commutator carries non-universal analytic terms, so it yields no inequality.

Casimir term. E_0 = e_0 L − πcv/(6L) is due to Blöte–Cardy–Nightingale and Affleck (1986, recalled). I found no rigorous inequality on the 1/L coefficient. Nobody states that c is paid for with velocity. The nearest remark is Pendry's per-channel limit on energy flow (quant-ph/0001007 [abstract]).

## Question 2, qubit chains

No classification or theorem exists. Latorre and Sierra, 2403.17242 [read], write that for d = 2 there is "no other model made with spins" with larger c than the Heisenberg chain, c = 1. Li, Lin and McGreevy, 2509.04596 [abstract], expect no CFT with c > d − 1 and find on four qubits only Ising, XX and Heisenberg. Alavirad and Barkeshli, 1910.00589 [abstract], give an LSM bound c ≥ n − 1 for Z_n × Z_n chains, the XYZ case at n = 2, so qubits sit at c = 1. Free-fermion NN qubit chains have c ≤ 1 by the zero count of Verresen, Jones and Pollmann, 1709.03508 [read].

Interacting Kitaev chain. Chepiga and Mila, 2206.11754 [read], study hopping, pairing and NN repulsion (their Eq. 1), the NN XYZ chain in a longitudinal field. They find Ising lines, floating phases with emergent U(1) and a c = 1 eight-vertex point, and no larger c. The c = 3/2 Ising plus Luttinger phase of Rahmani, Zhu, Franz and Affleck, 1504.05192 [abstract] and 1505.03966 [read], has an extra σ^xσ^x at distance two in spin form (their Eq. 1.7), so it is not NN. Dutta and Sen, cond-mat/0208216 [abstract], find the XYZ chain in a field gapless only on a line, with no c stated.

Other models. Uniform Dzyaloshinskii–Moriya terms gauge away to XXZ with a twist, so c = 1 (1902.09356, 1403.4568 [abstract]). Kitaev–Γ has an emergent SU(2)_1 phase, c = 1 (1910.14304 [abstract]), and Kitaev–Heisenberg–Γ has four Luttinger phases (2004.12954 [abstract]). No c above 1 is stated in 2509.21901 [abstract].

Two fluids. Gotta, Mazza, Simon and Roux, 2012.12393 [read in part], find a c = 2 coexistence of paired and unpaired fluids with NN hopping plus a pair hopping c†_{j+1}c†_j c_j c_{j−1}, which spans three sites. Ruhman and Altman, 1705.00635, and 2006.07128 [read in part] find c = 3/2 at pairing transitions with three-site or NNN terms. All have range three or more.

## Reproducing

`rank_k_ed.py N k L1 L2 ...` prints E_0, v_eff and the prediction. It takes seconds up to 10^6 states and about 28 minutes for (3,2) at L = 9. `rank_k_ratio.py` prints the ratio table and first violations. Both need only numpy and scipy.

## Relation to the rest of the repository

This tests claim 3 of `contrib/refined_conjecture`. The rank-k family refutes Σ c v ≤ 2π‖h‖ from (7,2), and "Where things stand" should say so.
