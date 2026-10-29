# Structural lemmas: what the easy mechanisms cannot do

Notation as in `README.md`. h is a Hermitian bond on C^d (x) C^d, M(h) = min_o (1/2) spread(h + o (x) 1 - 1 (x) o), A the Casimir coefficient, R(h) = A/(pi M log2 d). The spread is convex in o (it is the difference of the largest and smallest eigenvalue, a convex minus a concave function of the matrix), so every local minimum of o -> spread(h + o (x) 1 - 1 (x) o) is global. All four lemmas are proved; the numerical checks are in `lemmas_check.py` and `data/lemmas_check.json`.

## Lemma 1 (stacking is closed)

Let h_1 on C^{d_1} and h_2 on C^{d_2} satisfy A_i <= pi M_i log2 d_i, and let h = h_1 (x) 1 + 1 (x) h_2 on C^{d_1 d_2} (site = (factor 1, factor 2)). Then A = A_1 + A_2 (the ring Hamiltonians commute and act on different factors, so ground energies add at every L) and

    M(h) >= max(M_1, M_2),

hence A <= pi max(M_1, M_2)(log2 d_1 + log2 d_2) <= pi M(h) log2 d. The stacked chain satisfies the bound with the strict slack pi [max(M_1,M_2) log2 d_1 + max(M_1,M_2) log2 d_2 - M_1 log2 d_1 - M_2 log2 d_2] >= 0 from the cross terms, and with equality only when one factor has M = 0.

Proof of M(h) >= max(M_1, M_2). Fix any o on C^{d_1 d_2} and a unit vector psi in C^{d_2}. Compress h + o (x) 1 - 1 (x) o with |psi> on the second factor of both sites. The compression of h_1 (x) 1 is h_1, that of 1 (x) h_2 is the constant <psi psi| h_2 |psi psi>, and that of o (x) 1 - 1 (x) o is o' (x) 1 - 1 (x) o' with o' = (1 (x) <psi|) o (1 (x) |psi>) a Hermitian d_1 x d_1 matrix. The spread of a compression never exceeds the spread of the matrix (Cauchy interlacing), so spread(h + o (x) 1 - 1 (x) o) >= spread(h_1 + o' (x) 1 - 1 (x) o') >= 2 M_1. Symmetrically for M_2. QED.

Remark. M(h) <= M_1 + M_2 (take o = o_1 (x) 1 + 1 (x) o_2), and the survey found cases with M(h) < M_1 + M_2 strictly. The lemma shows the deficit can never be large enough to matter: a stack of saturating chains has R = (M_1 log2 d_1 + M_2 log2 d_2)/(M log2 d) <= 1. The same argument shows a bound of the form A <= pi M f(d) is closed under stacking whenever f is superadditive, f(d_1 d_2) >= f(d_1) + f(d_2), which log2 d satisfies with equality and the harmonic candidate H_{d-1} of the README satisfies strictly (H_{ab-1} >= H_{a-1} + H_{b-1} for a, b >= 2). A d-independent bound A <= K M has no such closure and is in fact false (SU(N) rank-2 symmetric chains).

## Lemma 2 (the gauge is trivial for reflection-symmetric bonds with simple extremes)

Let h be invariant under the swap S (S h S = h) and let its largest and smallest eigenvalues be nondegenerate. Then M(h) = spread(h)/2, the gauge optimization does nothing.

Proof. Let |t>, |b> be the top and bottom eigenvectors. Since S h S = h and the eigenvalues are simple, S|t> = +-|t>, so the two one-site marginals of |t> coincide, rho_1(t) = rho_2(t), and likewise for |b>. The first-order change of the top eigenvalue under epsilon (o (x) 1 - 1 (x) o) is epsilon [tr rho_1(t) o - tr rho_2(t) o] = 0, and the second-order change is a sum of |matrix element|^2/(lambda_top - lambda_n) >= 0. The bottom eigenvalue shifts by 0 at first order and by a nonpositive amount at second order. So o = 0 is a local minimum of the spread, hence global by convexity. QED.

Remark. For degenerate extremes the argument needs the first-order splitting. If the extremal eigenspaces are the symmetric and antisymmetric subspaces (any bond of the form a 1 + b S has this property, in particular every SU(N) permutation chain and the Heisenberg chain), the operator (o (x) 1 - 1 (x) o) restricted to the symmetric subspace is zero (it is odd under S), so again no first-order gain and M = spread/2. For bonds without reflection symmetry the gauge does help (control rows in `data/lemmas_check.json`, where M is as much as 40 percent below the half spread).

## Lemma 3 (twists of a symmetric bond)

Let w be a unitary on C^d and h_w = (1 (x) w) h (1 (x) w^dag). On an open chain the Hamiltonian of h_w is unitarily equivalent to that of h through U = (x)_i w^i whenever h commutes with w (x) w (in particular for any w if h is the Heisenberg or an SU(N) permutation bond, and for w in a U(1) symmetry of h in general). On a ring the equivalence holds up to the boundary twist w^L, so the bulk CFT data (c_a, v_a) and hence A coincide, while

    M(h_w) = min_o (1/2) spread(h + o (x) 1 - 1 (x) w o w^dag),

a different gauge set. If M(h_w) < M(h) for some w, the twisted chain would beat the untwisted one in R at the same A. For the Heisenberg bond this never happens: with X = o (x) 1 - 1 (x) o', o' = w o w^dag, write X = [(o - o')/2 (x) 1 + 1 (x) (o - o')/2] + [(o + o')/2 (x) 1 - 1 (x) (o + o')/2]. The second bracket is odd under S and vanishes on the triplet (top) and singlet (bottom) subspaces; the first is a symmetric one-site operator a (x) 1 + 1 (x) a with a = (o - o')/2 traceless, which splits the triplet at first order (it acts as a . S_tot, with eigenvalues +-|a|/2 and 0) unless a = 0, and shifts the singlet by zero. Either way the top eigenvalue does not decrease at first order, second order pushes both extremes outward, and convexity finishes the argument: M(h_w) = 1/2 for every w. The same reasoning applies to every SU(N) permutation bond (symmetric top, antisymmetric bottom). So the whole family of "twisted Heisenberg" bonds -(1 (x) w)|s><s|(1 (x) w^dag), w in SU(2), none of which is SU(2) invariant as a bond, saturates the conjecture: equality in A <= pi M log2 d holds on this three-parameter family at d = 2, not only at the Heisenberg point. Numerics: `lemmas_check.py` finds M_w = 0.5000 for eight Haar-random w and ring energies identical to Heisenberg for w of order 2.

## Lemma 4 (blocking loses a factor b^2)

Block b sites into one site of dimension d^b. The blocked bond is h' = h_{b,b+1} + (1/2) sum (internal bonds of the two blocks). Velocities in blocked units are v/b and c is unchanged, so A' = A/b. The blocked bond contains, after compressing all sites but the two at the junction onto product states, a copy of h + c (x) 1 + 1 (x) c' with some one-site operators c, c' (from the compressed internal bonds and the compressed gauge), so M' >= M_2(h) := min_{c,c'} (1/2) spread(h + c (x) 1 + 1 (x) c'). Hence

    R(h') = R(h) . M/(b^2 M') <= R(h) . M/(b^2 M_2(h)).

For the Heisenberg bond M_2 = M = 1/2 (the first-order argument of Lemma 3 with independent c, c'), and the b = 2 blocked bond has M' = 3/4 numerically, so R' = 1/6, which the ring fit reproduces (0.1676 from blocks of 4 to 7). Encodings that spread a chain over several sites of a larger alphabet are therefore the wrong direction; the ratio is largest at the smallest alphabet.

## Consequence for the search

Lemmas 1 to 4 kill the "many commuting fast modes" mechanism (i) and the encoding mechanisms in their simple forms: stacking, twisting and blocking never raise R above the maximum of the constituents. Anything that beats the bound must do so with a genuinely interacting bond on a small alphabet, which is what the numerical search targets.
