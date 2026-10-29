# C76 derivations and scope

All chains below have one fixed two-site bond. On-site terms are included by splitting them equally between the adjacent bonds. The scalar c denotes the sum of nonchiral sector central charges, so a compact boson contributes one even when its velocity differs from those of other sectors. Statements labeled “proved” are elementary deductions here, not new proofs of the quoted Bethe-ansatz continuum limits.

## 1. What q species and a hole actually count

Let the local orthonormal basis be |0>, |1>, ..., |q>, with X^{ab}=|a><b|, n_a=X^{aa}, and n=sum_{a=1}^q n_a. The hard-core constraint is

    n_0 + n_1 + ... + n_q = 1.

It follows as an operator identity that δn_0 = -sum_a δn_a. In particular, the (q+1)-by-(q+1) connected density covariance has the all-ones vector in its kernel. A hole density is not an independent conserved density. [proved]

Suppose the complete low-energy theory consists of q nondegenerate density liquids, with no extra neutral critical sector or extra Fermi pockets. Introduce q conjugate pairs φ_a, θ_a, normalized by

    δρ_a = -(1/π)∂x φ_a,
    [φ_a(x), ∂y θ_b(y)] = iπ δ_ab δ(x-y).

For parity-invariant couplings their quadratic Hamiltonian is

    H_LL = (1/2π) ∫dx [(∂xθ)^T V_J (∂xθ)
                       + (∂xφ)^T V_N (∂xφ)],

where V_J and V_N are positive definite. Set φ=V_J^{1/2} O Φ and θ=V_J^{-1/2} O Θ. This preserves the commutator. Choose O to diagonalize V_J^{1/2} V_N V_J^{1/2}, with eigenvalues v_a^2>0. The Hamiltonian becomes q independent scalar Hamiltonians, with θ stiffness one and φ stiffness v_a^2. Reciprocal rescalings of each canonical pair put both coefficients equal to v_a. Hence there are q oscillator towers, each with c=1. Unequal hopping changes the matrices and velocities, not their size. Mixed current-density terms can change left and right velocities but do not increase the number of canonical pairs. [proved within the complete Gaussian hypothesis]

At the SU(q)-symmetric point define the charge coordinate

    φ_c = (φ_1+...+φ_q)/sqrt(q),

and q-1 orthonormal coordinates perpendicular to (1,...,1), with the same transformation on θ. This is a change of basis in a q-dimensional space. Nonabelian bosonization identifies the relative coordinates with SU(q)_1, whose Sugawara central charge is

    c_flavor = dim SU(q)/(1+q) = (q^2-1)/(q+1) = q-1.

Thus c_total=1+(q-1)=q, not q+1 and not 2q-1. Left and right Fermi points are likewise the two chiral halves of one liquid, not two nonchiral liquids.

Independent U(1)^q conservation forbids vertices that change the individual species numbers. Density-locking vertices carry reciprocal-lattice mismatch determined by integer combinations of the densities. Jointly incommensurate densities exclude nonzero such combinations. Small symmetry-preserving changes of hopping and forward scattering then leave the q-mode Gaussian theory gapless, provided its stiffness remains positive. This is standard bosonization/RG input, not a proof of the Gaussian hypothesis for every microscopic bond.

The limitation matters. A rank-q charge algebra does not bound the central charge of a neutral coset. A neutral Ising or parafermion sector, several independently surviving pockets, or independently stable bound-complex fluids could invalidate the completeness hypothesis. They need a dynamical construction. Conservation of q charges alone neither supplies them nor rules them out. Consequently this derivation rules out the proposed extra mode from spin-charge separation, but is not a universal theorem c<=q for all color-conserving nearest-neighbor bonds.

## 2. Pure hopping freezes particle order

Consider

    h_hop = -sum_a t_a (X^{a0}⊗X^{0a}+X^{0a}⊗X^{a0}).

On an open chain, specify a basis state by occupied positions

    x_1 < x_2 < ... < x_N

and their ordered color word w=(a_1,...,a_N). Every allowed move exchanges a particle with a hole. It cannot interchange x_j and x_{j+1}. Therefore w is an exact conserved sector label, including for unequal t_a and arbitrary diagonal two-site interactions. On a ring only cyclic rotation of the word is possible. [proved]

For t_a=t and zero interaction, in every open-chain word sector the positional Hamiltonian is the same hard-core spinless hopping Hamiltonian. Its single-particle energies are

    ε_r=-2t cos(πr/(L+1)), r=1,...,L.

The many-body energies are sums over N distinct r. For fixed color numbers N_a every positional spectrum occurs N!/prod_a N_a! times. The propagating charge sector has c=1 at 0<N/L<1. The extensive frozen color degeneracy is not q-1 additional conformal modes. Its velocity is zero, so treating it as a CFT is invalid.

For L=8, N_1=N_2=2, the six words each have C(8,4)=70 positional states. Full ED gives six identical spectra, maximum discrepancy 1.155e-14. Setting t_2=0.7 changes the spectra but leaves exactly the same six connected components. These are finite-size algebraic checks, not an extrapolated c measurement.

Exchange terms X^{ab}⊗X^{ba} let adjacent particles pass and make flavor dynamics possible. This explains why the limit J→0 of a t-J liquid is singular. At each fixed J>0 a flavor velocity can be nonzero; setting J=0 before taking the infrared limit loses that sector. A small velocity does not increase c through a divergent 1/v thermal coefficient.

## 3. Explicit saturating hard-core bonds

Define

    P_occ = sum_{a,b=1}^q X^{ab}⊗X^{ba},
    T = sum_{a=1}^q (X^{a0}⊗X^{0a}+X^{0a}⊗X^{a0}).

The generalized supersymmetric t-J bond is

    h_tJ = -t T + t(P_occ - n⊗n).

For q=2 this is -t hopping + 2t(S_i·S_j-n_i n_j/4). If p_0=0 and p_a=1 for a>0, the graded permutation is

    P_g|ab> = (-1)^{p_a p_b}|ba>,
    -tP_g = h_tJ + t(n⊗1+1⊗n-1).

The identity follows separately on the hole-hole, mixed, and occupied-occupied blocks. Thus h_tJ differs from a rational graded-permutation bond only by a chemical potential and a constant. The matrix is Hermitian in the ordinary positive definite site inner product. The finite-density antiferromagnetic phase has one charge liquid and SU(q)_1 flavor modes, giving c=q at d=q+1. This continuum identification is literature input; the bond identity is proved here.

A second saturator avoids grading entirely:

    P_{q+1}=sum_{a,b=0}^q X^{ab}⊗X^{ba}.

Its antiferromagnetic chain is SU(q+1)_1. A staggered on-site hole/particle phase changes the sign of the mixed hopping on even rings, giving the uniform bond

    P' = -T + P_occ + X^{00}⊗X^{00}
       = h_tJ/t + 2n⊗n - n⊗1 - 1⊗n + 1.

All coefficients are order one. Generic small Cartan fields select interior, potentially incommensurate color densities while preserving the multicomponent liquid until a sea empties. These chains establish attainability of c=q, not an upper bound on all hard-core models.

For a mixed grading with at least one even color, the bare bond -P_g has lowest eigenvalue -1, and a uniform even-color product state attains -1 on every bond. Thus the zero-field chain has ground energy -L. A one-particle color excitation has discrete-Laplacian dispersion 2(1-cos k), hence z=2 at its threshold. Finite-density conformal t-J phases require the chemical-potential term implicit in the preceding identity. For +P_g, a uniform odd-color state gives the analogous conclusion. This proves why an unqualified central-charge assignment to every bare SU(N|M) permutation chain is wrong. [proved]

More generally, if all q species numbers are conserved, a two-site bond can only mix |ab> with |ba> for a≠b and act diagonally on |aa>. This follows because the unordered multiset of the two site labels is fixed by the charges. Its general form is a diagonal potential plus complex exchange amplitudes. Even this restricted finite-dimensional bond family is not proved here to obey c<=q outside its complete Gaussian regimes.

## 4. Interactions can exceed the free-fermion bound

C1 supplies c<=m for m canonical complex orbitals per site, d=2^m, when the Hamiltonian is quadratic and its zeros are isolated and linear. It does not apply to a general two-site operator on the same Fock space.

Label the occupation basis of m orbitals by a binary string A. The explicit interacting bond

    h_swap = sum_{A,B∈{0,1}^m} |B,A><A,B|

is the ordinary SU(2^m) permutation. It conserves each global orbital occupation, indeed each global local-configuration count. It is parity-even, Hermitian, and supported on two adjacent sites. Its continuum limit is SU(2^m)_1, so

    c=2^m-1.

For m=2, d=4, c=3>2. For m=3, d=8, c=7>3. This is a strong interacting realization, not a weak-coupling correction to free fermions. It contains correlated hopping and configuration exchange and is not being represented as a pure density-density Hubbard model.

For d=4 an explicit canonical-fermion expression can be given on an oriented bond with modes ordered a_i,b_i,a_j,b_j. Set

    D_a=1-n_ai-n_aj+2n_ai n_aj,
    D_b=1-n_bi-n_bj+2n_bi n_bj,
    Q_a=D_a+[a_i†(1-2n_bi)a_j+h.c.],
    Q_b=D_b+[b_i†(1-2n_aj)b_j+h.c.].

These are the swaps of the two occupation qubits for orbital a and orbital b. They commute and Q_a Q_b=P_4. The parity factors are the short Jordan-Wigner strings between the corresponding modes; dropping them changes the model. The formula is a finite fermionic polynomial, including higher than quartic terms. At a periodic seam the tensor-product prescription corresponds to the usual parity-dependent fermionic twist, which affects finite-size sectors, not the bulk c count. The simpler spin-orbital notation is

    P_4=(2S_i·S_j+1/2)(2T_i·T_j+1/2).

At weak coupling around m ordinary bands with one filled interval per band, linearization produces m right-left pairs. Forward scattering acts on their m canonical bosons. Relevant backscattering can remove modes; additional modes require a change of the low-energy degrees of freedom, not a re-labeling into charge and spin. This is a local weak-coupling statement, not a general interacting extension of C1.

## 5. Orthosymplectic Hermiticity audit

The rational vector OSp(n|2m) bond in the convention of Frahm and Martins is proportional to

    h=P_g + 2 E/(2-z), z=n-2m.

Choose the invariant matrix A=diag(I_n,J), where J is a direct sum of m real 2-by-2 antisymmetric blocks with J_12=1. In the ordinary tensor basis

    E_{ab,cd}=A_ab (A^{-1})_cd.

Let u=vec(I_n) and v=vec(J), embedded in the boson-boson and fermion-fermion sectors. Then

    E=(u+v)(u-v)^T,
    E^2=zE,
    E-E†=2(vu^T-uv^T).

The supports are orthogonal, ||u||²=n, ||v||²=2m. Therefore

    ||E-E†||_F²=16nm,
    ||h-h†||_F=8 sqrt(nm)/|2-z|.

For n,m>0 and z≠2 this is nonzero. It is not a boundary derivative: its matrix elements convert a boson-boson pair into a fermion-fermion pair with the opposite reverse sign. Such off-diagonal two-site transitions cannot come from an on-site telescoping term. Consequently these standard bonds fail the stipulated Hermiticity condition. Changing to a nonlocal metric does not exhibit a Hermitian nearest-neighbor bond on the specified tensor-product Hilbert space. The z=2 singular normalization requires a separate limit; its surviving E term still has the same Hermiticity obstruction for mixed grading.

This does not say that every Hamiltonian mentioned in connection with OSp is forbidden. It audits this rational family. Pure orthogonal and symplectic chains are different Hermitian cases; their antiferromagnetic signs and continuum sectors must be treated separately. In particular one must not apply a mixed-grading central-charge formula mechanically at n=0.

## 6. Finite-size estimators and their limitations

For several independent velocities,

    E_0(L)=e∞L-(π/6L)sum_a c_a v_a+o(L^-1),
    f(T)=e∞-(πT²/6)sum_a c_a/v_a+o(T²).

These determine different velocity-weighted sums. Dividing either by an arbitrary common velocity is not a c count. For a translation-invariant nondegenerate conformal ground state,

    S(ell)=(c/3)log[(L/π)sin(πell/L)]+s_0+corrections.

The half-chain size difference cancels s_0 at fixed density and compatible ground-state sectors:

    c_half(L,2L)=3[S_{2L}(L)-S_L(L/2)]/log 2.

For the SU(d) swap normalization used here, the Bethe-ansatz values used for calibration are

    e∞=1+(2/d)[ψ(1/d)+γ_E], v=2π/d.

These are reference inputs, not fitted to the ED data. They give the single-size estimator c_Cas=-6L(E_0-e∞L)/(πv). For SU(3) it decreases from 2.11346 at L=6 to 2.03514 at L=12. The half-chain estimate from L=6,12 is 1.98240. Both are consistent with c=2. They do not prove a thermodynamic limit by themselves.

Schmidt decompositions must use only charges actually conserved by the state. Swap and t-J sectors have fixed counts of every local label, allowing those charge blocks. Hubbard conserves N_up,N_down but not the number of doubly occupied sites. Its Schmidt blocks in the code are therefore indexed by the two orbital charges, not the four configuration counts. Free-Hubbard entropy at the closed-shell L=8 point is independently compared to the correlation-matrix answer in the audit script. Degenerate ground-state vectors and different shell sequences are excluded from c extrapolations.
