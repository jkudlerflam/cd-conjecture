# Second-order stability of the joint Majorana code

2026-10-27, resumed C8. The finite-dimensional construction and all effective-Hamiltonian coefficients below are proved. Their continuum interpretation uses the stated free-fermion operator expansion and perturbative RG. No finite-coupling microscopic CFT or gap is proved.

## 1. A graded joint code with controlled second moments

Split a q-dimensional port into grading eigenspaces of dimensions q_+,q_-. Choose m=2n<=2 q_+ q_- traceless Hermitian matrices F_a odd under the grading, normalized by tr(F_a F_b)=q delta_ab. Normalized real and imaginary matrix units between the two grading sectors suffice. Let Gamma_a,R_a=i Gamma_a P_log be the logical operators from RESULT.md. Use the joint completely positive map of that report, with only these odd F_a and a sufficiently small alpha, for example alpha=1/(4m sqrt(q)). Its Stinespring isometry V satisfies

    V^dagger F_La V=alpha Gamma_a,
    V^dagger F_Ra V=alpha R_a,
    V^dagger F_La F_Lb V=V^dagger F_Ra F_Rb V=delta_ab I,
    V^dagger F_La F_Rb V=0.                                    (1)

The second line follows because F_a F_b is grading even. Its only surviving term under the channel is its normalized trace. The last line follows directly from the joint channel, which has no cross-port terms.

This code can preserve an exact total parity. In the positive-square-root Choi dilation, write z for the diagonal port grading. The Choi matrix commutes with z_L z_R tensor P_log. The physical block grading

    Omega = (z_L tensor z_R) tensor (z_L tensor z_R tensor P_log)

then obeys Omega V=V P_log, and both physical port families are odd under Omega. Here the second parenthesis acts on the Choi environment. Thus an Omega-invariant block parent and the port bond sum_a F_Ra tensor F_La preserve product_j Omega_j. This construction does not require a large continuous microscopic symmetry.

For example, q=16 admits m=128 odd channels, with intended c=64. It exceeds d-1 for a d=4q selector by one unit, and exceeds d-1 for the d=3q+1 selector by sixteen units. These are resource comparisons, not microscopic Hamiltonians.

## 2. The flat parent, including adjacent bonds

Let P=VV^dagger, Q=1-P, and take an unrestricted block parent H_0=Delta Q. Couple identical blocks by lambda v, where

    v_(j,j+1)=sum_a F_Ra,j F_La,j+1,
    T_j=sum_a R_a,j Gamma_a,j+1.

The first-order target is lambda alpha^2 sum_j T_j. A finite-volume degenerate perturbation expansion on the product code gives

    H_eff^(2)=(lambda^2/Delta) {
       -M m(1/2+alpha^2) I
       +alpha^4 [(3/2)sum_j T_j^2
                          +sum_j {T_j,T_(j+1)}] }.              (2)

The displayed scalar is for a ring of M>=4 blocks. For an open chain it is proportional to the number of bonds, with only existing adjacent pairs included.

Proof for one bond. On two blocks the inverse excitation energy, after removing Delta, is

    R=P tensor Q + Q tensor P +(1/2) Q tensor Q.

Equation (1) gives Pv^2P=mI, while the intermediate projection onto both codes gives alpha^4 T^2. Intermediate states with only one block excited contribute, in total,

    2m alpha^2 I - 2alpha^4 T^2.

Combining the one- and two-excitation denominators, with the perturbative minus sign, yields

    -PvRvP=-m(1/2+alpha^2)I+(3/2)alpha^4 T^2.                   (3)

For adjacent bonds, the two outer blocks must remain in their codes, so only the common block can be excited. Its relevant matrix element is

    V^dagger F_La Q F_Rb V=-alpha^2 Gamma_a R_b.                 (4)

The two outer matrix elements supply alpha^2, the denominator is Delta, and the minus sign in the second-order formula cancels the minus sign in (4). The two orderings therefore give alpha^4 {T_j,T_(j+1)}/Delta. Disjoint bonds cannot return to the product code after creating an excitation, so they give zero. This proves (2).

Keeping only (3) misses the adjacent-bond term and gives the wrong sign in the staggered mass channel.

## 3. Continuum interpretation and its limitation

At the symmetric hopping point, each Majorana lattice field has linear modes at momenta zero and pi. Use the standard free-field expansion

    gamma_(j,a) ~ psi_Ra(x)+(-1)^j psi_La(x).

The alternating part of the bond T_j is a nonzero real multiple of

    (-1)^j Epsilon(x),    Epsilon=sum_a i psi_Ra psi_La.

Therefore a lattice correction

    a sum_j T_j^2 + b sum_j {T_j,T_(j+1)}

has a staggered mass-channel coefficient proportional to a-2b. Contact contractions alter scalar and lower-dimensional terms but do not change the coefficient of the cross-species quartic operator in this expansion.

For (2), a-2b=-alpha^4/2. Thus the coefficient of Epsilon^2 is negative. This is the attractive sign of the SO(m) Gross-Neveu interaction. One can identify the one-loop tendency without assuming a lattice mass theorem. Put epsilon_a=i psi_Ra psi_La and O=sum_(a<b) epsilon_a epsilon_b. The contraction of a common species in the product of two summands gives a positive O/|z|^2 term. For a given surviving pair, there are m-2 choices of contracted species. Hence the self-OPE coefficient is positive for m>2. With Euclidean action S_0+g integral O, integrating out a short-distance shell gives dg/dell=-C(m-2)g^2+O(g^3), C>0. Negative g runs to larger negative magnitude, while positive g runs toward zero at this order.

For m=2 the factor m-2 vanishes and the interaction is the usual Gaussian marginal deformation. For m>2 the flat parent is therefore not justified as a stable realization of the symmetric free target. Standard Gross-Neveu field theory associates the attractive flow with mass generation and spontaneous dimerization. In this report that mass-generation identification is assumed continuum input, and the conclusion about the actual microscopic chain remains heuristic. Higher orders, symmetry breaking, and shifted Fermi points have not been controlled.

## 4. Two-energy spectral weights reverse the coefficient

The negative sign in Section 3 is not forced by the joint code. It can be changed explicitly while preserving (1), the code gap, and total parity.

Let S_s=Q F_s V, with s running through the m left and m right channels. Write

    U_L=[Gamma_1 ... Gamma_m],    U_R=[R_1 ... R_m],
    U=[U_L U_R].

The syndrome Gram matrix is

    G=S^dagger S=I-alpha^2 U^dagger U.                          (5)

Since U U^dagger=2m I, this is positive when 2m alpha^2<1. Prescribe excited spectral weights on the syndrome source as

    M_e = p_e I + [ s_e U_L^dagger U_L   t_e U_L^dagger U_R
                    t_e U_R^dagger U_L   s_e U_R^dagger U_R ].   (6)

Their sum equals G provided

    sum_e p_e=1,    sum_e s_e=sum_e t_e=-alpha^2.                (7)

Positivity is transparent. On the dark subspace the eigenvalue is p_e. On the two bright copies of the D-dimensional logical space the eigenvalues are

    p_e+m(s_e+t_e),    p_e+m(s_e-t_e).                          (8)

Thus it suffices that p_e>=0 and p_e+m s_e>=m|t_e|.

Use two excitation energies and weights

    E_1=Delta,           E_2=3Delta,
    p_1=p_2=1/2,
    s_1=s_2=-alpha^2/2,
    t_1=0,              t_2=-alpha^2.                          (9)

The smallest weight eigenvalue is 1/2-(3/2)m alpha^2, which is strictly positive for the alpha chosen above. Their sum is (5).

For completeness these are spectral weights of an actual finite matrix parent, not merely formal coefficients. Form the positive normalized measure

    Mtilde_e=G^(-1/2) M_e G^(-1/2),    sum_e Mtilde_e=I.

Its square-root dilation W=stack_e sqrt(M_e) G^(-1/2) is an isometry. The physical frame S G^(-1/2) is also an isometry. Extend both to orthonormal frames, map one frame to the other, and assign energies E_e to the corresponding coordinate blocks. Padding states can have energy 3Delta. This is the same finite-dimensional spectral dilation as RESULT.md Section 4, now with a nonidentity syndrome Gram matrix. It leaves ker H_0=ran V and gap Delta. The dilation needs at most 4mD excited dimensions, and the original q^4 D Choi block has enough room for the stated odd-channel family.

Every M_e commutes with I_(2m) tensor P_log because its nonscalar entries are products of two odd logical operators. The dilation and its orthonormal completions can consequently be chosen separately in the two parity sectors. The resulting H_0 commutes with Omega. The supplied Python implementation performs these separate completions.

To calculate the effective coefficients, remove Delta from the energies and define

    C_ef=1/(E_e+E_f),
    p_inv=sum_e p_e/E_e,    s_inv=sum_e s_e/E_e,
    t_inv=sum_e t_e/E_e.

For general weights (6), the same-bond correction is a scalar plus A T_j^2 and the adjacent correction is B {T_j,T_(j+1)}, with

    A=-2alpha^2 s_inv-s^T C s,
    B=-alpha^2 t_inv,                                         (10)

and the scalar per bond is

    -m[p^T C p+2p^T C s+2alpha^2 p_inv].                       (11)

The first term in A comes from single-block virtual excitations. The second comes from two simultaneously excited blocks. Formula B is the cross-port resolvent on their common block. These contributions are obtained by the same insertion of spectral projectors used in (3)-(4).

For (9), in units Delta=1,

    s_inv=-(2/3)alpha^2,    t_inv=-(1/3)alpha^2,
    s^T C s=(7/24)alpha^4.

Hence

    A=(25/24)alpha^4,    B=(1/3)alpha^4,
    A-2B=(3/8)alpha^4>0.                                     (12)

The scalar is -m[7/24+(3/4)alpha^2]. The leading hopping remains lambda alpha^2 sum T_j. The sign of the singlet Gross-Neveu coupling has been reversed. No first-order splitting was added.

A useful more general one-parameter choice is t_1=(k-1/2)alpha^2 and t_2=-(k+1/2)alpha^2, with the same p,s. While (8) remains positive, the mass-channel coefficient is

    (A-2B)/alpha^4 = 4k/3-7/24.

It crosses zero at k=7/32. This illustrates that the adjacent-bond resolvent, not only the single-port spectral weights, controls the perturbative mass-channel sign.

The positive sign in (12) addresses one dangerous continuum coupling. It does not prove a stable lattice CFT. Higher orders can break the accidental continuous symmetry, generate other marginal couplings, and shift Fermi points. Exact total parity is retained, but it alone does not establish infrared stability.

## 5. Reproduction and checked values

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 c8_resume_checks.py

This writes resume_checks.json, spectral_parent_checks.json, and spectral_parent_q3_m4.npz. The NPZ contains the complete 324-dimensional H_0, code isometry, port matrices, grading, spectral weights, and second-order operators. The physical order is L,R,environment; a port matrix acts on its named q-dimensional factor and as identity elsewhere. Reordering the environment between the two ports makes them geometric endpoints, but does not make H_0 nearest neighbor.

For the flat-parent three-block checks, (q,m)=(2,2) and (3,4), the maximum isometry/compression/second-moment residual is 1.2e-15. The full second-order residual is at most 4.9e-14, and the adjacent-bond residual is at most 2.1e-17. These are direct transition-amplitude contractions in the subspace reached from the code by one port insertion. The flat parent preserves this subspace, so the reduction is exact for the coefficient being checked.

For the two-energy q=3,m=4 parent, the dimension is 324, the gap is 1 to numerical precision, and the smallest spectral-weight eigenvalue is 0.4921875. Kernel, parity, spectral-weight, same-bond and adjacent-bond residuals are all below 1.0e-14. The measured coefficient ratios are A/alpha^4=25/24, B/alpha^4=1/3, and (A-2B)/alpha^4=3/8.

One implementation check initially failed because a spectral factor of the Choi matrix was used while the parity formula assumed its positive square root. Replacing that factor by the positive square root makes the environment grading match the written construction. The corrected full run passes. This was a dilation-basis mismatch, not an adjustment of the desired spectral weights.

No central charge was estimated for the finite-lambda block chain. These checks validate the algebra and the sign-changing construction, not a counterexample to c-d or a proof of microscopic criticality.
