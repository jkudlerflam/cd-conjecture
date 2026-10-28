# C25. Second-order composition of the C18 cell with the mobile wrapper

2026-10-28. G3 of `../../contrib/routeA_assembly/README.md`.

## 1. Result and scope

[Proved] The requested C18 Section 8 cell, with a frustration-free clock
Laplacian, has an exact two-dimensional code after the tagged-vacancy
extension. Its four occupied sites give wrapped lengths b=5 and b+1=6.
The second-order bond is diagonal and has a strictly negative ZZ coefficient.
Both physical lengths, all four ordered pairs of lengths, all vacancy
excitations, and the helix excitations reached by the compensators are
included in the calculation below.

[Proved] The existing choice m_p=1/r_p equalizes the leading source
amplitude. It does not equalize the complete resolvent at finite vacancy
hopping. With s the common scale of the core Hamiltonian, the result is

    s K2_(p,q) = Kbar2_(p,q)(w/U) + O(s/Delta_wrap),
    Kbar2_(p,q)(w/U) = K2_bare + O((w/U)^2).

The first constant depends on the compensator norm and the finite core.
Delta_wrap includes the smallest vacancy tunneling gap, not just U. Both
limits are finite-dimensional statements proved in Sections 7 and 8.
Thus the mismatch can be made arbitrarily small by first choosing w/U
small and nonzero, then choosing s/Delta_wrap sufficiently small. No
uniform limit with respect to these two independent ratios is asserted.

[Proved] The mixed coefficient at order u lambda is zero for this C18
cell. Bell-syndrome conservation gives the selection rule. It need not
hold for a syndrome-mixing replacement of the port check.

[Numerical, dimensions specified in Section 10] At w/U=0.06 and
s/Delta_wrap=0.001, the mismatch of the traceless bonds is 0.001327 times
the traceless unwrapped bond norm. The limiting mismatch at that w/U is
0.001336. The full bonds also contain length-dependent scalar terms.
They are reported separately, since they act on the length degrees of
freedom and cannot be discarded from the mobile chain.

[Proved] The proposed G1 syndrome-mixing positive check in
`../../contrib/routeA_demonstrator/cell.py` is compatible with a
frustration-free code when its parameters satisfy the strict positivity
condition in Section 11. This verifies the local interface condition,
not its missing two-mode controller or a high-channel compilation.

[Heuristic implication] These results support the assembly's proposed
second-order composition mechanism for this test cell. They do not kill
route A. They also do not close G3 for a future noncommuting hopping cell,
or establish a conformal phase. The exact equalization claim at finite
w and s would be false for this implementation.

## 2. Inputs, changes, and inherited assumptions

[Assumed inputs, explicit attribution] The C18 transition table and its
complete legal path are taken from `../work_C18/RESULT.md`, Sections 2,
3 and 8, and `../work_C18/clock_checks.py`. The tagged-vacancy bond,
endpoint tags, length filter, arrow-sector selection and selector reset
are those of `../work_C8/tagged_vacancy.md`, `length_filter.md`,
`arrow_wrapper.md`, and `selector_wrapper.md`. The helix tracks and
compensators are those of `../work_C8/helix_equalizer.md`. Their global
sector-exclusion theorems are inherited, not independently re-proved.
Every matrix identity used for the invariant-sector calculation is
specified and proved below. No full-ring calculation is claimed.

[Proved construction] Replace each C18 hopping rule by

    J (|a> tensor I - |a'> tensor U)
      (<a| tensor I - <a'| tensor U^dagger).

This convention gives forward matrix element -J U and is positive.
Set the dimensionless clock hopping and input-pin strengths to one,
and multiply the entire core by s>0. The Bell checks are
Q_Bell tensor (I+Z_c/4), also multiplied by s. The code history is
now uniform, not the sine history of C18's hopping clock. Consequently
its unwrapped second-order numbers differ from the C18 table.

[Proved construction] To regard this as an open initialized q=64 core,
use positive endpoint penalties requiring P_L on the first occupied
site and P_R on the last, together with C18's positive forbidden-pair
penalties. Tagged boundary penalties in C8 implement these conditions
locally even when vacancies sit at physical endpoints. They vanish
throughout the sector considered here and have no source matrix elements.
The original half-onsite convention for the C0/R1 turn is retained.
No replacement of a gated two-site check by an always-active check is made.

[Assumed scope] Ordered tags and the legal clock sector are used for
exact source calculations. The C8 filter and selector exclude other
occupations and cell lengths under their stated hypotheses. This run
does not reassemble or count the whole repeated microscopic bond.
No calculation on C21/C22's compiler, C26/C27's gaps, or C8's readout
program is duplicated.

## 3. The four-site core in an explicit clock basis

[Proved] Write the occupied sites as P_L, bulk 1, bulk 2, P_R.
The legal clock has t=0,...,29. The only work gates are SWAP on
edges t=4,10,18,24. Let U_t be their cumulative product and gauge a
history vector by U_t^dagger. In this frame all propagation gates are
identity. The initialized work inputs are |i,0>, i=0,1.

Define L_t=(|t>-|t+1>)(<t|-<t+1|). Split the clock Laplacian
into its three physical occupied-bond contributions B_1+B_2+B_3.
The following lists give the coefficient of L_t. Unmarked entries
have coefficient one.

    B_1: 3,7,8/2,9,13,17,21,22/2,23,27
    B_2: 0,2,4,6,8/2,10,12,14,15/2,16,18,20,22/2,24,26,28
    B_3: 1,5,11,15/2,19,25

Here 8/2 means half of L_8, not L_4. The half terms are the original
one-site turn, split between its two neighboring bonds. Each B_j is
positive and annihilates phi=30^(-1/2) sum_t |t>.

[Proved] Let epsilon_t be one at t=5,...,10 and 19,...,24 and zero
otherwise. These are the times when the work bits have been swapped.
On the initialized work subspace the nonzero Bell-syndrome potentials
are diagonal in t and i:

    v_L,0(t)=v_R,0(t)=5/4,
    v_L,1(t)=5/4-(1-epsilon_t)/2,
    v_R,1(t)=5/4-epsilon_t/2.

Both Bell pairs in the code are Phi+. The input pin acts at t=1 on
the second input bit. It vanishes on |i,0>. Thus the full legal core
has dimension 30*4*16=1920 and its code has dimension two. Positivity
forces a zero vector to have a constant gauged history, both Bell pairs
Phi+, and second work input zero, proving completeness in this sector.

[Proved] A Pauli acting on either port maps its Bell pair to one of
three orthonormal nonzero syndromes. Therefore every traceless port
operator has zero code compression. A source insertion preserves the
gauged input bits because the clock gates touch only the two c bits.
The subsequent diagonal Bell potentials also preserve these inputs.

## 4. Exact tagged-vacancy matrices

[Proved construction] Put p=1 or 2 vacancies among N=4+p physical sites.
A charge configuration S is a p-element subset of {0,...,N-1}.
Let o_1<o_2<o_3<o_4 be the occupied positions and define

    D_j(S)=1 if o_(j+1)=o_j+1, and 0 otherwise,
    O_L(S)=1 if 0 is occupied, and 0 otherwise,
    O_R(S)=1 if N-1 is occupied, and 0 otherwise.

The charge Hamiltonian is

    H_ch(S,S)=U * number of occupied-vacant interfaces,
    H_ch(S',S)=-w if S' is one adjacent vacancy-particle exchange.

Its dimensions are 5 and 15. Denote its ground energy by e_p and its
normalized positive ground vector by chi_p. The configuration graph
is connected for w>0. Shifting a sufficiently large multiple of the
identity gives an irreducible entrywise nonnegative matrix, so the
maximal eigenvector is strictly positive and unique. This proves the
stated properties of chi_p. Reflection symmetry gives

    r_p=<chi_p|O_L|chi_p>=<chi_p|O_R|chi_p>.

[Proved] In a source block with active end a and logical bit i, the
complete charge-clock Hamiltonian, measured from the ground energy, is

    H_(p,a,i) = (H_ch-e_p) tensor I_30
        + s sum_(j=1)^3 D_j tensor B_j
        + s D_(edge(a)) tensor diag(v_a,i),                 (1)

where edge(L)=1 and edge(R)=3. Its dimension is 150 for p=1 and
450 for p=2. In particular a vacancy separating a port from its
partner removes that Bell check. An internal vacancy can also remove
clock propagation. These effects are retained in (1).

[Proved] Before a source insertion, remove the last term of (1).
Then chi_p tensor phi is annihilated, term by term after subtracting
e_p from H_ch. Restoring all Bell sectors and the input pin yields
exactly chi_p tensor the two-dimensional core code. To see completeness,
zero energy first requires the charge vector chi_p. Every D_j has
strictly positive expectation in chi_p, so positivity then forces each
original core term to annihilate the remaining vector. The argument
also includes the positive endpoint penalties.

[Proved] The occupied-site length filter can remain in its product
ground vector throughout a source excursion. Every gated filter term
annihilates that vector separately, vacancy hopping preserves the ordered
filter spins, and the port source acts as identity on them. Its omission
from (1) is an exact invariant restriction. The same is true of the
unused helix track and of the selector head for two port insertions.
The cut projector includes both X# and R#, so no head-wavefunction
factor multiplies K2.

## 5. Helix compensators without discarding their excitations

[Proved construction] Use helix pitch phi_h. For lengths 5 and 6, the
left endpoint of the B track has states

    u_1=|+(-4 phi_h)>,  u_2=|+(-5 phi_h)>,
    chi_h=sin^2(phi_h/2),  m_p=1/r_p,
    M=[m_1(I-|u_2><u_2|)+m_2(I-|u_1><u_1|)]/chi_h.        (2)

Taking expectations proves <u_p|M|u_p>=m_p. The right endpoint uses
M^*. The actual port operator is O_a F_a tensor M_a. Thus its source
vector has charge-clock factor

    z_(p,a)=(O_a chi_p) tensor phi.                         (3)

Its projection onto the charge ground vector has amplitude r_p, not
sqrt(r_p). The corresponding helix amplitude is m_p. Their product
is one, which is the proposed leading equalization mechanism.

[Proved] Untwist the B track and rotate its common |+> direction to
|0>. Its parent becomes sum of singlet projectors plus |1><1| on
the right endpoint. The all-zero vector is its ground state. Applying
an arbitrary one-qubit M to that vector creates only zero and one
flipped spin. The parent preserves flipped-spin number. Thus the full
reachable helix sector has dimension N+1, without approximation.

[Proved] On the one-flip sector the parent is the N by N matrix

    T_N = (1/2) L_path,N + |N-1><N-1|.                    (4)

Let h_l and v_l be its eigenvalues and normalized eigenvectors. With
u=|+(-(N-1)phi_h)> and u_perp=| - (-(N-1)phi_h)>, the helix source
measure has an atom of weight m_p^2 at zero and atoms

    |<u_perp|M|u>|^2 |v_l(0)|^2 at energy h_l.             (5)

Reflection and conjugation give the same measure at the opposite end.
The source spaces therefore have dimensions 150*6=900 and 450*7=3150.
The other helix track stays in its ground state. The complete many-spin
helix Hilbert spaces have not been replaced by expectation values.

## 6. Second-order bond and its strict ZZ sign

[Proved] Diagonalize (1) as H_(p,a,i) v_alpha=E_alpha v_alpha and set
c_alpha=|<v_alpha|z_(p,a)>|^2. Convolve this measure with (5), calling
the resulting energies and weights e_mu and a_mu. Each E_alpha is
strictly positive. Positivity plus the charge-ground argument would
otherwise require phi to have zero expectation in a nonzero Bell
potential, which is impossible.

[Proved] The second-order effective Hamiltonian is lambda^2 K2 with

    K2_(p,q)(i,j) = -3 sum_(mu,nu)
           a_(p,R,i),mu a_(q,L,j),nu /(e_mu+e_nu).        (6)

All off-diagonal logical entries vanish. The factor three is the three
orthogonal Pauli-created Bell syndromes. For example, in the basis
{Phi+,X Phi+,Y Phi+,Z Phi+} at each active Bell pair,
(sum_a sigma_a tensor sigma_a)|0,0> = sum_(a=1)^3 |a,a>.
No interference term is omitted. The resolvent is a Kronecker sum,
so its spectral eigenvalues are the sums in (6). This also proves
(6) directly from -P V Q (QH0Q)^(-1) Q V P.

[Proved] For an equivalent integral formula let
f_(p,a,i)(t)=<z_(p,a)|exp(-t H_(p,a,i))|z_(p,a)> and let h_p(t)
be the helix source heat kernel from (5). Then

    K2_(p,q)(i,j)=-3 integral_0^infinity
                  f_(p,R,i)(t) f_(q,L,j)(t) h_p(t) h_q(t) dt.

The only difference between i=0 and i=1 in (1) is a nonnegative,
nonzero diagonal reduction. The joint charge-clock graph is connected:
charge moves connect all configurations, and an end-vacancy configuration
allows every clock edge. The negative off-diagonal entries consequently
give strictly positive heat kernels for t>0. Duhamel's identity gives

    exp(-t H_1)-exp(-t H_0)
      = integral_0^t exp(-(t-v)H_1)(H_0-H_1)exp(-v H_0) dv.

Its expectation in the nonnegative, nonzero vector z is positive.
Therefore f_1>f_0 for t>0. Also h_p(t)>0. It follows that

    beta_ZZ=-(3/4) integral
       (f_(p,R,0)-f_(p,R,1))(f_(q,L,0)-f_(q,L,1))h_p h_q dt <0.

The positive source gaps prove convergence. This establishes the
nonzero interaction for both lengths without numerical inference.

[Numerical, 30-dimensional unwrapped source blocks] With s=1 the
uniform-history unwrapped bond has diagonal

    (-1.2000000000, -1.3776084330, -1.3176089428, -1.5396496795)

and coefficients in the order I, ZI, IZ, ZZ

    (-1.3587167638, 0.06991254734, 0.09991229243, -0.01110807592).

[Proved consistency check] For input 00, both potentials are 5/4 and
phi is a clock zero vector, so K2_bare(0,0)=-3/(5/4+5/4)=-1.2.
This matches the numerical normalization.

## 7. The limit of a fast wrapper, with a finite bound

[Proved] Fix w>0. Project onto the charge ground vector and both
helix ground vectors, leaving the whole clock space intact. Define

    d_(p,j)=<chi_p|D_j|chi_p>,
    Hbar_(p,a,i)=sum_j d_(p,j) B_j
                       +d_(p,edge(a)) diag(v_a,i).        (7)

The projected source is phi because r_p m_p=1. Thus the limiting
bond Kbar2 is obtained from (6) using the 30-dimensional operators
(7) and source phi. This is generally different at the two lengths.

[Proved bound] Here is a direct estimate that also accounts for source
components outside the wrapper ground space. In a fixed pair-source
sector write the full Hamiltonian as A+sW, with A the shifted charge
and helix Hamiltonians, P its kernel projector and Q=I-P. Put

    QAQ >= delta Q,  B=PWP >= a P,  ||W||<=L,
    C=QWP,  D=Q(A+sW)Q,
    x=P z_pair,  y=Q z_pair,  t=s/delta.

Both W and A are positive. For this cell one may use L=10.5, since
each gated clock has norm at most 4 and its active Bell potential
has norm at most 1.25. One can use

    a=(3/4)(d_(p,edge(R))+d_(q,edge(L)))>0.

The Schur complement and completing the quadratic form give exactly

    s <z_pair|(A+sW)^(-1)|z_pair>
      = (x-s C^dag D^(-1)y)^dag
          (B-s C^dag D^(-1)C)^(-1)
          (x-s C^dag D^(-1)y)
        + s <y|D^(-1)|y>.                                (8)

Since ||D^(-1)||<=1/delta, if t L^2<=a/2, the inverse in (8)
has norm at most 2/a. The resolvent identity bounds its difference
from B^(-1) by 2t L^2/a^2. Taking the difference of quadratic forms
in (8) therefore bounds its absolute error from <x|B^(-1)|x> by

    2t L^2 ||x||^2/a^2
      +(2/a)(2t L ||x|| ||y||+t^2 L^2 ||y||^2)
      +t ||y||^2.                                       (9)

Multiply by three for (6). The source norm is finite and bounded by
the product of the endpoint compensator norms. This proves the first
asymptotic statement in Section 1. A small ratio s/delta alone does
not supply a norm-independent error estimate for a large compensator.

[Proved] In this calculation delta can be the minimum of the two
charge gaps and the two helix gaps. The latter are bounded below by
1/[N(2N-1)], so the charge gap is the limiting scale in the scan.
Other excluded sectors of a complete selector construction can impose
smaller permissible microscopic perturbations. They are not additional
states reached in (6).

## 8. Small vacancy hopping and the missing-bond correction

[Proved] At w=0 the charge potential has two minima of energy U.
They have all vacancies in a single block at the left or right end.
Every D_j equals one on both configurations. All other configurations
have energy at least 2U. Write P_end for these two states and T for
the unit-hopping charge adjacency matrix, with ||T||<=N-1.
For chi=chi_p, E=e_p<=U follows by an end-state trial vector.
The Q_end eigenvector equation implies, whenever w(N-1)<U,

    ||Q_end chi||
        <= w(N-1)/[U-w(N-1)].                            (10)

Indeed Q_end H_ch Q_end-E >= [U-w(N-1)]Q_end, and the source
term is w Q_end T P_end chi. Because (I-D_j)P_end=0,

    0<=1-d_(p,j)<=||Q_end chi||^2=O((w/U)^2).             (11)

The bound is only useful when its right side is small, but its
asymptotic statement holds. The minimum extra interface cost is U,
not 2U as suggested by the assembly's informal description.

[Proved] Substitute (11) in (7). Every term is a fixed finite matrix,
so ||Hbar-H_bare||=O((w/U)^2). The two-cell source Hamiltonians have
positive limiting lower bound 3/2. The identity
B^(-1)-B0^(-1)=B^(-1)(B0-B)B0^(-1) proves the same bound for Kbar2.
It also proves the length mismatch is O((w/U)^2). Together with (9),
this establishes the ordered two-scale limit claimed in Section 1.

[Proved] Small w also makes the charge gap small. Any path between
the two end configurations needs 4p nearest-neighbor exchanges,
because the sum of vacancy positions changes by 4p. In the Schur
complement onto P_end, the off-diagonal matrix element consequently
first appears at order w^(4p). Its leading terms have one sign and
positive denominators from the states outside P_end. Reflection makes
the two diagonal entries equal. The resulting even/odd splitting
is positive and of order w^(4p)/U^(4p-1). This supplies a finite
perturbative proof of the exponents, independent of diagonalization.

[Numerical finite path sum] The script also sums shortest-path weights
on the 5-state and 15-state graphs. It gives leading coefficients 2
and 1/4 for p=1 and p=2. The precise coefficients are not needed for
the proof of either asymptotic bound. In particular the two-vacancy
gap scales as w^8, much faster than the missing-bond weight w^2.

[Proved scope qualification] The model contains no C19 Jacobi filters.
It does not test the assembly's additional cancellation from
sum_e X_ea=0. The bound above does not use that cancellation.

## 9. The mixed u lambda term

[Proved] Let P_x project onto the code at a marker configuration x.
Write T for the marker move, including vacancy-type conversion and
helix rotations, and V for the port exchange. At second order the
mixed matrix element between two code configurations is

    -u lambda P_y [T R V + V R T] P_x,                    (12)

where R is the reduced unperturbed resolvent in the intermediate
marker configuration. The same formula follows by expanding the
quadratic expression in uT+lambda V.

[Proved] Identify each cell's occupied particles in their order.
T transfers only a vacancy, so it preserves every port-partner Bell
state. The helix factors and selector head do not act on these Bell
pairs. Every unperturbed core check is diagonal in the Bell labels,
so R preserves them as well. V applied to a code creates one nonzero
Bell syndrome at each of two ports. With only one V these syndromes
cannot return to zero. Therefore each term in (12) vanishes separately.
The argument includes moves on neighboring cuts and marker moves
outside the two selected-length sectors. It does not assume that T
has no leakage.

[Proved limitation] The zero is particular to this syndrome-conserving
cell. Replacing its check by H_pair in Section 11 loses this argument.
Such a cell requires its own mixed-term calculation. No generic
improvement of the assembly's u lambda estimate is claimed.

## 10. Numerical results and verification

[Numerical specification] `composition.py` computes all four ordered
length pairs 11,12,21,22 at U=1, w=0.16,0.10,0.06 and
s/Delta_wrap=0.1,0.01,0.001. For each it uses both the positive
compensator (2) at pitches pi/2 and 0.2, and the smaller Hermitian
compensator of C8 at pitch 0.2. This is 27 parameter points and 108
second-order bonds. Every source eigenvalue and weight is retained.
The largest two-cell source product has dimension 3150^2=9,922,500
per Pauli channel. Its inverse is contracted through the exact sum
formula (6); no dense matrix of that size is formed.

[Numerical specification] The mismatch estimator is

    E_len=max_(p,q) ||tl[s K2_(p,q)-s K2_(1,1)]||
                       /||tl[K2_bare]||,
    tl[X]=X-Tr(X) I/4.

All matrices are diagonal, so the norm is the largest absolute diagonal
entry. The denominator is 0.1809329157. E_limit compares each s K2 to
its own Kbar2, with the same traceless norm and denominator.
The following table uses the positive compensator at pitch pi/2.

| w/U | Delta_wrap/U | s/Delta_wrap | E_len | E_limit |
|:--|:--|:--|:--|:--|
| 0.16 | 1.123e-7 | 0.1 | 0.055413 | 0.065154 |
| 0.16 | 1.123e-7 | 0.01 | 0.0088488 | 0.00089243 |
| 0.16 | 1.123e-7 | 0.001 | 0.0097320 | 0.0000092402 |
| 0.10 | 2.544e-9 | 0.1 | 0.061615 | 0.065357 |
| 0.10 | 2.544e-9 | 0.01 | 0.0028463 | 0.00089559 |
| 0.10 | 2.544e-9 | 0.001 | 0.0037326 | 0.0000092725 |
| 0.06 | 4.226e-11 | 0.1 | 0.064078 | 0.065414 |
| 0.06 | 4.226e-11 | 0.01 | 0.00076203 | 0.00089454 |
| 0.06 | 4.226e-11 | 0.001 | 0.0013266 | 0.0000092766 |

[Numerical] The limiting E_len values computed directly from (7) are
0.00974125, 0.00374192, 0.00133587. Their ratios to (w/U)^2 are
0.3805, 0.3742, 0.3711. These checks agree with the proved quadratic
bound. The finite-s mismatch need not decrease monotonically, because
its two corrections can cancel. The approximately quadratic E_limit
in this scan is an observation, not a replacement for the bound (9).

[Numerical] At w/U=0.10 and s/Delta_wrap=0.001, the coefficients of
s K2 in the order I, ZI, IZ, ZZ are

| length pair | I | ZI | IZ | ZZ |
|:--|:--|:--|:--|:--|
| 11 | -1.365598444 | 0.0702613804 | 0.1004132672 | -0.01116219287 |
| 12 | -1.366077769 | 0.0699944107 | 0.1004342704 | -0.01114222265 |
| 21 | -1.365934405 | 0.0702672161 | 0.1000173450 | -0.01114146730 |
| 22 | -1.367921193 | 0.0700007732 | 0.1000390613 | -0.01112164790 |

[Numerical] For the same w, the largest full-matrix difference from
Kbar2, without removing the scalar, is 0.63615, 0.073593, 0.0074780
as s/Delta_wrap decreases through the three scan values. It is much
larger than the traceless difference and decreases approximately
linearly. Formula (9) controls this term too. It is a length interaction
of order lambda^2/s times the displayed number, so a mobile application
must include it in its error budget. The full length mismatch at the
smallest ratio is 0.0029170 in the same scaled units. Writing p=1+n_left
and q=1+n_right, the scalar coefficient contains a density interaction
(c_22-c_21-c_12+c_11)n_left n_right. At this parameter point its
coefficient in s K2 is -0.00150746. Constant and linear length terms
can be shifted or retuned, but this term is an actual interaction.

[Numerical] At w/U=0.10 the mean adjacency factors are

    p=1: (0.9949540107, 0.9998021758, 0.9949540107),
    p=2: (0.9987279900, 0.9999935734, 0.9987279900).

These unequal factors locate the residual mismatch directly. The
endpoint occupation weights are 0.5051449013 and 0.5000126271.
At pitch 0.2 the positive compensator norm is 398.29, whereas the
smaller Hermitian choice has norm 2.0916. At pitch pi/2 the positive
choice has norm 6.7936. The traceless mismatch is unchanged to the
shown table precision for these three choices. This is numerical
agreement in the tested hierarchy, not exact compensator independence.

[Numerical stability and checks] The two charge reflection sectors
are diagonalized separately. A raw initial diagonalization failed the
reflection check by about 7e-9 at p=2,w=0.16, owing to the tiny even/odd
splitting. That implementation was corrected before the scan. A stable
block reduction is also necessary to avoid subtracting O(U) numbers to
obtain O(s) source energies. In the charge eigenbasis the dimensionless
source matrix is [[A,B^dag],[B,D]]. The graph of its low invariant
subspace obeys

    B+D X=X(A+B^dag X).

The code solves this equation to residual below 2e-11, normalizes
(I,X) and its orthogonal complement (-X^dag,I), and diagonalizes both
blocks. Both blocks and their source weights enter (6). Thus this is
a stable invariant-subspace computation, not low-block truncation.
The very smallest charge gap is resolved in double precision; only
four significant digits of that gap are reported. No certified
interval arithmetic is claimed.

[Numerical, independent comparisons] `checks.py` gives the following.

* Raw inversion of 450-dimensional charge-clock matrices agrees with
  the block spectral calculation for s times the resolvent quadratic
  form to 1.72e-9, at both ends, both inputs, and three spectral shifts.
  The wrapped zero-code residual is 1.66e-16.
* Full one-track helix matrices of dimensions 32 and 64 agree with the
  reduced 6- and 7-dimensional heat kernels at t=0,0.1,1,10 to 1.37e-12.
* The literal 16-dimensional pair-Bell source has precisely the three
  nonzero-syndrome components used in (6), with residual 7.70e-16.
* The literal G1 check and complete 32-dimensional demonstrator cell
  satisfy the positivity and code tests in Section 11. Its gap at
  D0=D1=1,kappa=0.5 is 0.5, its kernel dimension is two, and its code
  residual is 3.26e-17.

[Numerical scope] No central charge, finite-coupling many-cell spectrum,
or thermodynamic limit is measured. The 27-point scan took about
36 seconds with one BLAS thread. All checks use numpy and scipy.

## 11. Frustration freeness of the proposed G1 port terms

[Proved] In the notation of the demonstrator, write A=Phi-, B=|01>,
C=|10>, and let Phi+ be the zero-syndrome Bell state. The proposed
operator is

    H_pair=D0 |A><A| tensor I
       +D1 (|B><B|+|C><C|) tensor I
       +kappa (|A><B| tensor s+ - |A><C| tensor s- + h.c.).

The two states Phi+ tensor |i> are annihilated. On their orthogonal
complement the operator splits into two 2 by 2 matrices

    [[D0, kappa],[kappa,D1]],
    [[D0,-kappa],[-kappa,D1]],

and two uncoupled entries D1. Therefore it is positive for
D0,D1>0 and kappa^2<=D0 D1. If the inequality is strict, its only
kernel is Phi+ tensor C^2. Equality can add unwanted zero states
and must not be used to assert a complete initialized code.

[Proved] The left and right checks on a cell with separate Bell pairs
may share the logical work qubit and need not commute. Each nonetheless
annihilates the product of both Bell pairs and any logical vector.
Strict positivity on each nonzero Bell sector forces the intersection
of their kernels to be exactly that product code. Positive propagation
Laplacians with gates acting only on the work register annihilate its
history termwise, and compatible positive input pins preserve this
property. Hence syndrome mixing by this check does not force a
non-frustration-free core and does not violate the wrapper's hypothesis.

[Assumed remaining requirement] A future G1/G2 compiler must actually
supply those work-only gates, initialization conditions and complete
code. This calculation verifies the displayed check, not an unspecified
compiled Hamiltonian. The original unbounded_c/G1 abstract swap gadget
is a separate construction and is not established as a local FF core here.

## 12. Reproduction and handoff

[Numerical reproduction]

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B composition.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py

The files `composition_results.json` and `checks_results.json` contain
all bond diagonals, coefficients, source dimensions, scale values and
verification outputs. Imports of C18's rule table and the demonstrator
are read-only, with bytecode generation disabled. All generated files
are inside this folder. No network, cluster, subagent, commit or push
was used. The project STATE is outside the permitted write folder,
so the handoff is in `PROGRESS.md` and local `STATE.md`.

[Proved conclusion] G3's named small-cell calculation has a positive
answer in the ordered two-scale limit. The finite-w obstruction is the
mean-adjacency change of the core spectrum. It prevents exact second-order
equalization by m_p=1/r_p alone, but it is tunable and does not exclude
route A. The mixed u lambda term vanishes in this test. The next composition
calculation should use the actual hopping-producing G1 cell once supplied,
retain its noncommuting syndrome dynamics, and include the scalar length
interaction as well as the traceless payload mismatch.
