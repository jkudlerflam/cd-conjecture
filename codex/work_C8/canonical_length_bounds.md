# Canonical energy increments and physical long-cell tails

Status: proved estimates under the explicit atomic and coordinate
hypotheses below. They apply to the bare legal C8 Hamiltonian. They do
not establish a conformal limit or an all-order effective action.

## 1. Atomic assumptions and an equal-payload raising operator

Use the exact coordinates of unprojected_cell_reduction.md. On a cell,
let n=ell-b, with integer ell>=2, and let A>=0 commute with ell. Its
kernel consists of two isomorphic finite spaces P0 and P1 at n=0 and
n=1. The logical identification between them is an isometry R with

    R dagger R=P0, R R dagger=P1,
    [A,R]=0, [n,R]=R, ||R||=1.

R is zero on the rest of the cell space. It is a mathematical local
trial operator in cell coordinates, not an additional microscopic
interaction. In C8 the two code spaces carry the same logical payload.

Write P=P0+P1, Q=1-P, and assume the form inequalities

    A>=Delta Q,
    A|ell >= gamma(1+|ell-b|) for ell outside {b,b+1}.       (1)

The weaker A|ell>=gamma(ell-b) for ell>=b+2 will be used for the tail.
With

    C=max(2/Delta,1/gamma),

the following useful single-cell inequality holds:

    |n-P1|+Q <= C A.                                      (2)

To check it, on P0 and P1 both sides on the left vanish. On excited
states of lengths b and b+1 the left side is respectively 1 and 2.
On all other lengths it is |n|+1. Equation (1) covers all cases.

On M cells let H=sum_j A_j+V, with bounded onsite or neighboring-cell
terms V_X, conserving N=sum n_j, and

    v=sup_j sum_(X contains j)||V_X||.

The coupling is included in v. An additive scalar per cell is omitted.
In particular ||V||<=Mv. No upper bound on A, cell dimension, or
retained spectral gap is assumed. Open boundaries are allowed.

## 2. Canonical increments from local trial states

[proved] Let E_M(N) be the open atomic ground energy at charge N,
and 0<=N<=M. There is an all-code trial state in this sector, so

    -Mv <= E_M(N) <= Mv,
    <A> <= 2Mv                                             (3)

in any canonical ground vector. With rho=N/M, (2) implies

    <sum_j P1_j> >= M(rho-2Cv),
    <sum_j P0_j> >= M(1-rho-2Cv).                          (4)

For the first inequality, N=sum P1_j+sum(n_j-P1_j) and the absolute
value of the second expectation is at most C<A>. For the second,
sum P0_j=M-sum P1_j-sum Q_j. The joint bound (2), rather than two
separate estimates, gives the stated same constant.

Apply R_j to a ground vector psi of charge N. The sum of trial norms
squared is sum_j <P0_j>. Each trial lies in charge N+1. The variational
principle, summed over j, gives

    [E_M(N+1)-E_M(N)] sum_j <P0_j>
      <= sum_j <R_j dagger [H,R_j]>
      <= 2Mv.                                             (5)

Here [A,R_j]=0 and only terms touching j contribute to [V,R_j].
Thus ||[V,R_j]||<=2v. The identity in (5) is real in a ground vector;
bounding the absolute value of each term is a valid upper bound.
Using R_j dagger instead proves the lowering counterpart. Consequently

    E_M(N+1)-E_M(N) <= 2v/(1-rho-2Cv),
    E_M(N-1)-E_M(N) <= 2v/(rho-2Cv),                       (6)

whenever the corresponding denominator is positive.

[proved] If eta in (0,1/2), v<=eta/(4C), and both N/M and N'/M lie
in [eta,1-eta], iterating the two inequalities gives

    |E_M(N')-E_M(N)| <= (4v/eta)|N'-N|.                   (7)

All intervening integer sectors stay in the same interval. This is
a Lipschitz energy estimate, not a differentiability or compressibility
claim. It survives level crossings, degeneracy and a critical payload.
The same proof holds for a fixed ordinary atomic boundary twist, but
we use the open chain below so that no twist is silently carried
between different physical lengths.

## 3. The physical quotient and a pinned-marker trial state

The coordinate space at physical length L contains x in Z/LZ and
ordered cell words of total length L. Its reference-label action is

    U:(x;u1,...,uM) -> (x+ell1;u2,...,uM,u1),   U^M=1.

The physical Hilbert space is the U-invariant subspace, with normalized
orbit sums. On the full coordinate space the lifted H commutes with U.
All terms preserve x except transfers across the reference marker,
which change x by one. Each cell has ell>=2. These are established
properties of the C8 coordinate representation.

[proved] A normalized vector phi supported at x=0 has

    Psi=M^(-1/2) sum_(a=0)^(M-1) U^a phi,
    ||Psi||=1,  <Psi,H Psi>=<phi,H phi>.                    (8)

For a nontrivial reference change the new coordinate is a nonempty
proper sum of cell lengths, between 2 and L-2. Thus its support cannot
overlap x=0, or reach x=0 under one application of H, which changes x
by at most one. All cross terms in the norm and Hamiltonian expectation
vanish. This argument allows superpositions of cell lengths.

It gives the physical trial states needed below without choosing a
physical momentum or assuming that an atomic ring ground state has
the required cell momentum. It concerns the lowest energy over all
physical momenta in the fixed-(M,L) legal sector.

[proved] In particular, insertion of a length-b code vector at cell 1
into a ground vector of the remaining open chain yields

    E_phys(M,L) <= E_(M-1)(N)+v,  N=L-bM.                 (9)

All terms meeting cell 1 have total norm at most v. Removing them
leaves the open Hamiltonian on cells 2,...,M, independently of x.
Their excess length is N because the inserted cell has excess zero.
Equation (8) makes this an admissible physical variational state.
Similarly, -Mv<=E_phys(M,L)<=Mv when 0<=N<=M.

Deleting only the seam bond gives

    |E_phys(M,L)-E_M(N)| <= v,                             (10)

where the open chain keeps all onsite terms. The lower bound holds
on the full lifted space before imposing U-invariance; the upper
bound uses (8). Thus (7) also bounds physical energy increments,
with an additional endpoint cost 2v. No smooth equation of state is
deduced. Extra fixed payload-charge or fixed momentum constraints
would need separate admissible trial states; they are not asserted.

## 4. A fixed-length factorial tail

[proved] Suppose, in addition, every V_X changes an individual cell
length by at most one. Choose eta in (0,1/4], M>=1/eta, and

    2eta <= rho=(L-bM)/M <= 1-2eta,
    v <= min(eta/(4C), eta gamma/16).                     (11)

For the lifted physical ground vector Psi, let

    Q_R=1_(ell_1>=R),  a_R=||Q_R Psi||.

Then for every integer r>=2,

    a_(b+r) <= (2v/gamma)^(r-1)/r!.                       (12)

This is uniform in M, L, the upper cell energy and cell dimension.
For U-invariant states the squared norm is the mean long-cell
occupation over physical markers. The statement does not select an
absolute physical position as a distinguished marker.

Proof. Remove every bounded term touching cell 1. The result is
A1 plus an open (M-1)-cell Hamiltonian. In a cell-1 length block
ell=b+s the rest has excess length N-s. Thus

    H-E_phys >= A1+E_(M-1)(N-s)-v-E_phys                  (13)

as a lower form bound after this decomposition; the subtracted v
accounts for all removed terms, including off-diagonal length moves.

First consider 2<=s<=eta M. All excess lengths between N-s and N
have densities in [eta,1-eta] on M-1 cells. The lower endpoint obeys
N-s>=eta M. For the upper endpoint, N<=(1-2eta)M and M>=1/eta
give N/(M-1)<=1-eta. Equations (7), (9) and (13) therefore give

    H-E_phys >= [gamma s-(4v/eta)s-2v] I

for this part of the diagonal lower form. By (11), 4v/eta<=gamma/4
and 2v<=gamma s/4, so it is at least gamma s/2.

For s>eta M, use the elementary lower bound H>=A1-Mv and upper
bound E_phys<=Mv. They give gamma s-2Mv >= (gamma-2v/eta)s,
again at least gamma s/2. Combining the two ranges yields

    Q_R(H-E_phys)Q_R >= (gamma/2)(R-b) Q_R,
                       R>=b+2.                           (14)

More formally, both comparisons lower-bound the same diagonal
operator A1+H_rest-v-E_phys, so they can be combined by its length
blocks before applying Q_R. No assumption that H itself commutes
with Q_R is being made.

The length bandwidth gives

    Q_R H(1-Q_R)=Q_R V 1_(ell_1=R-1),
    ||Q_R V 1_(ell_1=R-1)||<=v.

Project the ground-state equation and invert the positive compression
in (14). This gives a_R<=2v/[gamma(R-b)] a_(R-1). Iteration from
a_(b+1)<=1 proves (12). The zero-coupling limit has zero wrong-length
occupation, as it should. No inverse of the global bare-Q Hamiltonian
has been used.

## 5. True low-energy windows and a spectral cutoff consequence

[proved] Let Pi_Omega project onto the physical legal-sector energies
in [E_phys,E_phys+Omega], and lift its range to coordinate space.
If 0<=Omega<=gamma/2, then

    ||Q_(b+r) Pi_Omega||
       <= a_r^*= (2v/gamma)^(r-1)/(r-1)!,   r>=2.          (15)

This controls arbitrary superpositions in the window, with no factor
equal to its dimension. To prove it, set T=Q_R Pi_Omega and use the
Sylvester equation

    (Q_R H Q_R)T-T H_window = -Q_R H(1-Q_R)Pi_Omega.

The two diagonal spectra are separated by at least
gamma(R-b)/2-Omega, by (14). The integral solution, obtained by
integrating the decaying semigroup on the left and the contracting
semigroup of H_window-E_phys-Omega on the right, has norm at most
the source norm divided by this separation. Its source is bounded
by v||Q_(R-1)Pi_Omega||. The denominator is at least
gamma(R-b-1)/2. Iterating proves (15). This is a true spectral-window
statement, unlike inversion of an extensive unperturbed code band.

[proved] A useful finite-volume consequence is a controlled cutoff
on every cell. Let C_R=product_j(1-Q_R,j), R=b+r, and

    delta_R=sqrt(M) a_r^* <1.

The union bound in quadratic forms gives
||(1-C_R)Pi_Omega||<=delta_R. C_R is reference-invariant, so it
defines a physical restricted Hamiltonian H_R=C_R H C_R. If e_k is
any eigenvalue of H in this low window, counting multiplicity from
the bottom, then its counterpart in H_R obeys

    0 <= e_k(H_R)-e_k(H)
      <= [Omega delta_R^2+Mv delta_R]/(1-delta_R^2).        (16)

For clarity, the cutoff is a restriction, with no virtual correction.
Here is the min-max proof. On the spectral subspace through e_k,
let D=Pi(1-C_R)Pi, so ||D||<=delta_R^2. Since
||H_window-e_k||<=Omega and H_window-e_k<=0, the real part of

    Pi C_R(H-e_k)Pi=(I-D)(H_window-e_k)

is bounded above by Omega delta_R^2. Inserting the second C_R
subtracts Pi C_R H(1-C_R)Pi, of norm at most Mv delta_R, since
[A,C_R]=0 and ||V||<=Mv. Also ||C_R psi||^2>=1-delta_R^2.
The projected subspace has the same dimension, so min-max gives the
upper bound. Restricting the variational space gives the lower one.

For fixed v,gamma,b and M proportional to L, Stirling's formula shows
that r proportional to log M/log log M can make (16) smaller than
any specified inverse power of M. For example taking the proportionality
constant larger than 5/2 makes the error o(1/M). This controls the
unbounded length sum at a conformal finite-size energy scale with a
slowly growing cutoff. It does not supply a fixed-cutoff effective
Hamiltonian or prove that the spectrum is conformal.

## 6. Scope for C8 and remaining limitations

The hypotheses use the existing flat-wrapper gap and length growth,
the identical payloads at b and b+1, and bandwidth-one marker transfer.
Thus take Delta=Delta_at, gamma=gamma_wrong, and v=epsilon J with
J<=2B_flat. The reference sectors M nearest L/(b+1/2) satisfy (11)
for eta, for example, 1/5 and sufficiently large L. The required
epsilon is positive and can be fixed before taking L to infinity.

The theorem strengthens the earlier physical 1/(R-b) occupation
bound to a factorial tail and extends it to an entire low-energy
window. No one-step dressing is needed. The dressed Hamiltonian
generally has larger length bandwidth, so the proof is not copied
to it. It also does not choose the ground marker number, compare
illegal selector sectors, integrate out the massive states, bound
Gaussian stiffness, or identify conformal towers. Those remain
separate obligations. The original c-d question remains open here.

Finite checks and their scope are recorded by canonical_length_checks.py
and canonical_length_checks.json. The physical ring control is a
scalar variable-length cell model, not a new microscopic counterexample.
