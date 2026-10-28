# Joining the complete cell wrapper to a mobile selector

Status: proved local Hamiltonian, full unperturbed ring ground space,
volume-independent gap above that entire space, and first-order projected
operators. The theorem assumes an initialized uniform frustration-free
payload core. It does not supply a high-channel core or prove a CFT at
nonzero coupling. No microscopic counterexample is claimed.

## 1. Inputs and the uniform lower bound

Take the wrapper bond w of arrow_wrapper.md for an initialized q-state
core of length k>=3. Put b=k+1 and D=dim ker H_core. Its physical alphabet
has size s=6q+3. The wrapper parameters are

    0<v<U/(k+1),   mu_wrap=2U,   K_wrap=6U.

Here v denotes the vacancy hopping, to distinguish it from the marker
hopping introduced below. All positive filter and payload terms are as
specified in that note. On open lengths b and b+1 the COMPLETE ground
codes W_1,W_2 have energies

    epsilon_b=-4U+e_1,    epsilon_(b+1)=-4U+e_2,
    0<=e_p<=U.                                           (1)

They have positive gaps Delta_1,Delta_2 given there. The charge ground
vectors are chi_p and their endpoint occupation probabilities are r_p.

One further consequence of the wrapper proof is needed for EVERY open
length l, including lengths that do not carry the intended payload.

    H_wrap(l)>=-4U I.                                    (2)

For l=1 the sum is zero. For l>=2, a sector with z bad wrapper bonds
has at most 2(z+1) endpoint tags. Positivity of all terms other than the
endpoint rewards gives

    H_wrap(l)>=6Uz-4U(z+1)=2Uz-4U.

Every bad-bond projector commutes with the wrapper Hamiltonian, so this
is an operator bound on each complete sector. Minimizing over z proves
(2). On a periodic wrapper chain, the stronger endpoint count 2z gives
H_wrap,ring>=2Uz>=0. These bounds do not depend on the unknown core
spectrum at other lengths or on the norm of its positive penalties.

## 2. One repeated bond, including the resets

The microscopic site is a selector letter tensored with a wrapper
letter, of dimension

    d=4s=24q+12.                                         (3)

The selector letters #,L,X,R and their arrows are those of C7. Write
Q for its forbidden-pair projector and M=|#><#|. Let P_int be the
diagonal selector projector onto allowed pairs whose RIGHT letter is
not #. Define the unperturbed bond, with the factors grouped by track,

    h_base = K_sel Q tensor I
       + (mu/2)(M tensor I+I tensor M) tensor I
       - t(|LX><XR|+|XR><LX|) tensor I
       + P_int tensor w.                                (4)

The tensor grouping in (4) is just notation for a two-physical-site
operator. In particular P_int tensor w is still a nearest-neighbor
term. No site-dependent coefficient or imposed marker sector is used.

All selector violations Q_i commute with (4), as do all marker
positions at this stage. A head hop leaves both the marker positions
and the legality pattern unchanged, so it leaves every P_int unchanged.
Consequently the selector head Hamiltonian commutes with the whole
wrapper Hamiltonian in each such sector.

On a legal marker cell of length l, counting the first # and excluding
the next #, exactly its l-1 internal bonds have P_int=1. The boundary
bond into the next # has P_int=0. Thus the wrapper is an OPEN chain on
that cell, with no wrapper constraint across its boundary. This is the
required arrow reset. The head ground energy is -2t cos(pi/l), and the
complete cell ground energy is

    E_l=mu-2t cos(pi/l)+epsilon_l,   l>=2,                 (5)

where epsilon_l=min spec H_wrap(l). The internal head and payload
factorization is exact.

## 3. Selecting exactly b and b+1 despite the payload energies

Use dimensionless quantities

    a_l=2(1-cos(pi/l)),
    delta_b=(b+1)a_b-b a_(b+1),
    e_c=a_(b+1)-a_b<0,    mu_c=2-delta_b>0,
    eta_b=min{(a_(b-1)-2a_b+a_(b+1))/2,
              (a_b-2a_(b+1)+a_(b+2))/3}>0.             (6)

Here b>=4. Convexity of a(x) on x>=2 proves positivity of the second
differences. The positive value of mu_c at b=4 and monotonic decrease
of delta_b prove mu_c>0 for every b>=4.

Choose, for example,

    t=2U[1/eta_b+b/mu_c+1/|e_c|],
    mu=t mu_c+4U-(b+1)e_1+b e_2,
    e=t e_c+e_2-e_1,
    K_sel=4t+9U.                                        (7)

Then mu>0, e<0, and E_b=eb, E_(b+1)=e(b+1). All parameters depend
on the chosen core but not on the total ring length.

Here is an explicit lower bound for every competing length, without
solving its payload problem. Set

    d_l=a_l-a_b-(l-b)(a_(b+1)-a_b).

For l>=b+2, writing s_j=a_(j+1)-a_j and using increasing s_j gives

    d_l=sum_(j=b+1)^(l-1)(s_j-s_b)
        >=(l-b-1)(s_(b+1)-s_b)
        >=eta_b(1+|l-b|).

For 2<=l<=b-1,

    d_l=sum_(j=l)^(b-1)(s_b-s_j)
        >=(b-l)(s_b-s_(b-1))
        >=eta_b(1+|l-b|).

Both final inequalities follow respectively from (m-1)/(m+1)>=1/3
for m>=2 and m/(m+1)>=1/2 for m>=1. By (1),(2),

    E_l-el = t d_l + epsilon_l-epsilon_b
                      -(l-b)(epsilon_(b+1)-epsilon_b)
            >=(t eta_b-U)(1+|l-b|)>0                   (8)

for every l other than b,b+1. In particular the bound remains positive
at arbitrarily large l. The unknown payload ground energy cannot win
by forming very long cells.

## 4. Complete ground space on all sufficiently large rings

Suppose first that the selector sector is legal and contains markers.
The cells are independent at u=lambda=0. Equations (5),(8) give H>=eL.
Equality requires every cell to have length b or b+1, every head to be
in its unique sine-wave ground state, and every wrapper to be in the
complete code W_1 or W_2. Thus all these tilings and their D-dimensional
payloads exhaust the ground space in legal marker sectors.

Legal marker-free selector words are all L or all R. Their hopping and
marker energies vanish, and their wrapper is periodic and nonnegative
by (2)'s cyclic version. They lose to eL<0.

Now take z>=1 selector violations. Their positions are conserved.
Cutting at them gives z legal open strings, because P_int and the
selector kinetic term vanish on every cut. Within a string, complete
marker-to-marker cells have energy at least e times their length.
There are at most two remaining fragments, the prefix and suffix.
Each fragment has at most one head, with kinetic energy at least -2t,
and one open wrapper chain, with energy at least -4U. Any unassigned
marker energy is nonnegative since mu>0. Since e<0, a fragment of
length l has energy at least el-2t-4U. Summing gives

    H_sector>=eL+(K_sel-4t-8U)z=eL+Uz.                  (9)

This proof includes overlapping violations, fragments of length one,
and arbitrary illegal wrapper sectors. It does not require a bound
on ||w||. The positive filter scale can be large without affecting
the estimate.

Every L>=b(b+1) admits a b/(b+1) tiling. Indeed choose an integer M
in [L/(b+1),L/b], whose interval width is at least one; then
r=L-bM lies between 0 and M. There are r long cells and M-r short
cells. Consequently the full ring ground energy is eL for every such
L, with the complete ground space just described. It has no externally
selected marker, wrapper-tag, occupation, or program sector.

The gap ABOVE THIS ENTIRE SPACE has the L-independent lower bound

    Delta_base >= min {Delta_1,Delta_2,
       2t[cos(pi/b)-cos(2pi/b)],
       2t[cos(pi/(b+1))-cos(2pi/(b+1))],
       t eta_b-U, U, |e|b(b+1)} >0.                    (10)

The large degeneracy inside the space is not a CFT. Equation (10)
does not assert a gap after adding interactions within it.

## 5. A two-site marker hop with scalar payload transport

Let C=|H_2><H_0| on the wrapper. Add the nearest-neighbor bond

    h_move=-u[ |R#><#L| tensor(C tensor I)
                    +|#L><R#| tensor(C^dag tensor I)].  (11)

The first term moves the marker on site i to site i+1. Site i is
transferred from the right cell to the left cell, so its vacancy
must change from H_0 to H_2. This explains both the location and
orientation of C. The same local move preserves all selector
violations. It need not preserve wrapper violations or the low-energy
cell space.

Let P project onto the complete ground space of Section 4. A projected
hop is nonzero only on neighboring lengths (b,b+1), which it changes
to (b+1,b). The selector head overlap is A_b^2, with

    A_b=sum_(i=1)^(b-1) sqrt(2/b)sin(pi i/b)
                            sqrt(2/(b+1))sin(pi i/(b+1)).

The converted wrapper overlap, proved in arrow_wrapper.md, is beta
times the identity on both logical payloads, with

    beta=[sum_(s=0)^k chi_1({s})chi_2({s,k+1})]^2>0.

Therefore the exact first-order matrix element is

    P h_move P = -u A_b^2 beta I_(D^2)                  (12)

on each allowed length exchange. Physical marker-pattern coordinates
retain the ring translation and relabeling constraints, as in
mobile_selector.md. A local term (11), without an imposed long-range
transport operation, produces the desired logical identity.

## 6. Projected payload coupling and the remaining gap

Suppose the input core has V^dag A_a(first)V=alpha Gamma_a and
V^dag B_a(last)V=alpha R_a. Extend A_a only on wrapper L tags and
B_a only on wrapper R tags, acting trivially on the helix. Write these
as Atilde_a,Btilde_a. Add

    h_link=lambda P_cut tensor sum_ab T_ab
                            Btilde_a tensor Atilde_b,    (13)

where P_cut projects onto the allowed selector pairs X# and R#.
The coefficient matrix T is real. This is Hermitian, uniform, and
nearest neighbor. On cells with p,p' vacancies, its compression is

    lambda alpha^2 r_p r_(p') sum_ab T_ab R_a tensor Gamma_b. (14)

There is no additional head-wavefunction factor, since P_cut sums
both possible selector letters preceding a marker. In particular
it cannot be replaced by R# alone without changing (14).

Equations (12),(14) give the full first-order moving-cell interface.
The length hopping is payload independent, while the payload hopping
depends on adjacent lengths through r_1,r_2. Those probabilities are
usually unequal at finite v. Their difference cannot be dropped when
claiming an exact product of CFTs.

This construction removes the earlier local arrow-reset and sector-
selection gaps in joining the wrapper to the selector. It leaves two
major tasks: a high-channel initialized uniform FF core, and control
of the coupled length/payload theory at fixed u,lambda>0. A finite-order
projection plus (10) does not control thermodynamic gaps of order 1/L.
No finite-coupling central charge follows from this note.

## 7. Checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 selector_wrapper_checks.py

At k=3,U=1,v=0.1, formula (7) gives t=89.98213207197,
mu=56.92794991641, K_sel=368.92852828788, and e=-18.33424471347.
The guaranteed wrong-length coefficient t eta_b-U is 1.69357335860.
Complete legal wrapper spaces at lengths 2 through 6 have dimensions
31,101,315,961,2903. Their spectra obey the lower and selection bounds;
the two selected energies differ from eb,e(b+1) by less than 3e-13.
Illegal wrapper sectors are excluded analytically, not by full ED.

The constructed physical bond at d=36 is a 1296-by-1296 sparse matrix,
with Hermiticity residual 5.5e-14. Exhaustive selector-neighbor checks
give no violation of the overlapping legality projectors.

A separate quantum test uses a graded qubit-swap core with D=2 and
k=3, then tensors on the helix filter. Its wrapper alphabet has 15
states and its selector composition d=60. The physical code supports
have sizes 128 and 320 at the two lengths. Both endpoint maps agree
within 3.2e-16. The actual two-cell H_0-to-H_2 conversion gives beta
times I_4 with error 8.7e-15; the full head-plus-wrapper hopping error
is below 7.2e-15. Unconverted transport is zero. The marker-pattern
graphs at (L,M,r)=(18,4,2),(27,6,3),(36,8,4) agree with their projected
ground-energy formulas below 1.4e-14.

The output is selector_wrapper_checks.json. No full-ring microscopic
diagonalization or finite-coupling central-charge estimate was performed.
