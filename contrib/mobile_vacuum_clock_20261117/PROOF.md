# Vacuum-clock control without a charged-mass hierarchy

Date: 2026-11-17. **[Proved below; not independently reviewed.]** The
results concern the exact legal ordered-material sector of the literal
22-state model. They do not select material lengths or exclude illegal
sectors. The reference ensemble is the positive length tilt used in the
[previous mobile contribution](../small_d_mobile_matching_20261116/PROOF.md).
Its fixed-vacancy minimizing blocks are covered; arbitrary prescribed
canonical densities are not silently substituted.

## 1. The exact vacuum clock and its common kernel

A cell with n bits occupies b=n+1 material particles, starting at a and
ending at marker m=a+n. Decode its vacuum space into clock positions
0,...,n and its 2^n payload bits. Put

    e_j=q|j>-|j+1>,  E_j=|e_j><e_j|,  j=0,...,n-1,
    Q_ff=sum_j E_j,
    Z_n=sum_(j=0)^n q^(2j),
    u=Z_n^(-1/2) sum_(j=0)^n q^j |j>,  P_0=|u><u|.

The exact contact-weighted edge Hamiltonian is

    J sum_(j=0)^(n-1) A_(a+j) E_j.                     (1)

Every E_j annihilates u. Consequently P_0 tensor the payload identity
commutes with (1), with the entire gap hopping Hamiltonian, and with
all other cells at zero port. The marker chemical term is

    J mu_sel C_m,   C_m=(A_(m-1)+A_m)/2,                (2)

which is identity on this cell's clock and also commutes with P_0.
The adjacent junction is identity on the marker state. Thus the common
code is exactly invariant before adding the terminal term, even with
the true surrounding interacting cells and all their contacts retained.
No averaging or independence of contacts is involved.
The single-cell statement uses its ordered root frame. The product of
these code projectors on all cells also descends to the physical cyclic
space; a root-frame cell label is not an additional physical state.

The nonzero eigenvalues of Q_ff are

    1+q^2-2q cos(r pi/(n+1)),  r=1,...,n.

This follows by taking the Gram matrix of its n edge vectors: it has
diagonal 1+q^2 and off-diagonal -q. Define

    d=1+q^2-2q cos(pi/(n+1)),
    c=1+q^2+2q cos(pi/(n+1)).                           (3)

Then d(1-P_0)<=Q_ff<=c(1-P_0), with d>=(1-q)^2.
These are clock identities, not a many-body gap assumption.

## 2. The necessary terminal term breaks exact invariance

The literal terminal energy is

    T=J q^2 C_m |n><n|.                                (4)

The contact half on the marker's outgoing bond must be retained.
Set

    p_n=q^(2n)/Z_n,   tau_n=q^2 p_n,
    epsilon_n=q^2 sqrt(p_n(1-p_n)).                     (5)

On the vacuum block,

    (1-P_0) H P_0 = J q^2 C_m (1-P_0)|n><n|P_0,
    ||(1-P_0) H P_0|| <= J epsilon_n                   (6)

at lambda=f=0. Equality holds when the gap space permits C_m=1.
For q=1/4, epsilon_1=1/68 and epsilon_2 is about 0.003776.
They are not made small relative to J by sending J/t to zero.

There is no alternative clock-only rank-one projector invariant under
all these contact coefficients. Independent variation of A_m first
requires invariance under |n><n|, and variation of the internal contacts
then requires invariance under every E_j. Adjacent edge vectors have
nonzero overlap -q. A common eigenvector of the rank-one E_j must
therefore annihilate all of them, hence be u. But |n><n| does not preserve
the line through u for n>=1 and q>0. For n=1 this is the same elementary
fact for two noncommuting rank-one projectors in dimension two.

The frozen eigenvector of Q_ff+q^2|n><n| consequently is not an exact
mobile code either. At small J the sea compression selects that frozen
vector, rather than u; fast vacancies do not by themselves remove this
bare-projector difference.

Deleting (4) is not a valid simplification of the selector. It makes
every zero-port cell energy equal to mu_sel. Negative mu_sel then selects
zero-bit cells, zero mu_sel makes lengths degenerate, and positive mu_sel
favors marker-free/infinitely long cells. No prescribed n>=1 is selected.

## 3. A channel that retains coherent clock gradients

Let P_perp=(1-P_0) on this cell's vacuum space, tensored with its payload.
Compress the full Hamiltonian to P_perp on this cell only. Other cells
retain all their vacuum and charged states and their mutual ports.
Incident ports vanish in the compression, since they change this cell
from vacuum to charged. Let H_s=H+s sum_i g_i, t/2<=s<=3t/4.

Use the following Kraus operators on this pinned space:

    K_(j,0)=A_(a+j) sqrt(E_j/c),
    K_(j,1)=S_(a+j)^* sqrt(E_j/c),
    K_* = sqrt(1-Q_ff/c).                               (7)

S_i raises gap g_i by one. Since A_i+S_i S_i^*=1, the Kraus squares
sum to the identity. Every E_j preserves P_perp, so the channel is
admissible in the pinned sector. It resolves the clock's *gradients*,
not its position; these gradients remain coherent rank-one vectors.
It does not assert that neighboring E_j are orthogonal.

For any input state define

    X=<sum_j A_(a+j) E_j>,  W=<Q_ff>,  Y=W-X>=0.

The channel removes one vacancy with probability Y/c. The exact
one-gap compression identities give a kinetic energy change bounded by

    (4t/c) sum_j sqrt(<A_(a+j) E_j><(1-A_(a+j))E_j>)
      <= (4t/c) sqrt(XY).                              (8)

This is Cauchy--Schwarz on the actual joint expectations. No product
<A_j><E_j>, scalar impurity substitution, or bath spectral gap is used.
It remains valid with a uniform gap-hopping phase.

## 4. The local error does not contain D, Delta or the cell length

All terms outside the target cell that commute with (7) are left intact.
The only source terms that can change in the pinned sector are:

- the target edges (1), whose norm is at most c;
- the terminal term, of norm at most q^2;
- the marker half mu_sel A_(m-1)/2, of norm at most |mu_sel|/2;
- the cell's readout, a positive operator of norm at most nf.

The other marker half mu_sel A_m/2 commutes with the channel: none of
its branches shifts gap m. Charged masses and charged internal moves
vanish on this vacuum cell. At neighboring cells they are unchanged.
The boundary junction is identity on the target marker; incident ports
have already vanished under the pin. This accounts for their actual
contact weights, rather than bounding every nearby full bond by D.

Write

    r=c+q^2+|mu_sel|/2+nf.                              (9)

The source energy increase is at most 2Jr. For a pinned ground state
the variational principle, (8), and the length tilt imply

    sY <= 4t sqrt(XY) + 2cJr.                          (10)

For any 0<theta<1, suppose

    Jr <= theta s d/(4c).                              (11)

On the pinned space W>=d. If X<=W/2, then Y>=d/2 and
2cJr<=theta sY. Equation (10) yields

    X >= d [(1-theta)^2 s^2]/[(1-theta)^2 s^2+16t^2].   (12)

If X>W/2, the same bound holds in the stated source window. Thus (12)
is uniform in volume and in all surrounding cells, including arbitrarily
long surrounding cells. The selected cell being pinned has fixed n.

## 5. A positive compressed clock-pin cost in the literal model

Replace the pinned cell by its zero-bit u clock, tracing out its old
internal state and preserving the entire remaining marginal. This
removes (1), leaves every marker term (including both contact halves)
unchanged, and lowers the nonnegative readout. The new terminal energy
is exactly J tau_n <C_m>, at most J tau_n. The old terminal contribution
is nonnegative. Incident port expectations are zero before and after
replacement. Therefore, with F denoting exact tilted ground energies,

    F_clock-pin(s)-F(s) >= J d_clock,
    d_clock=d [(1-theta)^2 s^2]/[(1-theta)^2 s^2+16t^2]
                                                       - tau_n. (13)

The same interacting bulk cancels; this is not a comparison of separate
extensive perturbation errors. No charge link has been cut.

A simple sufficient schedule, with no new Hamiltonian term, is

    0<q<=1/4,  theta=1/2,  t/2<=s<=3t/4,
    nf<=1/64,  0<J/t<=1/80.                            (14)

Use the parent's mu_sel=-(n+2)E_n and n>=1. Then
d>=9/16, c<=25/16, |mu_sel|/2<=3/512, and r<5/3<7/4.
The right side of (11), divided by tr, is at least 9/700>1/80.
Equation (13) gives

    d_clock >= d/257-tau_n > 0.                        (15)

For n=1, d>=1 and tau_1<=1/272, giving
d_clock>=15/69904. For n>=2,
d_clock>=9/(16*257)-1/4096>15/69904.
At the illustrative moderate choice q=1/4 the n=1 bound improves to
2/4369; all n>=2 have a larger bound. This is a strictly positive
volume-independent clock-pin comparison in the **true terminal-penalized
model**, not only in the frustration-free reference.

There is no additional restriction on lambda for (13): incident ports
vanish under this vacuum-clock pin, and other ports are unchanged by
the channel. This statement does not imply that a large-lambda ground
state stays close to the vacuum code.

## 6. A ground-state bound without assuming the cell is vacuum

Extend Q_ff by zero on the charged space and define K_* in (7) on the
entire cell. The charged space is then unchanged by the channel.
Charged masses and charged moves still make no source error: all
nontrivial channel branches are in the vacuum space. The two incident
ports now contribute at most 4JB to the source energy change, where

    B=sqrt(2+2 eta)|lambda|,   r_B=r+2B.

For a full tilted ground state, (10) becomes

    sY <= 4t sqrt(XY)+2cJ r_B.                          (16)

A second channel leaves charged states alone, resets the vacuum sector
to the zero-bit u state, and discards vacuum/charged coherences. Its port
expectations vanish. The original two port expectations have magnitude
at most 2JB in total. Ground-state minimality and positivity of the old
terminal and readout terms give

    X <= tau_n+2B.                                     (17)

Young's inequality in (16) implies

    W <= (1+16t^2/s^2) X + 4cJ r_B/s.

Let P_vac be the cell vacuum flag and let P_0 include all payload states.
Since Q_ff>=d(P_vac-P_0), we obtain the volume-independent bound

    <P_vac-P_0> <= [(1+16t^2/s^2)(tau_n+2B)
                                      +4cJ r_B/s]/d.   (18)

It treats coherent clock gradients and contact/clock correlations
explicitly. It bounds vacuum-clock excitation weight; charged weight
must be controlled separately by the charged-rail result. It is not a
full many-body spectral gap or a quasi-local dressed-code construction.

Under (14), the previously sufficient frozen allowance

    B <= (1-q^2)^2 q^(2n+4)/10                          (19)

already makes (18) less than 1/2 for every fixed n. No smaller port is
needed for this coarse clock bound. For n=1 use c=d, tau_1<=1/272,
r<=587/512, B<=1/40960 and d>=1: the right side is below 0.36.
For n>=2 use d>=9/16, c<=25/16, tau_n<=1/4096, r<5/3 and
B<=1/655360: it is below 0.492. These rounded constants have slack.
At q=1/4, f=lambda=0, the sharper unrounded bounds at J/t=1/80 are
about 0.338 and 0.241 for n=1,2. They are deliberately conservative.

The exact off-diagonal bound (6), augmented by Jnf+2JB, also gives
<P_vac-P_0><=[(epsilon_n+nf+2B)/d_clock]^2 from the pinned block.
For n=1,2 at q=1/4 this last bound can be trivial; (18) is stronger.
We do not hide that failure by requiring a smaller q.

## 7. Ensembles, physical gluing and density

If V_* minimizes E(N,V)+sV, then the pin comparison at that same V_*
follows from (13), exactly as in the previous contribution. The
Hamiltonians compared have the same N,V and the same bath. A proof at
every preassigned density would need additional canonical control.

For physical cyclic gluing, use the exact uniform-flux representation
and the 48t/N ground-energy bound proved in the previous contribution.
A physical tilted ground state, expressed in its rooted flux block,
has energy at most zeta_N=96t/N above that block's unrestricted tilted
ground energy. The two variational channels above then give

    <P_vac-P_0> <= [(1+16t^2/s^2)(tau_n+2B+zeta_N/J)
                         +4cJ r_B/s+2c zeta_N/s]/d.     (20)

The cell-average is root independent. This is one global allowance,
vanishing at fixed J/t>0 as N tends to infinity; it is not an error per
clock or an O(t) charge cut. The strictly finite rooted bound (18) and
the physical limit statement (20) should not be conflated.

Density control need not restore a large-D condition. At zero port the
contact-weighted charged terms before marker and readout energies are
positive: their absolute-row incidences per head are at most
A_left+(1+|g|)A_right, whereas Delta/2>=2 supplies each contact.
The vacuum edges and terminal energies are also positive. Thus in any
legal pattern

    H >= H_ch-JN(|mu_sel|+B).

A zero-bit frustration-free product trial has source energy at most
JNq^2, independent of Delta,D. Comparing tilted free energies as in
the earlier proof, the continuum exterior-density margin exceeds t/16,
while the material comparison costs at most
J(q^2+|mu_sel|+B) per material site. Under (14),(19), this is below
t/1000. The finite free-sea comparison error is at most 10t/N^2.
Consequently the minimizing labelled densities are in (1/3,2/3) for
N>=32. Including the global gluing allowance gives the same conclusion
for physical tilted minima for N>=64. These claims require Delta>=4,
D>=0, |g|<=1/2, f>=0 and bounded rail hoppings, as in the family.

## 8. Scope and the remaining obstacle

This establishes a positive vacuum-clock compressed pin cost and a
ground-state clock-excitation estimate, without adding states, terms,
or a D-dependent smallness condition. Other interacting cells and the
gapless vacancy dynamics are retained. Positive readout and nonzero
ports have explicit budgets (14),(18),(19).

The theorem does not prove a dressed product code, a dynamical massive
influence kernel, or an all-order action whose error is below every
payload scale. The bare projector is not invariant after the selector
term is restored. Integrating out its controlled local excitation weight
while retaining arbitrarily slow payload modes is still a distinct
matching problem. A local pin inequality alone does not accomplish it.

Wrong-length marker patterns, illegal sectors, and the complete
interacting multichannel CFT remain separate. The sufficient J/t=1/80
is a clock-control allowance, not a certified high-c point. Moderate
q=1/4 is an illustrative all-length selector/clock choice; the parent's
arbitrary-n band proof chooses additional much smaller parameters.
No total nonchiral Virasoro central charge is established here.
