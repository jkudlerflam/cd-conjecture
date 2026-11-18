# Nonzero-port ground-sector selection for the modified frozen d=9 bond

2026-11-18. **[Proved below; bounded literal corroboration; no independent
review.]** This transfers the bulk-subtracted cutting method of
[small_d_defect_stability_20261116](../small_d_defect_stability_20261116/PROOF.md)
after checking the actual two-site boundaries, zero-data cells and L1
readout. It uses no effective-band or conformal assumption.

## 1. Statement

Fix n>=1, b=n+1 and the literal modified nine-state Hamiltonian. Assume
q<=1/4, w<=1/2, Delta=3, D>=0, |g|,|t_d|<=1/4 and h_L>=0. Keep all
phases fixed, including on wrong-length cells. Use

    nu_head=1+q^2-(n+2)E_n,
    gamma=E_(n+1), B=sqrt(2)|lambda|,
    K>=0,

where E_m is the smallest eigenvalue of the clock Q_m in the prior
all-n derivation, with E_0=q^2. If

    gamma>5B,                                       (1)

then on every physical ring L=M b>=2 the full microscopic ground space
lies in the union of the b translated sectors consisting entirely of
one-head length-b cells. If S_M is their common exact ground energy,
every other sector with N_W nonselected boundary cells obeys

    E_sector-S_M >= N_W(gamma-5B)>0.                (2)

The only sectors without any boundaries have energy >=0>S_M and are
excluded separately. A nonselected cell includes both a live cell of the
wrong length and an inert headless cell. Different total head numbers
are therefore included. K=0 is allowed: no penalty gap or 1/K expansion
is used. The gap in (2) is a *sector* gap, not a logical excitation gap.

## 2. Why two-site boundaries give complete factors

For a word define

    B_i=[letter_i is not L][letter_(i+1) is not R].

Every B_i and every LR forbidden-pair projector is conserved by every
literal move. Fix their pattern. Between two successive B bonds, every
interior bond has its left letter L or its right letter R. The first
letter is not R and the last is not L. Exactly two cases are possible:

- A **live cell** L* H R* with one of the five heads. For m data sites
  it has the complete 5(m+1)2^m state space, including m=0.
- A **dead cell** L*R* with both runs nonempty, no head and precisely one
  LR transition. Its length is at least two. The LR position can be
  treated as a further conserved label.

There cannot be two heads in one such interval: the first head forces
only R letters to its right until a B boundary. A headless interval must
begin with L and end with R. This proves exhaustion, including the
formerly illegal sectors. If there are no B bonds on the ring, every
letter is L or every letter is R; these words are inert.

All terms except ports act within the complete live-cell factors. At
the left end the preceding letter is always in NOT_L, so the switch and
the W bias are identities on that preceding letter. At the right end
the following letter is always in NOT_R, so the closing link is likewise
intracell. This remains true when m=0 and either neighbour is a head.
There are no negative end-fragment energies to estimate.

Dead cells have only nonnegative readout and K energy, and no port can
touch them: their endpoint letters are data letters. A port between two
live cells has norm sqrt(2)|lambda| for every pair of lengths, including
two zero-data cells. The proof uses its full microscopic norm, not its
small code visibility.

## 3. All isolated-cell energies and costs

At h_L>=0, the full isolated live cell has exact minimum

    nu_head+e_v(m)=E_m-(n+2)E_n.                    (3)

The readout is positive and the all-zero-bit W path attains this energy.
The charged block has row-sum lower bound
nu_head+3-2-|t_d|-|g|, which lies above the W minimum for every m, including
m=0. At m=0 both end switches may join the same two head states, but
their combined absolute amplitude is at most 1+|t_d|; the same bound holds.

The isolated selected density is e_0=-E_n. A live m-cell has excess

    delta_m=E_m+(m-n-1)E_n.

The exact clock estimates

    (1-q^2)^2 q^(2m+2)<=E_m<=q^(2m+2),
    E_(m+1)/E_m<=r=q^2/(1-q^2)^2<=16/225

give delta_n=0 and inf_(m!=n) delta_m=E_(n+1)=gamma. For shorter cells,
r^(-k)-(k+1)>1 with k=n-m; for m=n+1 the cost is exactly gamma; the
remaining longer cells cost at least E_n. This includes m=0 exactly,
without an invented zero-data clock endpoint.

A dead cell of length ell has excess at least ell E_n+K>=2E_n>gamma.
Thus every nonselected boundary cell has a positive isolated excess
at least gamma, including sectors with altered head number.

## 4. Subtract the exact dressed bulk

Let A_j be the ground energy of j selected cells with open port boundaries,
retaining every internal head, clock, field and conversion state. Joining
two runs adds one port of norm B, so

    |A_(j+k)-A_j-A_k|<=B.

Subadditivity and superadditivity give a common limit u=lim A_j/j and

    |A_j-j u|<=B.

Define e=u/b. Closing a run changes its energy by at most B, giving

    |S_M-eL|<=2B,
    e<=e_0=-E_n<0.                                  (4)

The second inequality uses the all-zero vacuum product: every port has
zero expectation and the readout vanishes. For M=1 the port on the legal
one-cell space is zero; the same estimate holds. For M=2 the two oriented
ports are both retained. No extensive perturbative error is compared with
a finite defect energy.

In a competing sector delete all ports incident to a nonselected cell.
There are at most 2N_W deleted ports and at most N_W nonempty remaining
selected runs. Their boundary remainders cost at most another B N_W.
For an isolated wrong live cell use (3) and e<=e_0; for a dead cell use
its nonnegative minimum and e<=e_0. This yields

    E_sector-eL >= sum_wrong delta_cell-3B N_W.       (5)

Subtracting S_M<=eL+2B and using N_W>=1 proves (2). This is the step that
excludes a finite number of defects on every compatible finite ring; a
positive bulk density cost by itself would not suffice.

## 5. Schedule, checks and scope

An explicit sufficient choice is

    gamma_low=(1-q^2)^2 q^(2n+4),
    0<lambda<=gamma_low/(10sqrt(2)).                 (6)

Then every nonselected cell costs at least gamma_low/2. Lambda can be
decreased further for the local matching and velocity conditions. The
readout h_L=lambda^2 chi is nonnegative and automatically allowed; no
retuning of nu_head or extra scalar counterterm is needed. The exact
bulk e is an energy reference, not a new physical coefficient.

`checks.py` independently assembles all physical words on L=2 and L=3
rings (81 and 729 states), including headless and zero-data competitors.
It checks conservation of both boundary and forbidden-pair labels and
tests (2) in every competing block at a nonzero port satisfying (6).
The field .02 in this selector diagnostic is deliberately larger than
the band readout: the proof permits any nonnegative field. This numerical
check corroborates the exact argument; it is not an interval certificate.

This closes the **nonzero-port frozen ground-sector comparison for the
modified d=9 family on compatible rings**. It does not close interacting
CFT matching. Rings with b not dividing L cannot have the selected tiling;
the existing bounded-defect and growing translation-orbit issue of frozen
models persists. Those defects are not an extra conformal mode. No claim
for the original contact-weighted d=10 mobile model is made: its moving
vacancies destroy this fixed-factor cutting setup.
