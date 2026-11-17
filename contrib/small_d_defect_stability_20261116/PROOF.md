# Exact defect stability of the frozen 20-state chain

Date: 2026-11-16. **[Proved, not independently refereed.]** This is an
application and rederivation of the cutting argument in
[route A §§2–5](../counterexample_campaign/route_a/README.md), with the
new selector, negative marker energy, and readout checked explicitly.
No perturbative effective Hamiltonian or infrared assumption is used.

## 1. Statement and hypotheses

Use precisely the material bond of
[the 20/22-state construction §1](../smaller_d_unbounded_20261116/PROOF.md).
There are 20 states, nearest-neighbour range, a Hermitian bond repeated
with one-site translation symmetry, and size-independent couplings.
Fix an integer n>=1 and put b=n+1. Assume

    0<q<=1/4,  0<w<=1,  0<v<=1,  eta>0,
    D>=0,  |g|<=1/2,  Delta=4+E_n,
    mu=-(n+2)E_n,  f>=0.

Here E_k is the smallest eigenvalue of the (k+1)-site vacuum matrix Q_k
defined in §3 below. The phase pi/(2n) and all other bond conventions
are unchanged. Lambda is real and nonzero. Define

    gamma = E_(n+1),
    B = sqrt(2+2 eta) |lambda|,
    d_W = gamma-5B,
    d_I = K-2|mu|-7B.

**Theorem.** If d_W>0 and d_I>0, on every compatible physical ring
L=mb>=2 the full ground space is contained in the union of the b
translated all-length-n marker sectors. For m>=1 let S_m be the exact
ground energy of one such sector, with the entire cell Hilbert spaces
and all port and readout terms retained. Every legal marker sector with
N_W>=1 wrong cells satisfies

    E_sector - S_m >= N_W d_W > 0.                      (1)

Every sector with N_I>=1 forbidden bonds, and N_W complete wrong cells
between the cuts, satisfies

    E_sector - S_m >= N_W d_W + N_I d_I > 0.            (2)

Legal marker-free sectors have nonnegative energy, while
S_m<=-L E_n<0. No payload many-body gap is assumed or concluded. The
constants are independent of L, the number and arrangement of defects,
and the lengths of wrong cells or boundary fragments.

There is also a common, *exact dressed* bulk energy density e<0, defined
in §5, for which the stronger bulk-subtracted inequalities are

    E_legal - e L >= sum_wrong delta_k - 3B N_W,         (3)
    E_illegal - e L >= sum_wrong delta_k - 3B N_W
                             + (K-2|mu|-5B) N_I.       (4)

Here delta_k=E_k+(k-n-1)E_n; no approximate extensive energy is subtracted.

## 2. Exact conserved sectors and cell decomposition

The arrow alphabet is independently checked as follows:

| Letters | Incoming, outgoing arrow | Marker? |
| --- | --- | --- |
| L_a | 0,0 | no |
| R_a | 1,1 | no |
| V,D_(r,s),X_(r,s) | 0,1 | no |
| E_a | 1,0 | yes |
| M_v,M_(r,s) | 0,0 | yes |

For every off-diagonal rule, the two site marker indicators are
individually unchanged, the pair is legal both before and after the move,
and the two exterior arrows are unchanged. These checks imply that each
site's marker projector and each bond's *forbidden-match* projector
commutes with H. They do not assert conservation of each site's individual
incoming or outgoing arrow. The latter can change in a head move.

Every marker has outgoing arrow 0. Between consecutive markers a legal
word therefore has L's, at most one 0-to-1 head, and then R's ending at E;
or all L's ending at M. The displaced bit on E gives precisely
(9k+5)2^k states for a k-bit cell, including the five states at k=0.
The fixed-marker legal Hilbert space is the tensor product of these
complete cell spaces. No states are discarded.

All terms except the ports act within these factors. In particular the
boundary junction m D <-> m X is identity on the marker m and belongs to
the cell on its right; it is not an additional intercell coupling. A port
links the single source |M_v,V> to four mutually orthogonal targets,
with squared amplitudes 1,1,eta,eta. Its exact two-site norm is therefore
sqrt(2+2 eta). Compression to cells or fragments cannot increase it.

At a forbidden bond all non-diagonal terms and all vacuum edge
projectors vanish. Its K penalty is a scalar in the conserved sector.
After assigning each site's onsite energy once, cutting these bonds is
exact; no finite-K virtual excursion or 1/K correction exists.

## 3. Independent selector derivation, including all lengths

For k>=0,

    Q_k = sum_(j=0)^(k-1) (q|j>-|j+1>)(q<j|-<j+1|)
                                             + q^2 |k><k|.

Let B_k be upper bidiagonal with diagonal q and superdiagonal -1.
Then Q_k=B_k^* B_k, and (B_k^-1)_ij=q^(-(j-i+1)) for j>=i.
Consequently

    ||B_k^-1||^2 <= ||B_k^-1||_F^2
       = q^(-2k-2) sum_(s=0)^k (s+1)q^(2s)
       <= q^(-2k-2)/(1-q^2)^2.

The normalized vector proportional to (1,q,...,q^k) supplies the upper
bound. Thus, with r=q^2/(1-q^2)^2<=16/225,

    (1-q^2)^2 q^(2k+2) <= E_k <= q^(2k+2),
    E_(k+1) <= r E_k.                                   (5)

For the chosen mu the zero-port energy per site at length n is
e_0=(mu+E_n)/b=-E_n. The excess of any other cell is

    delta_k = mu+E_k-(k+1)e_0 = E_k+(k-n-1)E_n.

At k=n+1 this is E_(n+1). For k>=n+2 it is at least E_n, and
for k=n-m<n it is at least [r^(-m)-(m+1)]E_n>E_n.
The last inequality holds for every integer m>=1 since r^(-1)>=225/16
and successive increments of r^(-m)-(m+1) are positive. Hence

    delta_n=0,   inf_(k!=n) delta_k=E_(n+1)=gamma>0.       (6)

This proves the infinite-tail statement, not a cutoff test.
Deleting the first row and column of Q_k gives a Toeplitz matrix bounded
below by (1-q)^2. Interlacing and E_k<=q^2 imply an excitation gap at
least 1-2q>=1/2 when k>=1. At k=0 there is no excited vacuum level.

## 4. Cells and fragments with nonzero readout

At zero port but arbitrary f>=0, a complete cell has *exactly* ground
energy mu+E_k. The vacuum operator before the readout is Q_k tensor the
payload identity. The readout f N_(L_1) is positive and annihilates the
entire all-zero-bit vacuum path. It cannot raise this minimum.

An excited cell has one excited head (including terminal heads). Its
kinetic graph has absolute row sum at most 2: D hopping has modulus at
most 1, the express edges have modulus 1, and the junction joins the
endpoints of these two rails. Each basis state participates in at most
one conversion of magnitude |g|. Thus, before the marker potential,
the excited block is bounded below by

    Delta-2-|g| >= 3/2 > q^2 >= E_k.

This remains true for arbitrarily long cells. It also remains true for
truncated end fragments: restricting available head moves cannot increase
the absolute row sum. Vacuum terms are positive rank-one operators,
the terminal term and readout are positive, and a marker-free fragment
with no head has nonnegative energy.

Cutting N_I forbidden bonds leaves N_I legal intervals. Decompose each
into complete cells and at most two end fragments, choosing the prefix
to end at its first marker and the suffix to start after its last marker.
An interval without markers is one fragment. Each fragment has at most
one marker and at most one head, so its energy is at least -|mu|,
independently of its length. In fact this convention permits a sharper
one-negative-marker bound per interval; we only use the conservative
2|mu| per cut. This is the adjustment needed for this construction's
negative marker potential.

## 5. Dressed bulk energy without a payload gap

Let A_j be the ground energy of j complete selected cells with open
*port* boundaries. The junction at the left edge is retained as the
one-cell operator described in §2. This matches a run cut out of a ring.
Concatenating two runs changes their Hamiltonian by one port of norm
at most B. The variational principle gives

    |A_(j+l)-A_j-A_l|<=B.

The subadditive sequence A_j+B and superadditive sequence A_j-B have
the same energy-density limit u=lim A_j/j: their difference divided by j
vanishes. Equivalently, repeating a fixed run and taking its number of
copies to infinity gives the useful uniform estimate

    |A_j-j u|<=B,   j>=1.                               (7)

Set e=u/b. Closing a run changes its energy by at most B, so

    |S_m-e L|<=2B.                                     (8)

For m=1 the microscopic closing port is zero on the legal one-cell
space (a cell cannot contain both M_v and V); S_1=A_1 and (8) still holds.
For m=2 both oriented ports are retained, as in the physical ring.
The all-zero-bit product of vacuum ground paths has zero expectation
of every port and readout, whence

    u<=mu+E_n,   e<=e_0=-E_n<0,   S_m<=e_0 L.          (9)

The exact e includes all scalar dressing, readout effects, payload
interactions and auxiliary admixtures. Equations (7)–(9) use neither
an expansion in lambda nor a uniform spectral gap of the payload.

## 6. Surgery and closure of the boundary remainder

On a legal ring with N_W>=1 wrong cells, delete every port incident
to a wrong cell. There are at most 2N_W such ports and at most N_W
remaining nonempty selected runs. Their removal lowers the ground-energy
lower bound by at most 2B N_W. Apply (7) to the selected runs and §4 to
each isolated wrong cell. Since e<=e_0,

    mu+E_k >= e(k+1)+delta_k.

Adding all lengths proves (3). Comparing to S_m<=eL+2B on a compatible
ring gives

    E_sector-S_m >= sum_wrong delta_k -3B N_W -2B
                 >= sum_wrong (delta_k-5B),            (10)

because N_W>=1. This proves (1). Keeping only (3) with
gamma>3B would establish a positive cost relative to the bulk, but would
not by itself exclude a finite number of defects on every finite ring.
The stronger gamma>5B absorbs the entire boundary remainder. It rules
out every competing defect sector, including a single wrong cell.

For N_I>=1 cuts there are at most N_W+N_I selected runs, at most 2N_I
end fragments, and at most 2N_W+4N_I ports incident to wrong cells or
fragments. Use (7), the cell bound, the fragment lower bound -|mu|,
and e<0 to replace a fragment's zero baseline by e times its length.
Add K N_I. This proves (4). Subtracting (8), and using N_W+N_I>=1,
gives (2). These counts are conservative and also cover sectors with no
complete cells, and rings in which every bond is forbidden.

Marker-free legal rings are all L or all R: a 0-to-1 head could not
close its arrow circuit without a 1-to-0 marker. Their minimum is zero
for f>=0. Equation (9) excludes them directly, with no boundary loss.

## 7. Explicit finite schedule and compatibility with the band hierarchy

Define the closed-form lower bound

    gamma_low=(1-q^2)^2 q^(2n+4).

After fixing the parent construction's n,w,eta,u,v,nu,D,g, with
q=(wu)^2<=1/4, choose *positive* nu and any

    0<lambda<=gamma_low/[10 sqrt(2+2 eta)],
    f=lambda^2 g^2 |v_0 v_n|^2 D^(-2n-2) nu/q^2,
    K=1.                                               (11)

Here D>0 and g!=0 are chosen as in the leading-band argument, and v_j
in the field formula denotes the vacuum amplitude, not the bus ratio v.
The theorem allows any f>=0, so this required small positive readout is
automatically covered. Positive nu is allowed by the parent's finite
root-separation argument. If the finite-band approximation asks for a
smaller lambda, take the minimum of that allowance and (11); the defect
inequalities only improve. No ordering of the preceding parameters is
reversed and no lower bound on a payload velocity is introduced.

Indeed B<=gamma_low/10 and

    d_W>=gamma_low/2,
    2|mu|<=2(n+2)q^(2n+2)<=3/128,
    7B<=7/(10*4096),
    d_I>=1-3/128-7/40960>0.976.

All numbers are finite and strictly positive at each fixed n. This is
an explicit sufficient defect schedule; the parent's finite-band and
infrared conditions are additional. A smaller sequence lambda_s equal
to the displayed allowance times 2^(-s), s=0,1,..., leaves arbitrary
room for its later finite-order tolerances. No chemical-potential
retuning or extra scalar counterterm is required: mu remains exactly
-(n+2)E_n, and e is an energy reference, not a new Hamiltonian term.

## 8. Scope, incompatible rings, and remaining work

The bulk-subtracted estimates (3)–(4) hold on every ring; the comparison
to an all-n sector requires b|L. Frozen marker positions cannot supply
an all-n tiling when b does not divide L. The theorem therefore does
not repair the
[incompatible-ring obstruction](../decisive_attack_20261028/frozen/ALL_RING_OBSTRUCTION.md).
The same bounded-defect / translation-orbit argument applies with
delta_*=min(gamma-3B,K-2|mu|-5B)>0: inserting fewer than b sites into a
compatible selected trial costs at most a constant depending on b and
the fixed bond norm, so incompatible ground sectors have a bounded,
nonzero defect count. Their translation orbit, and hence ground-space
multiplicity, grows at least linearly in L. This concerns the all-ring
finite-size CFT interpretation, not a new propagating boson.

More explicitly, write L=mb+s, 0<s<b, and take m>=2. Inserting the s
sites changes at most s+2 bonds, so the compatible trial gives
E_0(L)<=eL+C with C=2B+2b||h||. Here |e|<=||h||. For sufficiently
large L this is negative and excludes the legal marker-free sectors.
Equations (3)–(4) then give 1<=N_W+N_I<=C/delta_*. If the conserved
pattern's translation orbit has size p, its defect count is L/p times
a positive integer. Thus p>=L delta_*/C; the translated sectors are
orthogonal and exactly degenerate.

The complete microscopic ground-energy density equals e. For example,
(3)–(4) and the marker-free bound give the lower bound eL-O(B), and
compatible selected sectors achieve e; existence and equality on all
length sequences follow by the same bounded-bond insertion argument.
An infinite-volume state with finitely many defects can share this
energy density. The positive finite-ring sector comparison is stronger
than a density statement, but is not a classification of all algebraic
infinite-volume ground states with imposed boundary conditions.

This settles the nonzero-port ground-sector and dressed-defect input
for the **frozen d=20 precursor**, at the sufficient schedule above.
It does not give the selected interacting chain's conformal limit.
In particular, c=n is still a conditional total nonchiral Virasoro
central charge, with a compact boson counted as one.

For **d=22**, vacancy motion makes material boundaries dynamical and
contact-weighted; the exact frozen cutting comparison is not a theorem
for that model. Fast-vacancy dynamical matching, positive dressed mobile
defects at fixed nonzero J/t, and the full multichannel conformal
identification remain separate. The
[2026-11-09 route audit](../route_audit_20261109/README.md) keeps these
model-specific gaps explicit. Older stronger headlines in STATUS.md
must not be read as proofs of those hypotheses for this new family.
No unbounded-c counterexample is established by this result.
