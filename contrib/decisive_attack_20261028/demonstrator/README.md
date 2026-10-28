# An arbitrarily weak XX reference inside the ten-letter demonstrator

## Claim and status

[proved] The second-order logical interaction of the existing ten-letter
nearest-neighbour, one-site invariant demonstrator can be made arbitrarily
close, in relative interaction strength, to the XX chain. This uses a fixed
nonzero microscopic ratio, chosen before the intercell coupling and before
taking the thermodynamic limit. The cell gap stays fixed. The precise finite
cell, parity and local reduction hypotheses of the existing finite-cell
massive-matching theorem hold in the conserved PBP sector.

[conditional on the repository's finite-cell matching theorem] That sector
therefore has an exact zero-temperature logical action equal to this weak
XXZ reference plus a small summable retarded correction. This does not prove
the subsequent RG theorem or identify the Virasoro central charge. In
particular, the source theorems examined below do not supply the full
conclusion requested in C47 by direct invocation.

This is an in-family simplification, not a c-d violation: d=10 and the
candidate infrared charge is c=1. It also does not prove selection of this
sector against all other type patterns.

Provenance: Codex GPT-6-family `demonstrator_ir` agent, 2026-10-28. The
algebraic proof below is checked by a separate exact rational resolvent
calculation in `verify_tuning.py`; no independent human review is claimed.
A separate campaign agent independently checked the resolvent formulas,
gap argument, parity and charge counts against the physical cell builder.
The original physical bond is in
[`routeA_demonstrator/cell.py`](../../routeA_demonstrator/cell.py).

## 1. Exact tunable effective bond

Set D0=D1=D>0, 0<|kappa|<D, and retain the physical port bond

    lambda [Jp (sigma+ sigma- + sigma- sigma+) + Jz Z Z].

All signs and normalizations here match the existing demonstrator. For one
active port, write p=Phi-, a=01 and b=10. The six excited states are

    p0, p1, a0, a1, b0, b1.

Their Hamiltonian has diagonal D and only two off-diagonal entries:
p0 <-> a1 with kappa, and p1 <-> b0 with -kappa. The other two states
a0 and b1 are singletons. The unused port-partner Bell state remains exactly
annihilated even when the logical spin changes. Thus a two-cell virtual
excursion is exactly the inverse of this six-state Hamiltonian's Kronecker
sum, acting on 36 excited states.

Since Z Phi+=p, sigma+ Phi+=a/sqrt(2), sigma- Phi+=b/sqrt(2), the source from
logical state c1,c2 is

    Jz |p c1,p c2> + (Jp/2)(|a c1,b c2>+|b c1,a c2>).

Put A=2D. The three diagonal resolvents encountered are

    r1 = 1/A,
    r2 = A/(A^2-kappa^2),
    r4 = (A^2-2 kappa^2)/[A(A^2-4 kappa^2)].

The last identity follows by diagonalizing the two independent 2-by-2
blocks: it is the mean of 1/(A+2 kappa), 1/A, 1/A, 1/(A-2 kappa).
The same-block two-flip entry has magnitude
kappa^2/[4D(D^2-kappa^2)]; the two block signs make it negative in the
opposite-logical-charge sector. The overall second-order minus sign then
gives

    H2 = a0 I + Jxy (sigma+ sigma- + sigma- sigma+) + Jzz Z Z,

    Jxy = Jp Jz kappa^2/[4D(D^2-kappa^2)],
    Jzz = 3 Jp^2 kappa^4/[32D(D^2-kappa^2)(4D^2-kappa^2)],
    a0 = -Jz^2 r4 - Jp^2(2r2+r1+r4)/8.

For completeness, the diagonal energies from the transverse source are
-Jp^2 r2/2 for equal logical spins and -Jp^2(r1+r4)/4 for opposite spins.
Their half-difference is Jzz because

    Jzz = Jp^2(r1+r4-2r2)/8.

The longitudinal source gives the same -Jz^2 r4 for all four logical states.
Mixed sources produce exactly the displayed hopping and no diagonal term.
These observations establish every entry of H2 and the formula for a0.

On a chain, second-order contractions on two distinct bonds vanish. An
outer cell is touched by only one port insertion; its other Bell pair stays
in the code, while the touched pair lies in the invariant excited
six-state space. No local resolvent restores that outer cell to the code.
The second-order chain is therefore the sum of these pair bonds.

After the usual staggered rotation if the hopping sign requires it, the
absolute XXZ anisotropy is

    |Delta_XXZ| = 2 Jzz/|Jxy|
                = 3 |Jp/Jz| kappa^2/[4(4D^2-kappa^2)].

Fix D, kappa and Jz nonzero. For every prescribed epsilon>0, choosing a
fixed nonzero Jp with

    |Jp/Jz| < 4 epsilon(4D^2-kappa^2)/(3 kappa^2)

puts the target inside the epsilon-neighbourhood of XX. This removes the
unjustified assumption that the original numerical value 0.14835 lies
inside an unspecified small-coupling radius of a rigorous RG theorem.
It does not assert such an RG theorem for the final retarded action.

## 2. The finite-cell matching assumptions actually hold

[proved] In the frozen PBP sector the one-cell Hilbert space has dimension
32. Its five qubits are pL,sL,c,sR,pR, and its code consists of the two Bell
pairs and an arbitrary c. Let QL and QR project off the corresponding Bell
states, including identity on c. The checks obey

    hL >= (D-|kappa|) QL,  hR >= (D-|kappa|) QR.

QL and QR commute, since they act on distinct port-partner pairs. Hence

    Hcell >= Delta(1-Pcode),  Delta=D-|kappa|>0,
    ker Hcell = ran Pcode.

The parity Gamma=ZpL ZsL Zc ZsR ZpR commutes with Hcell and restricts to
logical Zc. A transverse port operator is odd and Zport is even, so the
intercell bond is even. The Jordan-Wigner identification preserves
interval supports; its parity-dependent ring boundary condition must be
retained. Both parity sectors have dimension 16 and one code state, leaving
15 even and 15 odd excited states. The spectral embedding used in
[`FINITE_CELL_MATCHING.md`](../../direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md)
therefore uses 45 auxiliary modes and one logical mode, with auxiliary
one-particle energies at least Delta/2. These are analytic auxiliary
variables and do not alter microscopic d=10.

The exact charge is (ZpL-ZsL+Zc-ZsR+ZpR)/2. Its code charges are +/-1/2.
Taking code vacuum charge +1/2 gives a logical fermion charge -1; excited
states have integral charge differences from that vacuum, with parity equal
to that difference modulo two. Thus the charge-preserving version of the
spectral embedding applies as well. Auxiliary integration leaves exact
logical U(1). Charge conjugation and reflection can be transported through
the change of basis, although their action on auxiliary Fock coordinates
need not be linear.

The explicit finite-order local SW construction in
[`C45 PROOF`, Section 3](../../../codex/work_C45/PROOF.md) applies with
Fpert=0 and PVP=0. For an exponential support weight a, a sufficient input
bound is

    v0 = max(1, 2 exp(2a)(|Jp|+|Jz|)),
    G = 48 v0^2/Delta,
    C_SW = 2^16(v0 G/Delta + v0^3/Delta^2).

When lambda <= min(v0/G,Delta/(192v0)), its exact conjugation gives a locally
block-diagonal reference plus a remainder bounded by C_SW lambda^3 in the
weighted interaction norm. The compressed reference is precisely lambda^2
sum H2. A scalar density is retained in this equality.

Choose the field weight and the massive time-decay weight first. With
F=46 fermionic modes per cell, choose the physical support weight strictly
larger than a_out+46 log(8(H')^2), as required by the matching theorem's
coefficient conversion, and recompute C_SW at that weight. These constants
are very large but finite. The auxiliary block interaction has coefficient
norm O(lambda), and the remainder has norm O(lambda^3). The theorem's
smallness condition is then met by taking lambda sufficiently small.

[conditional on that matching theorem] For each fixed nonzero Jp, the
zero-temperature action correction has norm at most C lambda^3, including
the stated spatial and time moments. Relative to the hopping bandwidth

    W=|Jxy| lambda^2,

the correction is at most C lambda/|Jxy|. Thus the correct order is: choose
the desired small XXZ interaction through Jp/Jz; fix all finite coefficient
conversion constants; choose a fixed positive lambda sufficiently small
relative to |Jxy|; only then increase the number of cells. Setting Jp=0
would remove the reference hopping and is not an allowed substitute.

This verifies the matching hypotheses in one conserved sector. It does not
claim that a finite-order SW unitary exactly separates the whole low-energy
spectrum. The output is an action with summable retarded kernels.

## 3. Published results that still cannot be invoked verbatim

[read] Benfatto and Mastropietro,
[arXiv:0907.2837](https://arxiv.org/pdf/0907.2837), equation (1) and Theorem
1.1, treat nearest-neighbour XY exchange plus an exponentially decaying
instantaneous S3-S3 interaction. The theorem establishes relations among
critical exponents and the density susceptibility. It does not state a
finite-temperature pressure or Casimir theorem for arbitrary retarded,
higher-degree fermionic kernels.

[read] Mastropietro,
[arXiv:1303.5977](https://arxiv.org/pdf/1303.5977), equations (1)-(3) and the
theorem on page 2, treats a specified next-nearest-neighbour perturbation
with both the XXZ interaction and perturbation sufficiently small. Its
conclusion concerns zero-temperature Drude weight and susceptibility.

[unproved] Extending the single-channel RG induction to the matched
high-degree retarded action, its transported sources, and its scalar
finite-temperature term remains necessary. The finite-cell theorem's
zero-temperature cancellation of the block-diagonal auxiliary baseline
does not bound all finite-temperature source kernels. Its separate
exponentially small pressure comparison for the reference is not a
comparison between the complete interacting microscopic pressure and the
final logical RG pressure.

The present simplification makes that extension a weak-XX problem with a
finite microscopic cell, eliminating the unnecessary issue of controlling
a fixed anisotropy near 0.15. It does not by itself close C47 or the c-d
conjecture.

## Reproduction

Run `python3 -B verify_tuning.py` in this folder. It checks all entries of the
36-by-36 virtual resolvent against the closed formula in twelve rational
parameter cases and writes `tuning_checks.json`. Only the Python standard
library is needed. A separate floating-point comparison with the existing
1024-dimensional physical two-cell construction at kappa=0.2,0.4,0.6
agreed to roundoff; that comparison is supporting evidence, not the proof.
