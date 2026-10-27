# Two cell lengths with identical ports and coherent payload transport

Status: proved finite-dimensional code construction. No local parent or
finite-coupling CFT is asserted.

## 1. Padding an arbitrary code

Suppose V:C^D -> (C^q)^(tensor k) is an isometry, k>=2, whose first and last
site compressions are

    V^dagger F_a^(first) V = A_a,
    V^dagger F_a^(last) V = B_a.                          (1)

No commutation relation between A and B is required. Add one new local
vacancy symbol |0> orthogonal to C^q, so the new site dimension is q+1.
Extend F_a by zero on that vacancy. Regard V as living on k occupied sites.
For p=1,2 define

    W_p = (V tensor |0>^p + |0>^p tensor V)/sqrt(2).       (2)

The two summands have different, orthogonal occupation patterns. Hence
W_p^dagger W_p=I. In the first summand only the left physical endpoint is
occupied; in the second only the right endpoint is occupied. A one-site
operator cannot connect the two occupation patterns. It follows that

    W_p^dagger F_a^(first) W_p = A_a/2,
    W_p^dagger F_a^(last) W_p = B_a/2,    p=1,2.          (3)

Both lengths k+1 and k+2 therefore have the SAME boundary maps. This is
stronger than padding a single code on just one side, which would erase
one of its ports.

## 2. Moving the boundary preserves both logical states after projection

Place a short cell followed by a long cell on 2k+3 physical sites. Move the
boundary one site to the right and leave the spectator variables untouched.
The old and new isometries are W_1 tensor W_2 and W_2 tensor W_1. Then

    (W_2 tensor W_1)^dagger (W_1 tensor W_2)
         = (1/4) I_(D tensor D).                         (4)

To prove this, expand each side into its four occupation patterns. Label
the first site zero. The two cores on the input start at

    (0 or 1, k+1 or k+3),

and the two output cores start at

    (0 or 2, k+2 or k+3).

The only common occupied-site pattern has core starts (0,k+3). In it the
site transferred across the boundary is vacant. The overlap of each
unchanged core is V^dagger V=I, and the four factors 1/sqrt(2) give 1/4.
All other terms vanish by orthogonality of vacancy and occupied symbols.
The same statement holds for a hop in the opposite direction by taking
the adjoint. Adding a projector onto a vacancy at the transferred site
does not change (4).

This is a projected transition, not an exact invariant encoding. Its
amplitude is less than one; a physical hop also has components outside
the desired code spaces. It does not evade the earlier exact-product-code
theorem, whose hypotheses require that leakage to vanish.

## 3. Leading interaction is independent of cell length

For a boundary interaction sum_a F_a tensor F_a, every pair p,p' in {1,2}
has the same compression

    (W_p tensor W_p')^dagger h_boundary (W_p tensor W_p')
         = (1/4) sum_a B_a tensor A_a.                   (5)

If A_a=alpha Gamma_a and B_a=alpha i Gamma_a P, this is the previous
quadratic Majorana target with hopping coefficient lambda alpha^2/4.
If selector clock overlaps contribute a scalar A_b^2 as in the mobile
selector, the marker-hopping amplitude acquires the additional factor
1/4 from (4). The two length variables carry XX hopping tensor the
identity on the payload. The payload interaction (5) is independent
of both length variables.

Thus in each fixed-cell-number sector, and before higher-order corrections,
the bulk effective Hamiltonian is a sum of an XX length chain and the
logical target on the same cyclic cell index. On a physical ring their
translation and twist sectors are correlated by cyclic relabeling of the
cells. They are not unrestricted tensor products of finite-size spectra.
The local stress tensors would still add if this leading description
were the actual stable infrared theory. That last premise is not proved.

There is no local-parent theorem hiding in (2). The projectors
I-W_p W_p^dagger are arbitrary (k+p)-site operators. Even if V had a
nearest-neighbor parent, making the coherent two-placement ground spaces
in (2) is another problem. Hopping vacancies through the core introduces
additional positions and modifies its boundary moments. The construction
only removes a linear-algebra incompatibility between two lengths.

## 4. Dimension bookkeeping and a check

One may tensor this new spectator with the four-letter selector, at
d=4(q+1). For example the graded construction in the earlier C8 report
allows q=17 and m=2*8*9=144 odd channels. Its logical dimension is 2^72,
its formal target charge is m/2=72, and its Choi dilation fits in

    k = ceil(log_q(q^4 2^72)) = 22

spectator sites. The padding gives lengths 23 and 24 and d=72. The formal
payload alone would exceed d-1=71; a length mode would add one more charge
if it survived. These numbers establish no microscopic example. In
particular a uniform two-site parent for either padded code has not been
given, and the literal selector length counts must be aligned with the
chosen placement of its markers before constructing such a parent.

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 padding_transport.py`.
It uses the earlier q=2,m=2 simultaneous Choi code, with its environment
placed between its endpoints, as a k=5 core. It checks both padded
isometries, both sets of endpoint compressions, the full two-code overlap
(4), and all four interactions (5). These are matrix identities on the
physical spectator spaces, not tests of the formulas substituted into
themselves. The output is `padding_transport_checks.json`.
