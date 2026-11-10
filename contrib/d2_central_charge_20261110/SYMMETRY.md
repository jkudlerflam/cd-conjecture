# Exact antiunitary momentum pairing

This is an algebraic property of the fixed Hamiltonian, not a central-charge
identification. Let L be even, P=product_j Z_j, K be complex conjugation
in the Z basis, and

    V = product_j exp[i (-1)^j pi Z_j / 4],     A = V K.

For Jx=-Jy and the z-directed DM interaction used here, A H A^-1=H.
This follows either by rotating each Pauli operator, or from the real
on-site gauge in REAL_GAUGE.md. The field and ZZ term are invariant.
Moreover A^2=1, V^2=P, and one-site translation obeys

    T A T^-1 = P A.

To see the last identity, translation reverses the alternating rotation
angle, so T V T^-1=V^-1=P V; T is real and commutes with K.
The phase in V^2 is one for every even L, leaving precisely product_j Z_j.

If a state has parity eigenvalue p and T eigenvalue exp(i k), then A sends
its translation eigenvalue to p exp(-i k). Hence energies are paired as
k <-> -k in the even sector and k <-> pi-k in the odd sector. In particular,
odd k=0 and k=pi levels are exactly degenerate. Their orthogonality follows
from their different translation eigenvalues.

This explains the odd-ring doublet found in the complete spectrum scan and
why an optimization within that degenerate space can change interval
entropy at essentially unchanged energy. It does not prove a gap, a
conformal tower, or a number of gapless fields. An arbitrary real-space
superposition of degenerate momentum states need not share the entropy
profile of either momentum eigenstate.

`validate_symmetry.py` checks V H* = H V, V^2=P, A^2=1, and
T V=p V T in both independent full-space L8 parity blocks. It also checks
the predicted momentum pairings in all completed L16,20,24 spectra.
