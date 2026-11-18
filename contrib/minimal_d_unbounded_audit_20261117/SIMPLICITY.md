# Physical complexity and tuning audit

This is an eleven-state construction with demanding parameters, not a
claim that eleven states alone makes a model simple. The convention is
one repeated microscopic bond on the stated physical alphabet. No
blocked or encoded dimension replaces d=11.

## Independent real coefficients

Fix target n and the phases stated in PROOF.md. Express and junction
hops set the common positive energy scale to one. The boundary bias is
then fixed to one and Delta to three. The general modified bond has
nine remaining real rule coefficients:

    q, w, D, K, mu_total, t_d, g, lambda, h_L.

These coefficients are independent before a selector or band schedule
is imposed. The tree sets t_d=0. Selection fixes mu_total as a function
of q,n, and the band schedule trades h_L for the freely chosen nu.
This leaves seven freely chosen real numbers in the stated tree
schedule: q,w,D,K,g,lambda,nu. Inequalities restrict their region.
The original C84 rule family at fixed phases has eight independent
ratios: Delta,D,K,mu_total,t_d,g,lambda,h, after its common hop scale.
Its certified schedule prescribes those numbers more tightly.

The counts above are not just names for grouped formulas. For each
rule operator h, checks.py removes its constant and telescoping terms:

    h_can=h-c I-(A tensor I-I tensor A),
    c=tr(h)/d^2,
    A=[tr_right(h)-tr_left(h)]/(2d).

The ten modified rule operators, including the fixed unit-scale part,
have rank ten after this quotient. The nine original operators have
rank nine. This finite matrix rank check is recorded in checks.json.
The phases are prescribed n-dependent constants, not counted as free
phase controls. Allowing arbitrary phases would add controls.

## Literal monomials in the physical matrix-unit basis

These are raw counts for the bond convention that splits onsite
potentials evenly between its two ends. They are not counts of free
parameters, and they are not invariant under changing the operator basis.

| Transition family | Undirected transitions |
| --- | ---: |
| Vacuum hopping, including the final E edge | 3 |
| Charged data hopping | 4 |
| Charged express hopping | 4 |
| Express-to-marker edges | 2 |
| Data-to-marker closing edges | 2 |
| Marker-conditioned data/express switches | 6 |
| Charge conversion | 1 |
| Pair ports | 2 |
| **Total with a nonzero loop** | **24** |

Each transition and its Hermitian conjugate are two directed complex
matrix units. At generic nonzero coefficients the modified loop bond
has 48 directed off-diagonal and 115 diagonal units: 163 entries.
The tree omits the two closing transitions, giving 44+115=159.
The old generic C84 bond has 166 entries in the same convention.
Removing the R1 readout changes this raw count; it does not mean that
the new boundary interaction is free. That interaction changes three
already occupied diagonal entries. Grouping them as one projector sum
does not remove their physical content.

The counts are before the quotient above, since that canonicalization
can spread an operator over more matrix units. No lower raw entry count
is presented as an invariant simplicity improvement.

## Hierarchy, precision and robustness

For a coefficient comparison use the independent **physical rule
operator basis** just described, with its constant/telescoping pieces
removed and common hop scale fixed. The projector and readout rules
retain their original normalization; this is a stated basis, not a
basis-independent norm ratio. For the original n=11 certificate the
largest coefficient is K=2^2048 and the smallest is

    h=lambda^2 alpha^4 beta^2 g^2 2^-101,
    lambda=2^-480, g=2^-216,
    alpha=sqrt(2/14)sin(pi/14), beta=1/d_13,
    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1).

Thus log2(K/h)=3814.546756... and log10(K/h)=1148.292993....
These figures describe a microscopic rule-coefficient hierarchy;
they are not a count of literal terms or an emergent velocity ratio.
The original exact root enclosure proves finite tolerances around its
symbol. It does not make this hierarchy moderate.

The new all-n partial construction supplies no improved finite schedule.
It requires small q, large D, small g and lambda, and K chosen last so
the **absolute** finite-K splitting lies below the payload tolerance.
Its leading bandwidth scale

    W=lambda^2 g^2 |v_0 v_p|^2 D^(-2p-2)

can be extremely small even before the smallest Gram eigenvalue enters
the velocity. This emergent suppression is additional to the bare
coefficient hierarchy, not a replacement for reporting it.

The selector uses mu_total=1+q^2-(p+2)E_p. Since E_p is of order
q^(2p+2), this is a delicate cancellation of order-one energies. It
has an open tolerance, not just one tuned decimal: changing mu_total
by less than gamma/4 leaves the comparison intact if the finite-K
error is below gamma/2, where gamma=E_(p+1)/(p+2). But that guaranteed
interval shrinks rapidly with p. The finite simple-band margins also
give open regions at each n, with no n-independent width claimed.

The practical outcome is smaller physical d and a sharper proof gap.
It is not a simpler or more robust alternative at comparable evidence,
and it establishes neither an unbounded-c family nor a new CFT member.
