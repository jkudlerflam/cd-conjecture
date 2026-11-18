# What this seed does and does not settle

## The narrow result

The strongest conclusion is that the following implication is false at
the audited perturbative level:

> Reciprocal spectra with several opposite-momentum pairs must have a
> relevant pair-transfer coupling and cannot support several liquid modes.

The literal nine-state seed has exactly the resonance in question, a
nonzero generated coupling, and a positive leading one-loop irrelevant
margin. C86/C87 did not assert the implication above; their restricted
certificate exclusions remain intact. The seed also demonstrates why
indefiniteness of Herm(T) is not a mode-count obstruction once a physical
non-scalar onsite matrix is included.

The proposed continuum central charge is total unitary compact Virasoro
c=2 in a selected periodic frozen sector. It is not measured numerically,
and it is not inferred merely by counting bands. The calculations establish
a finite low-flavour effective mechanism, not arbitrary c at d=9 and not
a violation of c<=d-1.

## Ground-sector and limit ordering remain material

1. **The code is the legal, K=infinity code.** The 81 by 81 physical bond
   is defined on all states and respects charge and frozen markers there.
   However, its vacuum walk and added background swaps can leave the legal
   words at finite K. Their induced splitting is not covered by the code
   projection alone. A finite-K budget must make this splitting small
   compared with the narrow conducting scale and its interaction margin.
2. **Frozen is not mobile.** Each physical marker projector is conserved.
   No vacancy, vacancy motion, phonon or positional mode is included.
   One cannot add a mobile-mode contribution to the reported c=2.
3. **The zero-port selector is not the nonzero-port selector.** The
   inherited vacuum clock has energy `-2cos(pi/(b+1))` in a length-b cell.
   Its midpoint selector favors b=4 in the legal vacuum family. With the
   primary Delta fixed, charged cells at g=0 are above their corresponding
   vacuum by Delta, and conversion can lower them by at most |g|.
   This helps the zero-port comparison but does not compare the interacting
   bulk energy against dressed wrong-length, illegal or mixed sectors at
   nonzero lambda. That comparison must be redone for the added terms.
4. **The expansion is finite-cell, not uniform infrared matching.**
   Equation (2) of INTERACTIONS.md controls conversion order. It does not
   prove a volume-uniform Schrieffer-Wolff reduction across a gapless code.
   Multi-bond port processes first enter at higher port order and require
   an extensive local error bound.
5. **The hierarchy is substantial.** A useful asymptotic ordering is
   `g -> 0`, `lambda << g²`, then K large enough that its local errors
   are smaller than `lambda² alpha⁴ g⁴`. This separates higher port
   interactions from the computed leading conversion quartic. It is an
   ordering of sufficient small parameters, not a certified finite tuple.
   The order-one walk, conversion g, port lambda, large K, and controls
   of order lambda² g² are distinct microscopic scales. Small physical d
   does not make this a simple-coefficient Hamiltonian.
6. **The RG inference is conditional.** The primary leading margin is
   strict, so higher conversion corrections at fixed n become relatively
   small as g tends to zero. Promoting that to a microscopic CFT needs
   the usual analytic weak-coupling flow control and the matching/sector
   hypotheses just listed. The present output is not a rigorous all-scale
   CFT certificate.

## Where an arbitrary-flavour extension gets stuck

The two-mode onsite synthesis uses a special fact: the real symmetric
2 by 2 matrices have dimension three, exactly the span of the three
literal controls `N,N_L,LL+RR swap`. At n=3 the same control space still
has dimension three, while a general real symmetric matrix has dimension
six. At arbitrary n it supplies only

```
M_n=mu I+eta diag(p_i)+kappa B_n,
p_i=sum_(j>i) v_j²,  (B_n)_(i,i+1)=1-v_(i+1)².
```

It is incorrect to choose an arbitrary convenient M_n in the full Bloch
symbol and pretend this physical bond implements it. In particular,
the two-dimensional absolute-value/congruence root identity is not an
all-n construction. This is a concrete control-space limitation of the
current seed, not a universal dimension or mode bound.

There is a literal next test without more states. Give LL and RR swaps
independent coefficients and add the boundary density projector
`P_M tensor |L1><L1|`. In the legal code this supplies

```
M_n=mu I+eta diag(p_i)+b p_0 |0><0|
       +kappa_L B_L+kappa_R B_R,
(B_L)_(i,i+1)=sum_(j>=i+2) v_j²,
(B_R)_(i,i+1)=sum_(j<=i) v_j².
```

The boundary term probes the first payload bit after a frozen marker;
it is an actual two-site operator, not an individually addressed logical
field. At n=3, `p=(11/12,2/3,1/3)`, the two hopping profiles are
`(2/3,1/3)` and `(1/12,1/3)`, and their determinant is 7/36. The three
diagonal profiles are independent. These five literal controls therefore
span every real **tridiagonal** 3 by 3 onsite matrix. They do not supply
its missing (0,2) matrix element.

The sharply located continuation is to establish, within this actual
five-parameter structured family, 2n simple roots and a stable matrix of
generated resonant channels for growing n. The carrier and code already
exist at physical d=9 for every finite n, but neither the needed band
pattern nor interacting stability has been demonstrated for every n.
The carrier gap also decreases with n, so no n-independent perturbative
radius follows from the two-flavour result. This route is open; no
arbitrary-c certificate is claimed by this contribution.
