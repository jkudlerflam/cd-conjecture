# Exact face reductions, charges and boundary conditions

All Pauli matrices have eigenvalues +/-1. Onsite fields are counted
once; open-chain comparisons keep full fields on the endpoints. Bulk
statements use a fixed bond before taking L to infinity.

## 1. The uniform-charge face a=0

The charge is `Q=sum Z`. For b>0 set

```
G_j=exp[-i phi_j Z/2],   phi_j=pi/4-j pi/2.
```

Conjugating every open bond gives `b(XX+YY)+cZZ`; hZ is unchanged.
For b<0 use the opposite spiral, or a common spin flip before it.
On a periodic ring the closing XY bond has phase `exp(i L pi/2)`
up to the orientation convention. On L divisible by four it is the
untwisted XXZ chain; on L=2 modulo four its XY coupling changes sign.
The ZZ term does not acquire a twist. This is an exact unitary solution
map, not an assertion that its spiral is a physical uniform symmetry.

The established Bethe-ansatz infrared classification of the spin-1/2
XXZ chain in a longitudinal field gives a single compact boson at its
ordinary finite-velocity critical points. Gapped phases, the zero-field
ferromagnetic endpoint and saturation transitions are not extra CFTs.
The twist changes finite-size sectors, not this bulk charge count.
For b=0 the chain is diagonal. We distinguish this established physics
input from the elementary unitary identity proved here. Primary sources
and their mathematical scope are listed in SOURCES.

## 2. The staggered-charge face b=0

Rotate odd sites by X. On open chains and even periodic rings,

```
H(a,0,c,h) -> sum[a(XX+YY)-cZZ]+h sum(-1)^j Z_j,
sum(-1)^j Z_j -> sum Z_j.                              (1)
```

On odd rings the closing bond is a seam. Uniform-field XXZ Bethe ansatz
does not solve (1) when h is nonzero. At h=0 the ordinary zero-field
XXZ classification applies with anisotropy `-c/|a|`.

For h!=0, the rigorous exclusions already established in
[the staggered gap note](../d2_gap_reduction/staggered_slice.md) are

```
c<=-|a|:                         gap >= 2|h|,
c>r=sqrt(a²+h²/4), even L>=4:    gap = 4(c-r).           (2)
```

The second statement is the gap above exactly two Néel ground states.
At c=r and a!=0 there are L+1 ground states and an exact quadratic
defect branch, so this endpoint is nonconformal. For c=0, h!=0 the
free spectrum has gap 2|h|. The a=0 limits are classical and require
separate ground-degeneracy conventions, not a conformal interpretation.

For completeness, the same Suzuki–Fisher theorem also gives a useful
exclusion in the full four-coupling family: the transverse bond norm is
`max(|a+b|,|a-b|)=|a|+|b|`, so `c<=-(|a|+|b|)` and h!=0 imply a unique
ground state and gap at least `2|h|`. Its positive-field hypothesis is
essential; setting h=0 does not prove a zero-field gap.

Inside the remaining strip, small h perturbs the XXZ boson by a
staggered vertex of dimension `K=pi/[2 arccos(c/|a|)]` where the
zero-field XXZ liquid exists. Irrelevance for K>2 gives a controlled
small-field liquid, while relevance generally gaps that boson. This
does not classify every finite-h neutral sector. C153's near-boundary
argument assumes a positive dressed stiffness; its large-field
two-species locking argument is a scaling/bosonization hypothesis.
We retain these qualifications rather than treating charge conservation
alone as a theorem that the complete CFT has only one sector.

## 3. Complete free classification on c=0

Let `alpha_j=(prod_(l<j)Z_l)X_j`, `beta_j=(prod_(l<j)Z_l)Y_j`.
Then Z=-i alpha beta and, on an open bond,

```
XX-YY=-i(beta_j alpha_(j+1)+alpha_j beta_(j+1)),
XY-YX=-i(alpha_j alpha_(j+1)+beta_j beta_(j+1)).          (3)
```

The signed bulk Bogoliubov eigenvalues are

```
E_±(k)=4b sin k ± sqrt(4h²+16a² sin²k).                 (4)
```

The overall momentum orientation can reverse b. Equation (4), with
the particle-hole identification `(+,k)~(-,-k)`, is sufficient to
classify every zero. A full BdG zero count without that identification
would double the number of physical modes.

* If `b²>a²` and `0<|h|<2sqrt(b²-a²)`, the equation
  `sin²k=h²/[4(b²-a²)]` gives one physical Dirac sea: C=1.
* At h=0 and `b²!=a²`, two uniform Majorana chains have nonzero speeds
  `4|b+a|,4|b-a|`. Each contributes 1/2, so C=1. This includes pure
  pairing or pure DM when the nonzero coefficient is finite.
* At `|h|=2sqrt(b²-a²)>0`, the Fermi points coalesce with quadratic
  dispersion. This is a band-edge transition, not an additional CFT.
* At h=0 and `b²=a²`, one entire Majorana band is flat at zero
  (both are flat at a=b=0). The complete spectrum is not a compact
  finite-velocity CFT, even when a dispersing Majorana remains.
* All remaining cases are gapped. In particular h!=0 and `b²<=a²`
  cannot solve the zero equation.

The full finite open-chain many-body spectrum is checked by constructing
the antisymmetric Majorana matrix and enumerating its quasiparticle
occupations. Periodic Jordan–Wigner chains require parity-dependent
boundary conditions; no unconstrained periodic BdG spectrum is substituted
for the physical spin spectrum.

## 4. Zero field: exact asymmetric Ashkin–Teller representation

Put `A_j=(alpha_j+beta_j)/sqrt(2)`, `B_j=(alpha_j-beta_j)/sqrt(2)` and
`e^A_j=i A_j A_(j+1)`, `e^B_j=i B_j B_(j+1)`. Exactly,

```
H=-sum_j[(b+a)e^A_j+(b-a)e^B_j+c e^A_j e^B_j]
       +i h sum_j A_j B_j.                             (5)
```

For h=0, each species alone is a uniform Majorana chain. Pairing its
Majoranas into cells maps alternating e bonds to the transverse-field
and exchange terms of a critical Ising chain. The product term couples
the two energies on both kinds of bond equally: the asymmetric
self-dual quantum Ashkin–Teller chain, with Ising strengths b+a,b-a.
Open even-L chains have L/2 Ising spins of each colour and the usual
endpoint convention. On rings the two colour boundary twists and
parities are constrained by the original physical fermion parity.
They are not two independent periodic spin chains with freely chosen
boundary conditions.

For even L, the Majorana products `P_A=i^(L/2)prod_j A_j` and
`P_B=i^(L/2)prod_j B_j` commute with (5) at h=0. They are nonlocal
spin charges, not additive onsite U(1) charges. Their product determines
physical parity up to a conventional sign. A uniform nonzero h breaks
the separate parities by coupling A and B on each site.

## 5. An onsite four-site implementation of the duality

Define Hermitian unitaries

```
P=(X+Y)/sqrt(2), Q=(X-Y)/sqrt(2),
W_j=I,P,Z,Q for j=0,1,2,3 modulo four.
```

For every open bond, direct conjugation gives

```
(W_j tensor W_(j+1)) [a(XX-YY)+b(XY-YX)+cZZ] (...)†
     =b(XX-YY)+a(XY-YX)-cZZ,
W_j Z W_j†=(-1)^j Z.                                  (6)
```

Thus h=0 has the exact bulk duality `(a,b,c)->(b,a,-c)`.
This is the duality of C152, now implemented explicitly by a product
of onsite rotations. It preserves every bipartite entropy exactly.
The physical Hamiltonian remains the original uniform one.

On L divisible by four the periodic identity is exact without a seam.
For even L=2 modulo four, the transformed closing transverse bond has
the opposite sign from the other bonds; its ZZ coefficient still equals
-c. The missing twist explains the unequal untwisted six-/ten-site
spectra. Odd rings have a different seam and are not used here.
For h!=0, (6) produces a staggered field: it is not a duality within
the original four-parameter uniform-field family.

At `a=b`, (6) identifies c and -c, while the common rotation in
NORMAL_FORM gives the simple two-term Hamiltonian `2a XY+cZZ`.
Neither identity makes that interacting model quadratic or proves its
complete infrared content at order-one c/a.
