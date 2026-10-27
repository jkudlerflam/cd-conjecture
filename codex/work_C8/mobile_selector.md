# A mobile selector at a cell-length crossing

2026-10-27, resumed C8. The microscopic move, its conservation laws, and its first-order projection are proved. The projected-model spectra below are numerical. No conformal finite-u microscopic simulator, logical payload transport, or c-d counterexample is claimed.

## 1. An explicit local move

Use C7's four letters #,L,X,R and its allowed arrow pairs, with

    #=(R,L),   L=(L,L),   X=(L,R),   R=(R,R).

Add a marker hopping term to the same nearest-neighbor bond,

    h=K Q + (mu/2)(M tensor I+I tensor M)
      -t (|LX><XR|+|XR><LX|)
      -u (|#L><R#|+|R#><#L|),                                 (1)

where M=|#><#| and Q is the forbidden-pair projector. All coefficients are independent of ring length. The new move has identical external arrow endpoints on its two sides. It therefore preserves the allowed/forbidden status of its own bond and both adjacent bonds, just as C7's head move does. Every Q_i still commutes with H. Individual marker positions no longer commute with H, but total marker number sum_i M_i does.

The move #L -> R# shifts a marker one site to the right. The head in the cell to its right must not occupy the immediately adjacent site, since that site must be L. Its Hermitian conjugate moves the marker left.

## 2. Projection at a crossing

Fix b>=5 and set the unperturbed chemical potential to the C7 crossing between periods b and b+1,

    mu_c=2t-delta_b,
    delta_b=(b+1)a_b-b a_(b+1),
    a_l=2t[1-cos(pi/l)].                                      (2)

Then E_b/b=E_(b+1)/(b+1)=e_c<0. C7's convexity argument shows these are the only minimizing lengths. Consequently every tiling by b- and (b+1)-cells, with each head in its clock ground state, has energy e_c L. For fixed b the gap above this entire manifold is positive on ring lengths that admit such a tiling. It is the minimum of the internal clock gaps, the excess costs of other cell lengths, and the forbidden-sector gap. The manifold can be exponentially large; that does not invalidate the gap statement above the whole manifold.

For all sufficiently large L a tiling exists: choose an integer M between L/(b+1) and L/b. The number of long cells is then

    r=L-bM,    0<=r<=M.                                       (3)

In particular, every L>=b(b+1) admits at least one such M. Primality of L is no longer an obstruction to the tiling manifold.

The normalized head wavefunction in a cell of length l is

    phi_l(i)=sqrt(2/l) sin(pi i/l),    i=1,...,l-1.

A marker hop stays within the minimizing manifold only when it exchanges neighboring lengths (b,b+1) and (b+1,b). The left head overlap is

    A_b=sum_(i=1)^(b-1) phi_b(i) phi_(b+1)(i).                  (4)

On the right, shifting the marker removes its first site, so the overlap is sum_(j=2)^b phi_(b+1)(j) phi_b(j-1). Reflection of the sine wavefunctions makes this also A_b. Thus the projected off-diagonal matrix element is

    -u A_b^2.                                                 (5)

There is no diagonal first-order term. For b=5,8,16, A_b is respectively 0.9412156469, 0.9754576110, and 0.9934958020.

Within a fixed marker-number sector, label a cell by eta_j=0 or 1 according to its length b+eta_j. In the interior, (5) is the XX hopping of eta particles,

    H_gap=-J_b sum_j (tau_j^+ tau_(j+1)^-+h.c.),
    J_b=u A_b^2,    sum_j eta_j=r.                             (6)

On a ring one must also keep the reference marker position. A hop across the reference boundary shifts that position, so a fixed physical-translation sector gives a twist of the gap-chain boundary condition and a consistency constraint from relabeling the reference marker. Dropping this bookkeeping would lose finite-size levels. The script below instead builds the physical marker-pattern graph directly, with no reference-marker redundancy.

At physical momentum zero, the periodic XX ground wavefunction on gap sequences lifts to a positive wavefunction on physical marker patterns by cyclically identifying the reference. The usual positive fixed-r XX ground state is cyclically invariant, so the lift is well defined. It has energy

    E_XX(M,r)=-2J_b sin(pi r/M)/sin(pi/M),    0<r<M.             (7)

Every physical-pattern graph component admits this positive eigenvector. On a component where the hopping graph is connected, Perron-Frobenius makes it the unique component ground state. Connectivity is checked in the examples below; a general connectivity classification is not assumed in this note.

Equation (6) is a first-order coefficient, not a volume-uniform theorem for (1) at fixed u. Excited clock states and excluded lengths give higher-order terms. Their local control and the infrared phase need a separate argument.

## 3. Why a finite-density length fluid is different from one forced defect

In the thermodynamic approximation r/M=nu, the kinetic energy density obtained from (7) is

    e(nu)=e_c-(2J_b/pi) sin(pi nu)/(b+nu),    0<=nu<=1.         (8)

It is strictly lower in the interior than at either pure cell length. Its minimizing nu_* satisfies

    tan(pi nu_*)=pi(b+nu_*),    0<nu_*<1/2.                    (9)

To see uniqueness, differentiate sin(pi nu)/(b+nu). Its numerator is pi(b+nu)cos(pi nu)-sin(pi nu), whose derivative is -pi^2(b+nu)sin(pi nu)<0 on (0,1). It starts positive and ends negative, with its zero before 1/2. The minimum of (8) is therefore nondegenerate. Large b gives nu_*=1/2-1/[pi^2(b+1/2)]+O(b^-3).

This leading description has an extensive number of each cell length. It does not consist of an otherwise rigid crystal with one massive forced defect. That distinction is how it can avoid the finite-band defect obstruction in ring_obstruction.md. Changing L can be absorbed by small changes in M and r around their thermodynamic densities.

A plausible continuation is a c=1 Luttinger liquid of cell lengths. The XX description has exact marker-number conservation, and weak local corrections at an interior density have the usual spinless-fermion stability argument. Here that is only a continuum/RG proposal. Neither a thermodynamic local elimination theorem nor finite-u microscopic conformal estimators have been supplied. The first-order graph spectrum is evidence for the mechanism, not a proof of its completion.

## 4. Numerical graph checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 mobile_selector.py

The code generates every binary marker pattern whose consecutive spacings are b or b+1, applies the legal marker hops, and diagonalizes the resulting sparse graph Hamiltonian. Energies are in units J_b. It also exhaustively checks the legality indicators of both local moves against every choice of their neighboring letters. All local checks pass.

For b=5 and nu=1/2 the results are

| L | M | r | Graph dimension | Components | E0 | L gap |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 33 | 6 | 3 | 110 | 1 | -4.0000000000 | 2.385412 |
| 44 | 8 | 4 | 385 | 1 | -5.2262518595 | 2.340611 |
| 55 | 10 | 5 | 1386 | 1 | -6.4721359550 | 2.320290 |
| 66 | 12 | 6 | 5082 | 1 | -7.7274066103 | 2.309359 |

The ground energies agree with (7) within 3e-14. The displayed first gaps are consistent with a nonzero 1/L scale. They are not two independent central-charge estimators, and no central charge is inferred numerically from them alone. The density 1/2 is a convenient tested sector, not the optimal density (9) at mu=mu_c.

For comparison, keep only one long cell, r=1, with M=6,8,10,12. The physical graph is a cycle of length L=31,41,51,61. The ground energy is -2 and the gap is 4 sin^2(pi/L). Accordingly L^2 gap increases from 39.34345 to 39.44353, approaching 4pi^2. This is the quadratic isolated-defect behavior proved in ring_obstruction.md.

## 5. What this changes for the large-c route

The explicit move (1) offers a way to address the all-ring problem without keeping defects immobile or assigning an ordinary massive hopping band to only O(1) defects. It replaces a unique cell length by a finite-density quantum mixture of two lengths.

A logical payload is an additional problem. One needs compatible code isometries V_b and V_(b+1), both with the same logical space and the required boundary compressions. In addition, the marker hop must act between neighboring codes as a scalar multiple of the identity on their payloads, or as a known harmless transport unitary. A generic pair of abstract Stinespring codes does not satisfy that condition. Adding an idle physical site at one end of one code is not enough unless both possible directions of cell growth and shrinkage are handled.

Even with this transport property, the coupled length and payload modes require an infrared analysis. Their central charges cannot be added merely because both isolated leading targets are critical. This note resolves neither the repeated-bond code construction nor the coupled conformal spectrum. It supplies a concrete alternative to the frozen selector and a testable extra code-compatibility condition.
