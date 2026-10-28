# Odd-conversion qutrit candidate

This is a candidate, not a demonstrated counterexample.

## Microscopic bond

Local basis |0>, |1>, |2>, conserved occupation n=diag(0,1,2), and D=|2><2|. Orient each bond i,i+1 consistently. The Hermitian two-site term is

h = U/2 (D tensor I + I tensor D)
    - (|10><01| + h.c.)
    - tp (|20><02| + h.c.)
    - g (|20><11| + h.c.)
    + g (|02><11| + h.c.)
    + V11 |11><11|
    + V12 (|12><12| + |21><21|)
    + V22 |22><22|
    - t21 (|21><12| + h.c.)
    - mu/2 (n tensor I + I tensor n).

Every coefficient is real. Summing the same oriented h around a periodic chain is one-site translation invariant. The single-particle dispersion at mu=0 is e1(k)=-2 cos(k). The opposite conversion signs are physical, and are responsible for the nonzero pair-momentum minima. All reported energies omit mu unless a grand-canonical run is identified.

## Exact two-particle reduction

Write K for total momentum, c=cos(K/2), s=sin(K/2), and restrict 0<=K<=pi. The relative-coordinate basis consists of an onsite pair r=0 and separated atoms r>=1. A phase change of the r=0 state makes the off-diagonal entries real. The resulting semi-infinite Jacobi matrix is

H00 = U - 2 tp cos(K),
H11 = V11,
Hrr = 0 for r>=2,
H01 = -2 g s,
Hr,r+1 = -2 c for r>=1.

For V11=0, the surface Green function below the atomic continuum is

G(E,K) = [E + sqrt(E^2 - 16 c^2)]/(8 c^2),

where the square root is positive and E<-4|c|. Thus a bound state solves

E = U - 2 tp cos(K) + 4 g^2 s^2 G(E,K).

At E=-4 the self-energy becomes -2 g^2 s^2/(1+s). Fix 0<s0<1 and choose

tp = g^2 (2+s0)/[4(1+s0)^2],
U = -4 + 2 tp + g^2 s0^3/(1+s0)^2.

The Schur complement of H(K)+4 in the onsite-pair sector is then

F(s) = g^2 (s-s0)^2 [(2+s0)s+s0]/[(1+s)(1+s0)^2] >= 0.

For 0<K<=pi the separated-atom block Hsep(K)+4 is strictly positive. Positivity of its Schur complement proves H(K)>=-4, with equality only at s=s0. At K=0 the conversion vanishes, the atomic continuum touches -4, and the onsite level obeys H00+4=g^2 s0^3/(1+s0)^2>0. The bound-pair minima are therefore at ±Q, Q=2 arcsin(s0), and Epair(±Q)=-4. This establishes the two-body band minimum without relying on a momentum grid.

The relative-coordinate tail is x^(r-1), with

x = sqrt[(1-s0)/(1+s0)].

The onsite-pair residue is

Z = [1 + g^2 s0/(2(1+s0))]^-1.

Differentiating the Schur complement and using dE/dK=0 gives

Epair''(Q) = Z g^2 s0(3+s0)(1-s0)/[2(1+s0)^2].

At g=2 and s0=1/2,

tp=10/9, U=-14/9, Q=pi/3,
x=1/sqrt(3), Z=3/5, Epair''(Q)=7/15.

The atomic mass is 1/2 and the pair mass is 15/7. All units use the atomic hopping and lattice spacing as one. Direct relative-coordinate diagonalization agrees with the energy, momentum, and curvature.

## What would be required for c=3

A finite-density homogeneous phase containing atoms near k=0 and molecules near both ±Q could have three Luttinger modes. This conclusion requires all three species to have nonzero density, positive compressibility in all directions, and no relevant conversion or valley-polarization instability. The quadratic vacuum threshold itself has z=2 and is not a CFT counterexample.

In a dilute bosonization description, let theta_a, theta_+, theta_- denote the phases of the three fluids. A conversion a a -> p_+ carries momentum Q and is oscillatory away from density-dependent commensurabilities. The momentum-neutral process four atoms -> p_+ p_- has vertex 4 theta_a-theta_+-theta_-. At three decoupled Tonks fixed points, its scaling dimension is 16/4+1/4+1/4=9/2. For Q=pi/3, the three-pair valley conversion has vertex 3(theta_+-theta_-), also dimension 9/2 at that reference point. These power counts only show that a three-fluid fixed point is not immediately forbidden. The actual stiffness matrix can change the dimensions, and the microscopic model has not been shown to realize this fixed point.

Taking U=-14/9+delta changes the molecular threshold and can compensate finite-density preferences. Numerical scans use V11=0, V12=V22=20, t21=0, g=2, tp=10/9. Delta=0.9 gave large short-ring entropy fits and is the main finite-density test point. The rational form is U=-59/90 at delta=9/10.

## Tunable mass ratio, still only a vacuum construction

The exact threshold formula also allows the atom and pair masses to coincide. Setting s0=1/5 and g^2=180 gives

    g=6 sqrt(5), tp=275/4, U=269/2,
    Q=2 arcsin(1/5), Z=1/16, Epair''(Q)=2.

Hence mp=ma=1/2. This removes mass imbalance as an unavoidable feature of the two-body construction. It does not establish stability, miscibility, or a finite-density CFT for these stronger couplings. No many-body central-charge claim is made for this variant.
