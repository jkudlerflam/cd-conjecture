# Uniform low-frequency bounds for bare off-code sources

Status: proved operator estimates for the full microscopic Hamiltonian.
They control the size of retained source components, but do not construct
an effective action or prove a conformal infrared limit.

## 1. A local commutator inverse

Let H=A+V, where A=sum_j A_j has the atomic properties in
local_dressing.md and V is a bounded nearest-neighbor interaction with local norm

    v=sup_j sum_(Y contains j)||V_Y||.

Here v includes the weak coupling. Let O_X be Hermitian and supported
on a finite set X, with

    O_X=Q_X O_X P_X+P_X O_X Q_X.

Thus it is off diagonal relative to the local all-code subspace, not
necessarily relative to any interacting spectral subspace. With the
same reduced inverse as before set

    S_X=R_X O_X P_X-P_X O_X R_X.

[proved] S_X is bounded and anti-Hermitian and

    O_X=[A,S_X]=[H,S_X]-[V,S_X],
    ||S_X||<=||O_X||/Delta,
    ||[V,S_X]||<=2|X|v||O_X||/Delta.                          (1)

The first identity and the norm follow from the local block matrices.
Only interaction terms meeting X contribute to the last commutator.
These are exact identities, without a perturbation expansion.

## 2. True spectral windows, uniformly in ring length

[proved] Let Pi_I be the spectral projection of the full H onto any
energy interval I of width Omega. Then

    ||Pi_I O_X Pi_I||
       <=(Omega+2|X|v)||O_X||/Delta.                           (2)

Proof. Compress (1). Since Pi_I commutes with H, the first compressed
commutator is [H_I,Pi_I S_X Pi_I]. Subtract the midpoint of I from
H_I, whose norm is then at most Omega/2. The commutator norm is at
most Omega||S_X||. The other term is bounded by (1).

For a ground vector psi_0 of energy E0, the same identity gives the
one-sided spectral-weight estimate

    ||Pi_[E0,E0+Omega] O_X psi_0||^2
        <=[(Omega+2|X|v)||O_X||/Delta]^2.                     (3)

Indeed [H,S_X]psi_0=(H-E0)S_X psi_0. This proof applies also to a
degenerate ground space, one vector at a time. The possible one-point
component is included rather than silently subtracted. At Omega=0,
(2) bounds the source's complete ground-space compression by
2|X|v||O_X||/Delta.

Unlike a bare-Q resolvent, Pi_I is a true spectral projector. Its
nonlocality causes no volume factor in this estimate because it is
used only after the local identity (1). The estimate is valid at any
energy density, although the low-energy application is the useful one.

## 3. An exact smooth high/low operator split

Choose a fixed real even smooth function f_0 with f_0(x)=1 for |x|<=1
and f_0(x)=0 for |x|>=2. One explicit choice uses eta(x)=exp(-1/x) for
x>0 and zero otherwise:

    f_0(x)=eta(2-|x|)/[eta(2-|x|)+eta(|x|-1)].

The ratio is well defined for every x; at zero the function is constant
on a neighborhood, so the absolute value introduces no nonsmoothness.
Let g_0 be its inverse Fourier transform, with convention

    f_0(omega)=int g_0(t) exp(i omega t) dt,
    C_0=int |g_0(t)|dt, C_1=int |g_0'(t)|dt.

Both constants are finite because f_0 is smooth with compact support.
For Omega>0 define g_Omega(t)=Omega g_0(Omega t) and

    O_low=int g_Omega(t) exp(itH)O_X exp(-itH) dt,
    O_high=O_X-O_low.

[proved] These bounded operators satisfy

    ||O_low|| <=(C_1 Omega+2C_0|X|v)||O_X||/Delta.             (4)

To prove this, insert (1). Integrating the [H,S_X] term by parts gives
i int g_Omega'(t)tau_t(S_X)dt, of norm at most
Omega C_1||S_X||. The other term costs C_0||[V,S_X]||.
The boundary terms vanish. This is an operator-norm bound on the full
Hilbert space, not only on ground-state matrix elements.

[proved] In a ground state, O_high has no spectral weight below Omega,
including in the ground eigenspace. Its Euclidean two-source function
therefore obeys

    0<=<O_high exp[-tau(H-E0)] O_high>
                          <=||O_high||^2 exp(-Omega tau).     (5)

This follows from the spectral multiplier 1-f_0((E-E0)/Omega), which
vanishes for 0<=E-E0<=Omega. Equation (5) concerns the high-filtered
source. It does not assert exponential decay for O_X itself.

The split is quasi-local in cell coordinates. On-site A evolution does
not propagate support, while the bounded interaction V has a finite
Lieb-Robinson velocity proportional to v. Approximate tau_t(O_X) by
its evolution in a ball of radius R. Split the filter integral into
|t|<=R/(2v_LR) and its complement. The first contribution is bounded
by the exponential Lieb-Robinson tail; the second by the tail of
g_Omega. Since g_0 decays faster than every inverse power, for each
fixed power p one obtains a spatial error bounded by a constant times

    ||O_X|| [exp(-aR)+C_p(1+Omega R/(2v_LR))^(-p)].

Constants may depend on the fixed support X and the chosen filter,
but not total ring length or the upper on-cell energy. No assertion
of exponential spatial decay for a compact-frequency filter is made.

## 4. What the split does and does not accomplish

For v/Delta small, choose for example Omega=sqrt(v Delta). Then
O_low has norm O(sqrt(v/Delta))||O_X||, and O_high has a strictly
positive spectral threshold at the intermediate scale Omega, much
larger than v. All scales remain fixed before the thermodynamic limit.
At Omega of order v, the low part is O(v/Delta), but the high-source
threshold is then of the same order as the retained dynamics.

These estimates explain quantitatively the algebraic tails in
massive_source_tails.md. Their amplitudes can be small while their
asymptotic decay is slow. They also show why subtracting a one-point
function does not supply the required operator separation.

The definition uses the full interacting dynamics. It is not an
explicit change of degrees of freedom or a convergent expansion of an
irreducible effective action. The high-filtered source need not create
one independently countable massive particle, and overlapping filtered
sources are not independent variables. Thus (4)-(5) cannot be substituted
for the connected-kernel construction still missing in C8.

Exact symmetries are preserved by this filtering: conjugating by a
symmetry commuting with H can be taken inside the integral. This is
useful for a future construction, but does not by itself establish its
Gaussian stiffness or dressed momentum bounds.

## 5. Finite checks

[numerical] The last part of massive_source_tails.py tests (2)-(3)
on the literal four-cell spin Hamiltonian, dimension 256. Besides the
J=.2, theta=.1 control, it uses J=.02, theta=.02, Delta=1. The exact
perturbation local norm is bounded by

    v=Delta |sin(theta/2)|+2J.

There is one on-site difference and two incident XX bonds. Spectral
windows of widths J,2J,Delta/2 are formed from the full interacting
eigenvectors. Every compression and ground-source weight obeys the
bound. For the smaller-coupling fixture, the width-.5 window has
source norm .0199986667 and weight .0003999467, below the respective
bounds .5999996667 and .3599996000. These constants are conservative;
the point is independence of volume in the proof, not tightness at
this one size. No spectral filter constants are estimated numerically.
