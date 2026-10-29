# C102 result

2026-10-29

[Proved conditional on the standard XXZ Bethe/CFT formula] A is not convex or subadditive. The failure occurs inside one critical phase, with c=1 throughout and constant gauge-optimized M. On the critical XXZ cone the opposite statement holds: A is concave and superadditive, strictly away from proportional couplings. Gauge optimization cannot change either conclusion.

## The counterexample

Use S=Pauli/2 and h_delta=Sx tensor Sx+Sy tensor Sy+delta Sz tensor Sz. For -1<delta<1,

    A(delta)=pi sqrt(1-delta^2)/(2 arccos delta),
    M(delta)=1/2.

[Proved conditional] Since h_0=(h_-1/2+h_1/2)/2,

    A(h_0)=1 > 9 sqrt(3)/16
                 =[A(h_-1/2)+A(h_1/2)]/2.

Equivalently A(h_-1/2+h_1/2)=2>9 sqrt(3)/8. The midpoint defect is 0.0257214207425. The velocity input is [Sirker and Bortz, Eq. (13)](https://arxiv.org/pdf/cond-mat/0511272). The curvature proof is in [PROOF.md](PROOF.md), Sections 2 and 3.

[Proved conditional] For h(J,K)=J(XX+YY)+K ZZ in spin units, J>0 and |K|<J, write A=(pi J/2)f(K/J). The full derivative calculation proves f'>0 and f''<0. Therefore every convex subcone containing two anisotropies violates subadditivity. Restricting to nonnegative antiferromagnetic couplings does not help. The segment delta=0 to delta=1/2 already violates convexity.

## The surviving comparisons

[Proved conditional] A increases with either J or K separately while remaining in the XXZ cone. At fixed M it increases towards Heisenberg. Concavity gives the supporting-plane bound

    A(J,K)<=pi(2J+K)/6<=pi J/2=pi M.

This proves the proposed bound for XXZ by monotonicity towards the permutation endpoint. It does not bound the gradient norm, which still diverges near the ferromagnetic endpoint as described in the round-one update.

[Proved] A is convex on any cone where the bulk ground-energy density e is linear and the relevant Casimir limits exist on a common sequence. A common infinite-volume minimizing state is sufficient for this linearity. The identity behind this statement is

    Jensen_defect(A_L)=(6/pi)L[L Jensen_defect(e)-Jensen_defect(E_L)].

Without bulk linearity, the two nonnegative energy defects have no ordered difference. No nontrivial interacting critical cone satisfying the sufficient condition was identified.

[Proved] On the restricted complex-hopping subspace h(t)=t|01><10|+conjugate(t)|10><01|, A=2|t| and M=|t|. Thus A is convex there, and checking the extreme circle of the restricted unit disk proves A<=2M<pi M. Periodic Casimir comparisons require compatible hopping phases and lengths. This is a small explicit subcase of the already established quasi-free bound.

[Proved] Global concavity also fails. The critical hopping bonds t=exp(+/-i pi/3)/2 both have A=1, while their midpoint t=1/4 has A=1/2. All three admit the same vacuum sequence L divisible by six. No global operator-order monotonicity can survive scalar gauge invariance either.

## Why extreme bonds do not solve the general problem

[Proved] Supremizing A over scalar and coboundary gauges gives A itself. The quotient M-unit ball is the image of the Hermitian contraction ball, so every extreme quotient point has a two-level representative. That geometric statement is valid. Convexity of the physical A is the missing ingredient.

[Proved] There is a direct failure using exposed extreme points. The Pauli bonds X tensor X and Y tensor Y have M=1 and A=0 on even rings. Their midpoint has M=1 and A=2. Thus even the actual A values at two exposed extreme bonds do not control their interior by convex interpolation. Any convex majorant must overestimate A on some extreme bonds. General projector and two-level bonds remain unclassified; no universal bound is proved here.

## Numerical checks and their limits

[Numerical] [NUMERICS.md](NUMERICS.md) records XXZ through L=16, full-space XYZ through L=12, spin-1 bilinear-biquadratic chains through L=10 on even sequences or L=9 on period-three sequences, four seeded rotated and gauged XXZ triples, and two seeded critical Ising bonds. There are 136 recorded eigenpairs. The largest matrix dimension is 12870 and the largest residual is 2.29e-11. All validation assertions passed.

[Numerical] The XXZ midpoint defect at L=16 is 0.025578090905, with the correct sign and close to the Bethe limit. All four seeded XXZ triples also have positive midpoint defects, analytically and by ED. The spin-1 critical segment b=1,1.5,2 in h_b=S.S+b(S.S)^2 gives positive defects 0.084608969 from L=6,9 and 0.073220042 from a three-size fit. These support the same direction but do not establish its thermodynamic sign. TB and ULS provide Bethe calibrations; AKLT provides an exact zero-energy-subtraction check. Massive XYZ and Haldane samples show substantial misleading short-ring coefficients.

[Scope] Bethe checks used closed velocity formulas and bulk-energy quadrature, not finite-size Bethe-root solutions. Random bonds were constrained to known critical families. No generic random criticality classification, large spin-1 extrapolation, or arbitrary-projector proof is claimed. The exact XXZ counterexample already closes the proposed general convexity route. Only this folder was modified; the project-level STATE.md was left untouched under the explicit confinement rule.

[Files] Reproduce with `python3 checks.py` and `python3 summarize.py`. The outputs are numerics.json and NUMERICS.md. PROOF.md contains the full arguments, and PROGRESS.md and STATE.md record completion.
