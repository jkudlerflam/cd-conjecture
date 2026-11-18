# A nonzero-J clock response that keeps the vacancy spectrum

2026-11-18. Status: exact finite-volume identities and a curvature bound
for a defined second-cumulant kernel; an explicitly conditional acoustic
asymptotic calculation. No full microscopic relative-curvature theorem
or L-independent selector allowance is claimed.

## 1. Specify the approximation before differentiating it

For the original zero-port vacuum Hamiltonian write

    H_w(J,eta)=H_ch+J alpha_N H_mat,w
                        +eta J sum_j(A_j-alpha_N)h_j,   (1)

where H_mat,w=sum_j h_j. Eta=1 gives the literal contact-weighted model.
Eta!=1 is an analytic reference in ordered material coordinates, not a
proposed additional nearest-neighbour physical Hamiltonian. All original
microscopic coefficients remain fixed with L. Fix the spectator bits.

At eta=0 the ground state is Omega times the product of the sine clocks
v_b. The first eta derivative vanishes. Let omega>0 be an exact bath
excitation energy, and let delta_(b,r)>0 be a dimensionless frozen clock
excitation energy above v_b. With

    z_(b,r,s)=<r|h_(b,s)|v_b>,

the exact coefficient of eta^2 is the sum of a scalar term and independent
single-clock terms:

    E_w^[2](J)=-J^2 sigma_w^T G sigma_w-sum_c K_(b_c)(J),
    K_b(J)=sum_(bath states beta, r>0)
          J^2 |sum_s z_(b,r,s)<beta|A_(a+s)|Omega>|^2
                       /[omega_beta+J alpha_N delta_(b,r)]. (2)

This is the exact second eta coefficient for every fixed finite N,V and
J>0. It retains the full free-bath spectrum and the actual clock gap
proportional to J. It does not replace the denominator by a uniform
vacancy gap or expand that denominator in J/omega.

Different excited clocks are orthogonal. Therefore all unchanged clocks
cancel in the difference between a defect word and the target; the scalar
term is the full difference

    -J^2[2d^T G sigma_A+d^T Gd].                         (3)

The unchanged-background cross term must be retained. Equations (2)-(3)
are meaningful for a local balanced patch, a change of head number, or
either patch inside an arbitrarily long unchanged background.

## 2. A uniform-in-J curvature bound for this kernel

For omega>0 and d>=0 the elementary identity is

    d^2/dJ^2 [J^2/(omega+Jd)]
                    =2omega^2/(omega+Jd)^3 <=2/omega.   (4)

Every spectral weight in (2) is nonnegative. Consequently

    0<=K_b''(J)<=2 R_b(0),
    R_b(0)=sum_(beta,r>0) |source_(beta,r)|^2/omega_beta. (5)

This is exact, not a small-J Taylor estimate. It holds at all finite
volumes, with all the gapless bath modes present. It supplies a
volume-independent bound whenever the SOURCE-PROJECTED static R_b(0)
has one; an unweighted bound on the entire bath inverse is unnecessary.

For a replacement, the absolute clock curvature is at most twice the
sum of R_b(0) over removed and inserted clocks. Together with (3) and
SUSCEPTIBILITY.md (15), this proves the corresponding relative-curvature
condition for the resummed second-cumulant expression. The constants
depend on changed sources, not the number of spectator cells.

The distinction from the original condition (24) is crucial: higher
eta cumulants, at eta=1, have not been bounded. Formula (5) controls the
specific source of nonuniformity caused by sending the clock mass to
zero with J. It does not control dressing of the bath by the periodic
background, repeated clock flips, or connected processes involving
several clocks. A finite-volume eta series has no established
volume-independent convergence radius reaching eta=1.

As a literal check, the coefficient in (2) is compared with finite
eta ground energies for the actual (3,3) vacuum word at N=V=6,
J/t=0.2 and the n=9 value of nu. The 4158-state calculation gives
E^[2]/t=-0.00463503379648; a symmetric small-eta calculation agrees
within 1.5e-9. Its curvature magnitude falls from 0.241292924716/t
at J=0 to 0.214694749188/t at J/t=0.2, as (5) requires. This validates
the nonzero-J formula; it is not an extrapolated thermodynamic selector
comparison and does not change n=9's physical target.

## 3. An exact markerless dipole identity

Let H_b=nu 1 minus the open-path adjacency matrix, with ground energy
epsilon_b=nu-2cos theta, theta=pi/(b+1). On head positions 0,...,b-1 put

    X=diag(r-1/2),
    D_b=sum_(s=-1)^(b-1) s h_(b,s).

The actual endpoint contact halves in the original model give the exact
identity

    D_b=(X H_b+H_b X)/2,
    <r|D_b|v_b>=(epsilon_r+epsilon_b)<r|X|v_b>/2.         (6)

Because sum_s z_(b,r,s)=0 for r>0, a long-wavelength clock transition
has amplitude

    sum_s exp(iqs) z_(b,r,s)
                   =iq <r|D_b|v_b>+O(q^2).             (7)

No bias, readout field, or new operator has been added. The cancellation
comes from the true eigenvector of the original sine clock.

There is a closed expression for the spectral moment controlling the
first acoustic logarithm. For b>=3,

    M_b=sum_(r>0) delta_(b,r)^3 |<r|D_b|v_b>|^2
       =8 sin^2(theta)/(b+1)
                [cos(theta)(1+2epsilon_b^2)+2epsilon_b]. (8)

Proof: set K=H_b-epsilon_b, C=sqrt(2/(b+1)),
w_j=v_(j+1)-v_(j-1), with v_(-1)=v_b=0, and e=|0>-|b-1>.
Writing a=2C sin theta, the path recurrence gives

    Kw=a e,  KXv=-w,  KDv=-epsilon_b w-(a/2)e,
    <w,Kw>=2a^2 cos theta,
    <e,Kw>=2a,  <e,Ke>=4cos theta.

Thus M_b=<KDv,K KDv> is exactly (8). The closed formula and the direct
eigenbasis sum agree in the script for b=9,10,11.

## 4. The first acoustic logarithm is J^5 log J in this channel

Use the Gaussian acoustic part of the exact free fluid as a stated
low-energy input, with material-coordinate sound velocity v and

    <A_q(tau) A_(-q)(0)>_connected
                        =s |q| exp(-v|q|tau)+higher terms,
    s=[alpha'(r)]^2 K_gap/(2pi),
    K_gap=rho^-2, v=2t rho sin(pi rho).                  (9)

At rho=1/2, v=t and s=0.07857579749228010.... These thermodynamic
coefficients are fixed, not adjustable matching constants. Existence
of a uniformly controlled microscopic remainder in (9) is NOT proved
by the static calculations. The following result is conditional on
that acoustic description; it is not a theorem about every term of
the full energy.

Combining (7) and (9), the low-frequency part of one term in (2) is

    -J^2 s |D_r0|^2/pi integral_0^Lambda
                          q^3/[vq+J alpha delta_r] dq. (10)

For m=J alpha delta_r the integral is

    Lambda^3/(3v)-m Lambda^2/(2v^2)+m^2 Lambda/v^3
                       -m^3/v^4 log[(vLambda+m)/m].    (11)

The first nonanalytic term of this OFF-DIAGONAL SINGLE-CLOCK ACOUSTIC
CHANNEL is therefore

    +s alpha^3 M_b/(pi v^4) J^5 log(t/J).               (12)

It is not J^2 log J or a curvature divergence. A calculation that froze
the clock but divided by a fictitious bath gap would miss both (5) and
(12). The analytic J^2,J^3,J^4 coefficients in (11) depend on the
ultraviolet completion; the coefficient multiplying log(t/J) does not.
Other connected microscopic channels can generate additional terms;
we do not identify (12) as the complete model's leading nonanalyticity.

At n=9, with the original selector nu fixed for all three lengths, the
dimensionless coefficients of J^5/t^4 log(t/J) are:

| cell length | M_b | acoustic logarithm coefficient |
| ---: | ---: | ---: |
| 9 | 0.0550388827961 | 0.0000361947687787 |
| 10 | 0.0406699708189 | 0.0000267454591236 |
| 11 | 0.0309584903717 | 0.0000203589779410 |

For the three local replacements, subtraction of unchanged clocks gives

| Replacement | coefficient in E_B-E_A |
| --- | ---: |
| (10,10) -> (9,11) | +0.0000030628284726 |
| 9 length-10 -> 10 length-9 | +0.0001212385556751 |
| 11 length-10 -> 10 length-11 | -0.0000906102709491 |

This is another explicit reason not to transfer a favorable sign from
the modified nine-state biased clock or from one order to all orders. These signs do
not reverse the previously calculated J^2 terms.

At the original selector for large b, epsilon_b=O(b^-2), so (8) gives
M_b=8pi^2 b^-3+O(b^-4). The balanced logarithmic excess is O(b^-5);
either head-number patch gives O(b^-3). The general-density prefactor
in (12) is explicit in (9). It is not uniform as rho approaches one:
v vanishes and the prefactor scales as (1-rho)^-4. The previous
one-vacancy calculations therefore cannot be substituted for this
half-occupancy limit.

## 5. What is still required at eta=1

For a defined patch let Delta E be its exact energy difference and
Delta E_res the first-order energy plus the resummed second-cumulant
difference (2)-(3). The remaining microscopic task can be stated as

    |d^2/dJ^2 [Delta E-Delta E_res]|
                         <= C_rem(b,rho)/t             (13)

over some fixed 0<J/t<=j_0, uniformly in the number of unchanged cells,
or a sufficient integrated second-difference version. One also needs
the projected free-response bounds in SUSCEPTIBILITY.md (15).
If both held with explicit constants, integrating them would give a
genuine selector allowance against the known J alpha Delta_patch.
Neither (13) nor a full all-momentum static bound is established here.

For clarity, the quantitative conditional allowance is explicit. Suppose
tau>=sup_N t T_N, p>=sup_N t P_N and C_rem in (13) have been proved,
and the integrated remainder has no additional term linear in J.
With a_patch,b_patch from SUSCEPTIBILITY.md, set

    C_patch=a_patch p+b_patch tau+C_rem/2.

For compatible N>=10, alpha_N increases to alpha_infinity from below;
the explicit lower bound is alpha_min=alpha_10
=1/2-[200 sin^2(pi/20)]^(-1). Integration would give

    Delta E >= J alpha_min Delta_patch
                                      -(J^2/t) C_patch,
    0<J/t<min(j_0,alpha_min Delta_patch/C_patch)

as a sufficient positive-margin interval for that comparison. The
source sums are Delta_balanced=0.00400720928955... and
Delta_head=0.0200360464477... in either head-number direction. The
missing tau,p,C_rem are not supplied by the small-ring observations;
the conditional expression is not a numerical microscopic allowance.

The next decisive step is thus concrete: bound the determinant sums
for the local response and period-10 harmonics, and then bound the
connected eta>=3 remainder of THIS original clock/contact Hamiltonian.
The scalar mean reference and the conserved-charge-field repair do
not supply those missing bounds. No charged excitation gap, finite-port
thermodynamic selection, logical 2^-35 dressed-band allowance, or
compact unitary Virasoro CFT claim follows from this contribution.
