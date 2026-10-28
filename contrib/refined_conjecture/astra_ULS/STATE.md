# SU(N) ULS crossover calculation, 2026-10-27

Closed-book calculation in this directory. No web or literature consulted. The parent project STATE.md was read. The parent is outside the writable workspace, so this local file records the result.

Exact results derived from recalled nested Bethe equations, using the usual thermodynamic string description:
- H=sum P, v=2*pi/N exactly.
- e0=1+2[psi(1/N)+EulerGamma]/N = -1+pi^2/(3N^2)-2*zeta(3)/N^3+pi^4/(45N^4)+...
- Rank-r hole density rho_r(lambda)=sin(pi*r/N)/(N[cosh(2*pi*lambda/N)-cos(pi*r/N)]).
- With alpha_r=pi*(N-r)/N, hole dispersion epsilon_r(p)=(2*pi/N)[cos(p-alpha_r)-cos(alpha_r)]/sin(alpha_r), 0<=p<=2*alpha_r. Bandwidth (2*pi/N)*cot(pi*r/(2*N)). Rank r transforms in wedge^r fundamental.
- Smallest dispersion-curvature scale is momentum 1/N, energy 1/N^2. Broad r=1 branch has O(1) bandwidth but already deviates from its endpoint tangent at energy O(N^-2).

Numerical thermodynamics:
- tba_check.py solves the partially decoupled rank/string TBA, with x=log Y and kernel 1/(2 cosh(pi lambda)). It differentiates the TBA to obtain entropy, rather than differencing free energies.
- roots_tba.py finds the first relevant 10% crossing. scan_tba.py samples smaller T to distinguish the upper from lower crossing.
- Cutoff M=40 (48 for N=2), 1024 rapidity points, |lambda|<=max(24,6N). Residuals <2e-9. These are numerical estimates, not certified bounds. Quote about three significant digits.
- N=2: T*=0.47516548, ell=6.61158, sigma=.251357 (upper crossing).
- N=3: T*=0.21323571, ell=9.82197, sigma=.213505 (upper crossing).
- N=4: T*=0.12772011, ell=12.29874, sigma=.202687 (upper crossing).
- N=5: T*=0.08808426, ell=14.26631, sigma=.200676 (upper crossing).
- N=8: T*=0.09948689, ell=7.89449, sigma=.401882 (lower crossing). The entropy overshoot is less than 10%, so the definition permits a longer interval than for N=5.
- N=2 cutoff check M=80, K=2048, R=40 at T=.47516548 gives s/s_CFT=1.099998648. At T=100, s=.693109434456424, consistent with log(2)-3/(8*T^2)+O(T^-3).
- N=8 at t=N^2*T=2.5: M=32,R=40 gives ratio1.09259108; M=64,R=48 gives1.09245282. Thus the absence of a 10% overshoot is robust to these cutoffs.
- N=8 at t=2 gives ratio1.08139286; N=16 at t=2 gives1.07395687, at t=6 gives.91374248.

Large-N status, do not upgrade to a theorem:
- A candidate scaling limit is N ideal fermion flavors, density1/N per flavor, quadratic dispersion k^2, with the total density fixed. This matches both fundamental particle/hole dispersions, leading e0, and leading CFT entropy.
- ideal_gas.py computes its canonical scaling entropy S(t), where 1=(1/pi) integral_0^infty dq/(1+exp((q^2-m)/t)). S=(1/pi) integral h_F(n_q)dq.
- S(t)=t/6+t^3/(45*pi^2)+... at small t. Its maximum ratio 6S/t is1.06946707 at t=2.45739154. Its first 10% failure is the lower crossing at a=6.107952120559406.
- Conditional predictions: T*=a/N^2, ell*=1.02868935*N, sigma*=.916192818/log N.
- Conditional intermediate asymptotic: s=log(N*sqrt(T))+1.5-log(2*sqrt(pi))+o(1), for N^-2 << T <<1 in a joint low-temperature, large-N limit. Constant=.2344878765.
- This calculation has NOT derived the uniform large-N thermodynamic reduction from the full TBA. Finite-number spinon scattering becoming weak does not suffice because color/string sums grow. N^-2 and the numerical prefactor are well-supported estimates, not a controlled large-N theorem here.

Finite volume:
- For the identity-sector CFT comparison use L divisible by N; other congruence classes have a different charge sector.
- Fixed conformal levels require L/N large to suppress dispersion curvature. Uniformly controlling every state through Delta~N/4 requires, heuristically, L/N^2 large: one descendant may put m~N units into a single momentum, whose relative curvature correction is O(m*N/L). Typical states and primary families can have less stringent requirements. Marginal logarithms are a separate correction.

Additional cutoff checks:
- N=8, M=64, K=1024, R=48 at T=.0995 gives ratio .89993090576. The M=40 root is stable at the quoted three significant digits.
- N=32, M=24, K=1024, R=128 at t=6 gives ratio .91157046581, entropy .88308388875, equation residual <1e-9. String-cutoff convergence must be checked before treating this as an extrapolation point.
- A Schur-Weyl segment argument alone cannot prove the thermodynamic scaling. Cutting the chain into length-L blocks loses O(1/L) energy per site; for L~N this is much larger than the proposed N^-2 thermal scale.
- The N=32,M=40 follow-up was interrupted during the entropy tangent solve to stay within the 20-minute budget. Do not use N=32 as a converged large-N benchmark. All reported finite-N onset calculations and their stated checks completed.
- Ideal-gas free-energy checks: F(2)=2.9450300204030384, S(2)=.3542649211507684; F(6)=.29449561044496875, S(6)=.9059908101588092.
