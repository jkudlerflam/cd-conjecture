# C93 result

[Proved, partial resolution] A block representation resolves the payload-parity obstruction, including at finite K. A fixed-cutoff construction also controls the ultraviolet Hamiltonian remainder in a prescribed strong spatial coefficient norm. The leading color fields can be absorbed into the logical chemical-potential matrix with an explicit preservation budget.

[Assumed, unresolved main target] A dressed mobile threshold and a controlled matched action for the d=13 chain have not been established. The contact-clock problem remains. The action in PROOF.md section 7 is explicitly conditional, with its known reference coefficients and missing inputs separated. This report does not claim that all four C91 obstructions have been overcome.

## The new block grading

[Proved] Root the markers and assign the terminal marker to the cell on its left. Let E_j count all D, X and M letters in cell j, including illegal configurations, and let Q_j be its microscopic charge. Then

    S=sum_j (-1)^j E_j,
    Nhat_j=(Q_j-(-1)^j E_j)/2

have the following properties. Every Nhat_j is integer-valued. Every literal Hamiltonian term preserves S and sum_j Nhat_j on an open chain or an even-cell ring. The port creates or destroys one charged letter in each of two adjacent cells, so the alternating contributions cancel. On the vacuum code, Nhat_j is the payload number.

Consequently (-1)^Nhat_j supplies the cell parity needed for the graded auxiliary embedding and the parity-twisted Gaussian extension of SEPARATE_JW.md. It survives finite-K grammar corrections. It is not a diagonal grading assigned to microscopic letters independently of their cell label. An odd-cell seam and the physical cyclic root quotient still require their boundary/frame bookkeeping.

Every virtual trajectory from the vacuum code remains in S=0. In this sector the total charge is Nhat=Q/2 and its parity is independent of the root frame. This resolves the total-parity descent for code-connected amplitudes, without claiming that other physical sectors are energetically excluded.

## Fixed cutoff and color fields

[Proved] PROOF.md equations (13)-(21) give an explicit schedule at delta=1/2 and material-site spatial exponent a=256. With output exponent one, field weight four and thirteen logical modes, the coefficient conversion costs at most

    C_blk=15 exp(179.8022806534...).

For any epsilon>0, the stated finite T,R and positive J/t make the ultraviolet remainder's operator-valued logical coefficient norm at most epsilon W_ref/2, where

    W_ref=J alpha_c lambda^2 E_unit,
    alpha_c=1/2-2/pi^2,
    E_unit=[sqrt(1/8)sin(pi/16)]^4 d_15^-2 (2^-272)^2,
    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1).

This is a strong spatial Hamiltonian bound. It does not control the retained filtered contact or furnish exponential decay in relative time for the desired influence action.

[Proved reference coefficient] At lambda=2^-608 the leading color term is D_col=diag(t^2 w_r/(JK)), with

    w_r=(1-alpha_c) sum_(l=0)^r (1/8)sin^2(pi(l+1)/16).

In H_onsite=-psi^dagger mu_log psi it is absorbed by mu_log,new=mu_log,old-D_col. Taking

    K>=max(2^2560,2^2243 (t/J)^2)

puts this coefficient below one quarter of C88's additional normalized band allowance. Absorption does not create thirteen physical tuning knobs. The centered contact coupling, clock off-diagonal terms and higher grammar corrections remain; explicit bounds for the centered bare-clock projection are in PROOF.md equations (28)-(29).

## The contact obstruction that remains

[Proved] The normalized physical block contact in the half-filled free sea obeys a mean-square bound 44/ell. Free Fermi-patch compression also gives the operator bound ||P_L(B_i-<B_i>)P_L||<=14 nu_L, where B_i=n_i n_(i+1) and nu_L is the retained momentum fraction. These statements control particular free-charge observables. They do not replace the contact after a specified material particle when it multiplies a dynamical clock operator.

[Proved reference derivative; heuristic continuum consequence] At mean vacancy number u=1 per material letter,

    d alpha/du=-1/4-1/pi^2.

Thus a smooth contact fluctuation can produce a dimension-two charge-payload density coupling. The deviation cannot all be called irrelevant. The frozen clock still has order-J leakage under a single contact, even at zero port. A dressed threshold relative to the interacting bath ground state is the missing dynamical result.

[Proved conditional constants; assumed application] If that result supplies C48's product auxiliary representation with gap Delta_* and heat coefficient at most 3*8^14, its constants can be

    L=115.7258872224..., Theta=2L,
    kappa=Delta_*/32, beta>=4L/Delta_*,
    A_dec=2^22 exp(L)*16386.

C48 then requires k+r<=Delta_*/(64Theta) and bounds the influence correction by A_dec r. Delta_*, k and r for the actual mobile chain have not been obtained. Assigning r only the ultraviolet remainder would omit the leading contact coupling and would be invalid.

[Numerical and exact finite controls] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. Exact enumeration checks 2070 directed moves in two unrestricted 576-state material sectors and all 8192 code bit strings. A sixteen-dimensional CAR and auxiliary-parity check has zero recorded residual. Small free-sea variance and constant checks pass. No parent certificate rerun, many-body d=13 diagonalization, dressed-gap calculation or central-charge measurement was performed.

[Proved, deliverables] PROOF.md contains the derivations, PROGRESS.md the work record, checks.py and checks.json the controls, and STATE.md the local handoff. No files outside work_C93 were modified. The next required work is the zero-port mobile gap-clock threshold and its material-label contact kernel estimate.
