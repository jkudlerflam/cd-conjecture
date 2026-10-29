# C101 proofs and explicit assumptions

2026-10-29. Set Q=A/(pi M log_2 d). Velocities use the physical site spacing. A is the vacuum Casimir coefficient, not an arbitrary finite-ring effective central charge. No claim below bounds the unweighted sum of central charges at fixed d.

“Proved conditional” means that the displayed implication is proved, given the explicitly stated continuum identification or Bethe input. It does not mean that the lattice-to-CFT limit has been proved here. Finite-dimensional and rational-certificate statements are independent of those inputs.

## 1. The gauge optimization

Write w(X)=(lambda_max(X)-lambda_min(X))/2. The identity w(X)=inf_a ||X-aI|| makes w convex and unitarily invariant. For the ordinary site exchange S and a reflection-symmetric bond ShS=h, put D(o)=o tensor I-I tensor o. Then

    S D(o) S=-D(o),
    h=[h+D(o)+S(h+D(o))S]/2,
    w(h)<=w(h+D(o)).

Thus M=w(h), even when o is an arbitrary Hermitian matrix, not symmetry restricted. An additive scalar never affects M. All ungraded bonds and graded permutation matrices used below obey this reflection property.

For the Hubbard bond in a local Fock basis use fermionic reflection F. First average the gauge over total fermion parity P tensor P, which fixes h and replaces o by its even part. For even o, F(o tensor I)F=I tensor o. Applying the same averaging argument proves M=w(h). This deals with unrestricted gauges, including parity-odd gauges, without incorrectly equating ordinary and fermionic reflection.

A useful lower bound follows from a two-dimensional principal submatrix. A transition changing both site labels cannot be altered by D(o). In the examples with two equal diagonal entries and hopping t, the compression has eigenvalue separation 2|t|. The full width is at least the compressed width by the min-max principle. Hence M>=|t| whenever reflection or the preceding parity argument makes the displayed bond optimal.

## 2. Families already settled in C35

The continuum inputs, fusion spectra, and algebraic comparisons in [C35](../work_C35/RESULT.md), Sections 3 and 5, are inherited, not re-derived. Its normalization is as follows, with H_k=sum_{j=1}^k 1/j.

| Family | Bond normalization or spectrum | d | Exact M | A |
| :-- | :-- | :-- | :-- | :-- |
| SU(N) fundamental | P | N | 1 | 2pi(N-1)/N |
| SU(N), symmetric rank k | channel energies 2H_{k-j}, j=0,...,k | binom(N+k-1,k) | H_k | 2pi k(N^2-1)/(N(N+k)) |
| Spin S TB | sum_{j=0}^{2S} H_j Pi_j | 2S+1 | H_{2S}/2 | (pi/2) 3S/(S+1) |
| SO(N) vector, N>=3 | P-2K/(N-2), K=N times singlet projector | N | N/(N-2) | pi N/(N-2) |
| Zero-field XXZ | (XX+YY+Delta ZZ)/2 | 2 | 1, for -1<=Delta<=1 | pi sin(gamma)/gamma, Delta=cos gamma |

The symmetric-representation inequality established there is

    2k(N^2-1)/(N(N+k)H_k) <= log_2 binom(N+k-1,k).

It holds for every integer N>=2,k>=1. Equality requires (N,k)=(2,1). For k>=2 its Q is strictly below 3/4. The next largest discrete symmetric-family value is Q(3,1)=4/(3 log_2 3)=0.8412396714. Spin S TB is the N=2 subfamily after the factor-of-two normalization in the table. SO(N) has Q=1/log_2 N. XXZ has Q=sin(gamma)/gamma and approaches one at gamma=0. The Delta=-1 endpoint has quadratic low-energy dispersion and is outside the CFT premise.

Reflection proves that every width in this table is already gauge optimized. The spin-1/2 Heisenberg chain h=S.S has M=1/2 and A=pi/2. This is the equality case in conventional spin units.

## 3. XYZ and Sp(2r)

### XYZ

For h=J_x S_x tensor S_x+J_y S_y tensor S_y+J_z S_z tensor S_z the four Bell energies are

    (-Jx-Jy-Jz)/4, (-Jx+Jy+Jz)/4,
    (Jx-Jy+Jz)/4, (Jx+Jy-Jz)/4.

Pairwise differences are half the signed sums or differences of two couplings. If a>=b>=c are the sorted absolute values, the largest difference is (a+b)/2. Therefore M=(a+b)/4, for all real couplings. This identity was independently checked on random two-site matrices.

The critical zero-field XYZ loci are the XXZ loci and their images under axis permutations and staggered pi rotations. On these loci the two largest absolute couplings coincide. These local unitary transformations preserve the two-site spectrum and the bulk velocities; on rings one must use a compatible even-length vacuum sequence. Consequently Q=sin(gamma)/gamma<=1. A generic anisotropic XYZ chain is gapped and has no nonzero infrared CFT A. The ferromagnetic quadratic endpoints are excluded, rather than assigned a spurious CFT charge. The critical-locus input is described in [Fendley and Hagendorf, after Eq. (2)](https://arxiv.org/pdf/1006.0237).

### Sp(2r)

Let J be the real antisymmetric symplectic form, |J>=sum J_ab |ab>, and K_J=|J><J|. Then <J|J>=2r and P|J>=-|J>. For

    h=P+K_J/(r+1),

the symmetric, antisymmetric traceless, and singlet eigenvalues are

    1, -1, (r-1)/(r+1).

For r>=2 both extremes occur, so M=1. For r=1 the antisymmetric traceless representation is absent, giving M=1/2. This exceptional rank is the rescaled Heisenberg chain.

The published multivelocity identification assigns c=1 to r elementary sectors, with r-1 speeds 2pi/r and one speed pi/(r+1). The remaining r-1 Bethe branches are composite. Thus, conditional on completeness of that identification,

    A/pi=2(r-1)/r+1/(r+1)<2<=log_2(2r),   r>=2.

At r=2, Q=2/3. For r=1, A=pi/2 and Q=1. The elementary-branch dispersions are analytic, but the assignment of the complete Casimir content for general rank is argued from finite-size Bethe evidence in [Martins, Eqs. (14)-(16)](https://arxiv.org/pdf/cond-mat/0012159). It is not a uniform rigorous CFT theorem in that source. Even counting each of its composite branches once more with c=1 would give A/pi=2(r-1)/r+r/(r+1), still below log_2(2r): check r=2,3 directly, and use the upper bound 3 for r>=4. This last observation is a robustness check, not a proof excluding arbitrary additional sectors.

## 4. A positive stiffness bound

This is the tool used for Hubbard and t-J. Couple conserved on-site charges to uniform bond phases phi_a. For a smooth finite-ring ground-state branch, second-order perturbation theory gives the stiffness matrix

    D_ab = lim_{L->infinity} (1/L) partial_a partial_b E_L(phi)|_{phi=0}
         = lim (1/L)[<H''_ab>-2 sum_{m>0}
             Re(<0|H'_a|m><m|H'_b|0>)/(E_m-E_0)].

The derivative is taken before the thermodynamic limit. This is not the curvature of a thermodynamic energy minimized over all winding sectors. Degenerate shells can be avoided by a compatible length or boundary-twist sequence. For any real vector x the second term contracted with x_a x_b is nonnegative, hence

    D <= Kdia,       (Kdia)_ab=lim <H''_ab>/L

in the positive-matrix order. Consequently tr(W D)<=tr(W Kdia) for every positive W. This is just the optical sum-rule inequality, with the positive regular spectral weight retained on the correct side. The convention D here is twice the stiffness called D_c in [Stafford and Millis, Eq. (2.5)](https://stafford.physics.arizona.edu/pdf/mott.pdf).

For a complete balanced SU(q) Luttinger liquid, use an orthonormal basis in the q flavor-charge space. Its uniform charge direction is e0=(1,...,1)/sqrt(q). The other q-1 directions generate the Cartan currents of SU(q)_1. In this normalization the stiffness eigenvalues are

    D_charge=v_c K_c/pi,     D_spin=v_s/pi (q-1 copies).

The central charges are 1 and q-1, so A=v_c+(q-1)v_s. These statements are continuum inputs, including completeness. They do not hold for an unspecified extra neutral sector or for a general multivelocity dressed-charge matrix.

For q=2, twisting both physical flavors by phi gives D_total=2v_c K_c/pi; twisting them by opposite phases gives D_opposite=2v_s/pi. This fixes the normalization and agrees with free spinful fermions, for which both velocities are 2t sin(pi n/2).

## 5. Hubbard at every filling, zero magnetic field

Take t>=0 and

    h=-t sum_sigma(c†_{1sigma}c_{2sigma}+h.c.)
      +(U/2)(n_{1up}n_{1down}+n_{2up}n_{2down})
      -(mu/2)(n_1+n_2).

Chemical potential is part of the bond, not an allowed telescoping gauge. Its full sixteen-dimensional spectrum is

| Particle sector | Energies, with multiplicities |
| :-- | :-- |
| 0 | 0 |
| 1 | -mu/2 +/- t, each twice |
| 2, triplet | -mu, three times |
| 2, odd doublon | U/2-mu |
| 2, remaining singlets | U/4-mu +/- R, R=sqrt(4t^2+U^2/16) |
| 3 | U/2-3mu/2 +/- t, each twice |
| 4 | U-2mu |

The one- and three-particle sectors are two-by-two hopping blocks. In the remaining singlet block, the singly occupied singlet and even doublon have diagonal entries -mu and U/2-mu and off-diagonal magnitude 2t, yielding the displayed roots. This also proves the multiplicity count 1+4+6+4+1=16. Define Emax and Emin by the maximum and minimum of this finite list. Then

    M(U,mu,t)=(Emax-Emin)/2 >=R>=2t.

This max/min formula is an exact closed-form answer for arbitrary chemical potential and either sign of U, with gauge optimality from Section 1. At the particle-hole symmetric chemical potential mu=U/2, it simplifies to M=R. Adding a constant to use U(n_up-1/2)(n_down-1/2) changes no width.

For the repulsive balanced metallic phase the standard Bethe/Luttinger inputs are K_c>=1/2, K_s=1, one charge boson and one spin boson. At half filling the charge boson is gapped for U>0, leaving only the spin boson. The parameter conventions and Bethe results are in [Schulz, Section 4.3](https://arxiv.org/pdf/cond-mat/9412036) and his [original correlation-exponent result](https://doi.org/10.1103/PhysRevLett.64.2831). We use this parameter range as a known Bethe input, not a new proof of dressed-charge monotonicity.

Let T denote the full hopping Hamiltonian and tau=-<T>/L. Charge and opposite-spin twists have the same diamagnetic operator -T. Section 4 gives

    2v_c K_c/pi<=tau,       2v_s/pi<=tau.

Therefore, away from half filling,

    A=v_c+v_s <= (pi/2)tau(1/K_c+1) <=(3pi/2)tau.

For every interacting state with balanced filling n, its kinetic expectation is bounded below by the free hopping ground energy in that particle-number sector. Filling the free dispersion -2t cos k up to k_F=pi n/2 gives

    <T>/L >= 2 integral_{-k_F}^{k_F} (-2t cos k) dk/(2pi)
            =-(4t/pi)sin(pi n/2).

This variational bound does not assume that the interacting momentum distribution is a Fermi step. Combining the inequalities yields

    A<=6t sin(pi n/2),
    Q=A/(2pi M)<=3/(2pi)<1.                         (H)

It applies to every 0<n<2, n!=1, at finite repulsive U, conditional on the stated continuum inputs. At n=1 retain only v_s and obtain A<=2t, which is stronger. At U=0, A=4t sin(pi n/2), and the exact maximum is Q=1/pi at half filling. At n=0,2 the ground state is gapped or at a quadratic threshold, rather than a finite-velocity conformal liquid.

For balanced attractive Hubbard, the Bethe input is one gapless charge boson with K_c>=1 and a spin gap. The same argument gives A<=2t sin(pi n/2), hence Q<=1/(2pi). This attractive parameter range is kept explicit as an input; it was not independently derived in this run. The half-filled charge velocity is related by particle-hole transformation to the repulsive spin velocity. Nothing in this proof covers arbitrary spin polarization, where the charge and spin modes mix and the scalar K_c argument must be replaced by a matrix estimate.

Numerical checks use the repulsive zero-field Lieb-Wu reduction, t=1, u=U/4:

    R_u(x)=integral_0^infinity cos(omega x)/[pi(1+exp(2u omega))] d omega,
    rho(k)=1/(2pi)+cos k integral_{-Q}^{Q} R_u(sin k-sin k')rho(k') dk',
    epsilon(k)=-2cos k-mu+integral_{-Q}^{Q} R_u(sin k-sin k')cos k' epsilon(k') dk',
    epsilon(Q)=0,    n=integral rho(k)dk,
    v_c=epsilon'(Q)/(2pi rho(Q)),
    v_s=-[integral cos k exp(pi sin k/(2u))epsilon(k)dk]
            /[4u integral exp(pi sin k/(2u))rho(k)dk].

The spin formula follows from the exponential tail of the spin convolution kernel 1/[4u cosh(pi lambda/(2u))]. It agrees with integrating epsilon' against the same exponential and using epsilon(+/-Q)=0. These are the zero-field dressed-energy equations used in [Yu and Fowler, Eq. (4.8) and following](https://doi.org/10.1103/PhysRevB.45.11795). The computation is independent evidence, not the proof of (H).

## 6. Supersymmetric t-J with q equal flavors

The local space is a hole plus q flavors, d=q+1. For q>=2 define

    T=sum_{a=1}^q (|a0><0a|+|0a><a0|),
    P_occ=sum_{a,b=1}^q |ba><ab|,
    h=-tT+t(P_occ-n tensor n)-(mu/2)(n tensor I+I tensor n).

The full spectrum consists of

    0;  -mu/2 +/-t;  -mu;  -mu-2t.

The last two are occupied symmetric and occupied antisymmetric channels; the antisymmetric channel is present for q>=2. Multiplicities are 1,q,q,q(q+1)/2,q(q-1)/2. Thus

    M=[max(0,-mu/2+t,-mu)-min(0,-mu/2-t,-mu-2t)]/2 >=t.       (T1)

This is exact for every mu and q>=2, with reflection proving gauge optimality. For q=1 the antisymmetric channel is absent and the model is free hopping.

The continuum input in the interior unpolarized phase is U(1) charge times SU(q)_1 spin, with A=v_c+(q-1)v_s. We need only K_c>=1/q. One can obtain this lower bound from the electron momentum-distribution exponent in [Kawakami](https://doi.org/10.1103/PhysRevB.47.2928). The SU(q) Luttinger vertex gives

    theta=(K_c+1/K_c-2)/(2q).

The published range 0<=theta<=(1-1/q)^2/2 implies

    K_c+1/K_c<=q+1/q,
    (K_c-q)(K_c-1/q)<=0,
    1/q<=K_c<=q.

Only its lower endpoint is used. This argument avoids assuming that every graded phase has the same dressed-charge range. For q=4 the independent charge normalization and Drude relation are displayed in [Sirker, Damerau and Klümper, Eq. (5)](https://arxiv.org/pdf/0809.3427). Their K_c equals q times ours.

In flavor-charge coordinates put W=I+(q-1)e0 e0^T. It is positive, with eigenvalue q along charge and one along spin. Hence

    A<=q v_c K_c+(q-1)v_s=pi tr(W D)<=pi <sum_j B_j>/L.

For hopping a flavor past a hole, the transported charge vector is e_a, with e_a^T W e_a=2-1/q. For exchanging two distinct occupied flavors it is e_a-e_b, with squared W length 2. Differentiating the Peierls phases twice therefore gives the explicit bond

    B=t(2-1/q)T-2t P_off,
    P_off=sum_{a!=b>=1}|ba><ab|.                         (T2)

The hole-flavor blocks have eigenvalues +/-(2-1/q)t. Distinct occupied flavor blocks have eigenvalues +/-2t. Identical occupied flavors and two holes give zero. Thus B<=2t I. For q>=3,

    A<=2pi t<=pi M log_2(q+1).                          (T3)

For q=2 this bound is insufficient. Use the SU(2)-invariant ground-state expectation. Since P_off=2(S_x S_x+S_y S_y), averaging its spin rotations replaces it by (4/3)S.S. Consequently <B>=<Bbar>, where

    Bbar/t=(3/2)T-(8/3)S.S.

The following finite operator inequality is certified by exact rational arithmetic:

    Bbar_12+Bbar_23+Bbar_34 < (9/2)t I                 (T4)

on the full 3^4=81 dimensional space. To specify the certificate independently of floating eigenvalues, form C=6 sum Bbar/t. On each occupied neighboring pair C contributes diagonal -4 for equal spins and +4 for opposite spins, with exchange -8 for opposite spins. Hole-flavor exchanges contribute +9. All other entries vanish. The integer matrix 27I-C is block diagonal in (N_up,N_down). `tj_certificate.py` constructs all 15 blocks and factors each as L diag(p_j) L^T using exact fractions. Every p_j is positive, and every reconstructed entry is checked for equality. The full list of rational pivots is in `tj_certificate.json`. The largest block has dimension 12. This proves (T4), rather than estimating its validity from a floating spectrum.

Translate (T4) around a ring. Every Bbar bond occurs in three windows, so

    <sum Bbar_j>/L <=(3/2)t,
    A<=(3/2)pi t<pi t log_2 3<=pi M log_2 3.            (T5)

The strict numerical gap is unnecessary: log_2 3>3/2 follows from 3^2>2^3. Thus the q=2 result is an analytic consequence of a reproducible finite exact certificate and the stated continuum inputs. It is uniform in filling and chemical potential throughout that unpolarized conformal phase. At full filling the charge sector is absent and the SU(q) fundamental bound applies. At vacuum thresholds quadratic modes invalidate the CFT premise.

## 7. General graded permutations and what is not proved

For parities p_a in {0,1}, real fields f_a, and real J, set

    h=J P_g+(f tensor I+I tensor f)/2,  f=diag(f_a).

The exact eigenvalue list is

    J(-1)^{p_a}+f_a                       for |aa>,
    (f_a+f_b)/2 +/-J                     for each a<b.

Therefore M is half the range of this list. Ordinary site exchange fixes both terms, so no coboundary lowers it. This settles the gauge optimization for arbitrary grading and diagonal fields.

It does not settle A. A bare mixed-grading -P_g chain with an even color has a uniform even-color ground state of energy -L and a quadratic one-particle threshold. The +P_g counterpart has the analogous odd-color state. These are not complete unitary z=1 CFT examples. Finite-density phases need the fields that select them, and those fields must be retained in the displayed M.

Section 6 covers one even hole and q equally populated odd flavors, and relabelings or unitary equivalents. SU(1|1) is free hopping, with M=max(t,|mu|/2) after centering and A=sqrt(4t^2-mu^2) inside its band. Arbitrary numbers of both parities with flavor splitting, or several independent dressed seas with a nondiagonal stiffness matrix, remain open here. A name such as SU(N|M) does not specify the sign, fields, occupied seas, or whether a quadratic sector survives. No all-gradings theorem is claimed.

## 8. Adversarial constructions and the function of d

Products of reflection-symmetric chains obey M_product=sum_i M_i because their spectra add and reflection fixes the full bond, including arbitrary entangling coboundary gauges. A_product=sum_i A_i and log d_product=sum_i log d_i. Therefore any such product of verified factors obeys the conjecture. For m identical Heisenberg chains,

    d=2^m, M=m/2, A=m pi/2, Q=1/m.

Products make the ratio smaller. For general factors without reflection, the normalized-partial-trace argument in C35 Section 6 gives M_product>=max_i M_i, already sufficient for closure. Adding spectator local states also increases log d without adding A.

The spin-orbital screen takes the physical rung as d=4, with h=aP_s+bP_t+cP_sP_t. Its exact M is half the range of a s+b t+c s t for s,t=+/-1. This is a literal one-site invariant nearest-neighbor family containing independent fast legs, unequal velocities, coupled ladders, and the SU(4) point. The completed local screen and its limitations are in NUMERICS.md. It finds no two-size ratio exceeding one, but cannot certify the continuum limit of unfamiliar points.

Bethe nesting levels and strings are not independently countable CFT modes. The Sp example shows this explicitly. Staggered transfer-matrix constructions and deformations by higher conserved charges require a separate interaction-range audit. If they contain three-site or longer densities, they are not examples of the fixed nearest-neighbor problem. Grouping b original sites into a rung changes d to d^b and the velocity in rung units to v/b. Conversely, an effective logical chain on frozen b-site cells has physical velocity b v_logical. These factors cannot be omitted.

Within the verified families, the largest Q is one, attained by Heisenberg and its scale/unitary equivalents. XXZ can approach it continuously, so there is no isolated runner-up over continuous parameter families. SU(3) is the runner-up only in the discrete symmetric SU(N) family.

C35 Section 3.4 proves that symmetric SU(N) rank-k chains have a dimension-constrained envelope A/M=Theta(log d/log log d). On k=N it gives

    A/M ~ [pi/(2 log 2)] log d/log log d.

Thus a dimension-independent bound, or any replacement growing asymptotically slower than log d/log log d, is impossible. No example here requires the full log d asymptotically. The coefficient pi at d=2 is forced by Heisenberg. These facts neither prove a smaller universal envelope nor disprove the proposed log_2 d form. There is no evidence from this work requiring a replacement, and no universal proof.
