# Symmetry and momentum ledger

Tags: microscopic identities are **[proved]**; low-energy statements below
are **conditional derivations** in an explicitly stated two-valley field
content. Symmetry alone does not establish that content for V4.

## 1. Exact microscopic actions, including phases

Use the original frame
H0=sum[(XX-YY)/2+(XY-YX)/16+ZZ]-(7/4)sum Z.
The uniform U=exp[-i(pi/8)sum Z] gives the displayed Hamiltonian. Put
n_j=(1-Z_j)/2, b_j=(X_j+iY_j)/2 (annihilator), P=(-1)^N. On an even
periodic chain or in the infinite bulk define T b_j T^-1=b_{j+1},
I b_j I^-1=b_{-j}, K as Z-basis complex conjugation, and

R=exp[-i sum_j (pi/2+pi*j)n_j], A=R K, Theta=I K, Sigma=R I=A Theta.

Their actions are

| Symmetry | b_j image | Type |
|---|---|---|
| T | b_{j+1} | unitary |
| P | -b_j | unitary |
| A | i(-1)^j b_j | antiunitary |
| Theta | b_{-j} | antiunitary |
| Sigma | i(-1)^j b_{-j} | unitary |

All preserve H0. A²=Theta²=1, **Sigma²=P**, and A T A^-1=P T.
The displayed frame has their U-conjugates. Open boundaries do not retain
bulk translation, so these are bulk selection rules, not statements that
boundary pair sources or boundary mixing vanish.

For O_r=sum_j b_j^dagger b_{j+r}^dagger,
Sigma O_r Sigma^-1=(-1)^(r+1) O_r. Thus uniform *undressed* even-r pair
sources are forbidden, including a complex coefficient; odd-r sources
are allowed. Density-assisted sources have their own symmetry action:
this is not an all-orders ban on every operator containing an even-r pair.
This refines the scope of [C179 section 2](../../codex/work_C179/PROOF.md).

## 2. Bosonic valleys and field convention

Assume b_j contains exp(ik_A j) Psi_A + exp(ik_B j) Psi_B, with two distinct
incommensurate centres. A or Sigma exchanges these components, so
k_A+k_B=pi modulo 2pi. In an unbroken state their effective densities obey
rho_A=rho_B=rho/2. These are *dressed Gaussian species densities*; rho is
not asserted to equal the microscopic expectation <N>/L.

Use Psi_alpha proportional to exp(i theta_alpha) times integer density
harmonics exp[2im(pi rho_alpha x-phi_alpha)], and
[phi_alpha(x),partial_y theta_beta(y)]=i pi delta_alpha,beta delta(x-y).
The smooth species density is rho_alpha-partial_x phi_alpha/pi.
Set phi_c,s=(phi_A +/- phi_B)/sqrt(2), likewise theta.
A compatible field origin gives:

| Symmetry | phi_A(x) | theta_A(x) | A/B exchange |
|---|---|---|---|
| T | phi_A(x+1)-pi rho_A | theta_A(x+1)+k_A | no |
| P | phi_A(x) | theta_A(x)+pi | no |
| A | phi_B(x) | -theta_B(x)-pi/2 | yes |
| Theta | -phi_A(-x) | -theta_A(-x) | no |
| Sigma | -phi_B(-x) | theta_B(-x)+pi/2 | yes |

The B images follow by exchange. Shifting field origins rephases amplitudes
but cannot remove an allowed perturbation. The constants fix Sigma²=P;
ignoring them risks incorrectly assigning a phase to charge-changing vertices.
For example Sigma: theta_c -> theta_c(-x)+pi/sqrt(2), phi_s -> phi_s(-x).

## 3. Complete local-vertex and resonance rule

In this **bosonic** compactification, label local vertices by integer vectors
n=(n_A,n_B), m=(m_A,m_B):

V(n,m)=exp[i(n dot theta-2m dot phi)],
q(n,m)=n_A k_A+n_B k_B+2pi(m_A rho_A+m_B rho_B) modulo 2pi.

The Hamiltonian requires n_A+n_B even and q=0. This includes products
and density-dressed processes, not just bare pairing. Let S exchange A,B.
The other symmetry actions, with complex coefficients conjugated for the
antiunitary operations, are

A: V(n,m) -> i^(n_A+n_B) V(Sn,-Sm),
Theta: V(n,m;x) -> V(n,m;-x),
Sigma: V(n,m;x) -> i^(n_A+n_B) V(Sn,-Sm;-x).

Consequently a resonant vertex must be combined with its Sigma partner.
In this field origin Theta permits real cosine coefficients. For even total
charge N, the allowed orbit combination is proportional to
cos V_phase(n,m) + (-1)^(N/2) cos V_phase(Sn,-Sm).
If the two cosines coincide and N=2 modulo 4, it vanishes; that candidate
is forbidden. This orbit projection, parity, and the momentum equation are
the general selection algorithm, including resonances. Derivatives do not
cancel a nonzero lattice momentum. Allowed derivative descendants have
larger dimension; linear gradients shift density/current rather than lock
fields. The generic primary cosines below are already symmetry allowed,
so descendants create no additional stability condition when all of those
primaries are irrelevant.

At symmetry and balance set k_A=k, k_B=pi-k. Then
q=(n_A-n_B)k+n_B*pi+pi*rho*(m_A+m_B).
For generic independent k and rho, q=0 identically requires
**n_A=n_B=2r, m_A=-m_B=s**, r,s integers. Their orbit phase is one, so all
cos[2sqrt(2)(r theta_c-s phi_s)] are allowed. This agrees with C179/C180
on the generic integer lattice, after fixing the operator frame.
It means emergent relative charge conservation and total charge modulo four
in this generic continuum theory, not an extra exact microscopic U(1) or Z4.

## 4. Leading processes and special resonances

The table lists representatives; symmetry partners and Hermitian conjugates
are understood. M=m_A+m_B. All dimensions, including density dressing, must
be evaluated by the full covariance formula in GAUSSIAN.md.

| Process | n | m for simplest term | Momentum at balance | Simplest dimension |
|---|---|---|---|---|
| Relative backscattering | (0,0) | (1,-1) | 0 | 2 G(phi_s,phi_s) |
| Quartet | (2,2) | (0,0) | 0 | 2 G(theta_c,theta_c) |
| Mixed quartet/backscattering | (2r,2r) | (s,-s) | 0 | x(r,s) |
| Single valley transfer | (-1,1) | (0,0) | pi-2k | G(theta_s,theta_s)/2 |
| Two-particle valley transfer | (-2,2) | (0,0) | -4k modulo 2pi | 2 G(theta_s,theta_s) |
| Inter-valley pair | (1,1) | (0,0) | pi | G(theta_c,theta_c)/2 |
| Intra-valley pair | (2,0) | (0,0) | 2k | [G(theta_c,theta_c)+G(theta_s,theta_s)]/2 |
| Pure charge umklapp | (0,0) | (p,p) | 2pi*p*rho | 2p² G(phi_c,phi_c) |

Dressing any row by arbitrary m adds pi*rho*M to its number-changing
momentum. In particular, an inter-valley pair can resonate when
rho*M is an odd integer, intra-valley pairing when 2k+pi*rho*M=0 modulo
2pi, and pair transfer when -4k+pi*rho*M=0 modulo 2pi. Density-only
processes require rho*M an even integer. At rho=1/3, for example, the pure
charge p=3 umklapp is allowed, but other dressed charge-changing vertices
must also be tested; it is not an exhaustive one-cosine description.
Some candidates are killed by the orbit projection. Listing a resonance
condition is not a claim that its amplitude is nonzero or that V4 satisfies it.

For a relevant nonzero detuned vertex, oscillations cut off the flow at a
scale set by 1/|q|; a locking tongue can extend around exact resonance.
Finite window peaks do not bound that tongue's width. Symmetry-breaking
valley polarization invalidates the balanced ledger; reapply the general
q equation with the actual densities and populated fields.

## 5. Bosonic peaks are not fermionic Fermi endpoints

The saved physical boson peaks near k_A/pi=-0.05884 and
k_B/pi=-0.94116 have sum -1 and folded separation 0.88232; their window
resolution is about 0.056pi. They support the symmetry-related-centre
assignment, not a controlled four-endpoint fermionic fit.

If an emergent *fermionic* two-sea description is independently justified,
use its actual p_alpha,R/L and define K_alpha=(p_R+p_L)/2 and
kF_alpha=(p_R-p_L)/2. In that convention vertices
exp[i(n dot theta-p dot phi)] satisfy p_alpha=n_alpha modulo 2 and have
q=n dot K+sum p_alpha kF_alpha. The symmetric endpoint relations are
p_A,R+p_B,L=pi and p_A,L+p_B,R=pi. They do not imply equal left/right
speeds within one valley. Identifying kF_alpha=pi*rho_alpha needs the
appropriate conserved emergent charge and reference sea.

A physical Jordan-Wigner fermion multiplies b by the *total* density
string, introducing both species' phi fields and momentum harmonics.
It is not automatically a fermion of one valley with only that valley's
string. The fermionic parity constraint on p must not be imposed on the
bosonic operator lattice while retaining bosonic peak centres. C180's
special-resonance table used the fermion frame; it cannot be numerically
specialized using the saved boson peaks and raw density. C189/C191 already
identify the frame issue, though their exponents at other V and density
are not measurements of V4.
