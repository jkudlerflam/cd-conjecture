# Equalizing the payload ports without freezing cell lengths

Status: proved local construction and exact first-order identities.
It builds on selector_wrapper.md and does not construct the missing
high-channel core. Its finite-coupling infrared behavior is not proved.

## 1. Problem and construction

The wrapped codes of lengths N=b,b+1 have port maps r_p alpha Gamma
and r_p alpha R, p=N-b+1. At finite vacancy hopping, r_1 and r_2 are
generally different. Thus their projected intercell interaction carries
r_p r_p' and couples the length variables to the logical payload.

This dependence can be removed exactly at first order using two qubit
tracks per physical site. The total local dimension becomes

    d_equal=4(24q+12)=96q+48.                             (1)

Choose 0<phi<pi, and put

    |+(theta)>=(|0>+exp(i theta)|1>)/sqrt(2),
    h_phi=|s_phi><s_phi|,
    |s_phi>=(|01>-exp(-i phi)|10>)/sqrt(2),
    a=I-|+(0)><+(0)|.                                   (2)

On a cell of N physical sites, track A has the uniform bond h_phi
and pin a on its first site. Track B has the same bond and pin a
on its last site. Their complete joint ground state is the product

    |G_N>= product_(j=0)^(N-1)
       |+(j phi)>_(A,j) |+((j-N+1)phi)>_(B,j).           (3)

The single-pinned helix theorem in length_filter.md proves uniqueness
and the gap bound 1/[N(2N-1)] for each track. This is a genuine
repeated NN parent, not a length-dependent Hamiltonian. Only its
ground state depends on N.

In the full selector chain, put both helix bonds on P_int. On the
cut projector P_cut=X# + R#, add pin a on track B of the left site
and pin a on track A of the right site. Each pin is a one-site
operator within the two-site bond. Thus every legal cell has exactly
the two endpoint pins specified above. All added terms are positive,
so the exclusion of selector violations and other lengths in
selector_wrapper.md is unchanged. They have zero ground energy at
every cell length. The new complete selected codes are W_p tensor G_N,
and the base gap is further bounded by the smaller helix gap at
N=b+1. No extra imposed boundary or marker condition is introduced.

## 2. An explicit positive endpoint compensator

At the first physical site, use the B-track states

    u_1=|+(-(b-1)phi)>,   u_2=|+(-b phi)>.

They are distinct and obey

    chi=1-|<u_1,u_2>|^2=sin^2(phi/2)>0.

For any desired positive expectation values m_1,m_2, define

    M_L=[m_1(I-|u_2><u_2|)
                 +m_2(I-|u_1><u_1|)]/chi.              (4)

Then M_L>=0, ||M_L||<=(m_1+m_2)/chi, and

    <u_p|M_L|u_p>=m_p,  p=1,2.                         (5)

Indeed the term orthogonal to u_p gives zero, and the other
projector has expectation chi. This proves the formula without
requiring orthogonal length flags. Choose m_p=1/r_p. At the last
site the A-track states are u_p^*, so use M_R=M_L^*.

If Atilde_a,Btilde_a are the original wrapper port operators, replace
them by physical single-site operators

    Ahat_a=Atilde_a tensor I_A tensor M_L,
    Bhat_a=Btilde_a tensor M_R tensor I_B.                (6)

Tensor factorization of the cell code and (5) give, at BOTH lengths,

    (W_p tensor G_N)^dag Ahat_a(first)(W_p tensor G_N)
          =alpha Gamma_a,
    (W_p tensor G_N)^dag Bhat_a(last)(W_p tensor G_N)
          =alpha R_a.                                  (7)

Consequently the microscopic cut bond
lambda P_cut sum_ab T_ab Bhat_a tensor Ahat_b compresses to

    lambda alpha^2 sum_ab T_ab R_a tensor Gamma_b,        (8)

independently of both adjacent lengths. All compensators are fixed
matrices once b,r_1,r_2 are chosen; they do not depend on total L.
Both ancilla states have even grading, so an existing exact payload
parity is preserved.

A smaller Hermitian compensator is often useful. Write rho_p=|u_p><u_p|,
mbar=(m_1+m_2)/2 and delta=m_1-m_2. Then

    M_small=mbar I+[delta/(2chi)](rho_1-rho_2)            (8a)

has the same two expectations, since
<u_1|rho_1-rho_2|u_1>=chi and the second expectation is -chi.
Its norm is mbar+|delta|/[2sin(phi/2)]. It is positive whenever
mbar>=|delta|/[2sin(phi/2)], but positivity is not needed for a
Hermitian boundary observable. Formula (4) is a positive choice
without any additional condition. Either gives (7),(8).

## 3. Local marker transport remains nonzero and payload independent

Consider a marker move taking neighboring lengths (b,b+1) to
(b+1,b). The transferred site starts as the first site of the right
cell and ends as the last site of the left cell. On that site the
two helix phases change as

    (A,B): (0,-b phi) -> (b phi,0).

Let U_b=diag(1,exp(i b phi)). Acting with U_b on each qubit makes
this change exactly. At the b old sites of the left cell only
track B changes, with every phase decreasing by phi. At the b
retained sites of the right cell only track A changes, again by
-phi. Thus the overlap from all unchanged physical sites is

    [exp(i phi/2)cos(phi/2)]^(2b)
          =exp(i b phi) cos(phi/2)^(2b).                (9)

Replace the forward marker conversion C_wrap in selector_wrapper.md
by the local operator, on the transferred site only,

    C_new=exp(-i b phi) C_wrap tensor U_b tensor U_b.     (10)

Its adjoint supplies the reverse move. Equations (9),(10) show that
the projected hopping remains a positive real scalar times the
identity on BOTH logical payloads, with coefficient

    J=u A_b^2 beta_wrap beta_helix,
    beta_helix=cos(phi/2)^(2b)>0.                       (11)

No global change of the helix states is applied by the Hamiltonian.
Their nonzero overlap is what allows a local move. Orthogonal length
labels at untouched endpoints would instead kill this transport.

There is no necessary exponential suppression from the new tracks.
For phi=1/b,

    beta_helix=exp[-1/(4b)+O(b^-3)],
    ||M_L||<= (1/r_1+1/r_2)/sin^2(1/(2b))=O(b^2)

when r_p are bounded away from zero. The compensation costs polynomial
operator norm rather than vanishing transport. At every fixed b all
norms are finite, and the microscopic perturbation coefficients can
be chosen below the corresponding base gap. This is resource
accounting, not a thermodynamic perturbation theorem.

## 4. Internal readout fields

The same two-state interpolation equalizes any physical operator
supported near one endpoint whose compression is a positive scalar
s_p times a fixed logical operator. For example, a two-site density
readout at the first two occupied payload sites, lifted to the first
two physical sites, acquires the probability s_p that both sites are
occupied with the appropriate tags. Tensor (4), with m_p=1/s_p,
onto track B of the first site. This remains a two-site physical
operator, and its compression is the same logical density at both
lengths. The selector projector must enable it on the first internal
bond, beginning at #, and on both possible following letters L or X.

This statement requires the original initialized core to have that
endpoint readout. It does not create the independent logical density
fields of classical_readout.md without its missing controller.

## 5. What is now exact and what is not

Within the complete base ground manifold, the first-order Hamiltonian
has a payload-independent length hopping (11) and a length-independent
logical bond (8). At fixed M its local terms are those of the gap XX
chain plus the chosen logical chain. On a ring their cyclic reference
identification still acts jointly on gaps and payloads. This is a
global momentum constraint, not a new local interaction, and it must
be retained when assigning individual levels.

If the logical target is a free n-flavor fermion chain, its oscillator
towers and the one length oscillator give the expected total c=n+1
in this FIRST-ORDER model. The marker constraint shifts current and
charge sectors; it does not justify discarding the length boson.
The joint spectrum with an XX payload is now checked in
joint_fluid.md, retaining the combined ring constraint. Its infrared
stability under the microscopic higher-order terms is not proved.

The construction removes the previous r_p r_p' obstruction exactly,
at a fixed factor-four alphabet cost. It does not give an exact
invariant product code: both the compensated bond and the marker
move have leakage. The earlier exact-encoding no-go therefore
remains intact. A high-channel initialized uniform core and control
at fixed nonzero u,lambda are still missing.

## 6. Checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 helix_equalizer_checks.py

At k=3, b=4, the two tested helix pitches are pi/2 and 1/4. The
wrapper has r_1=0.505953234954 and r_2=0.500012721419. Both complete
one-track parent spectra at N=4,5 obey the stated gap bound; the
smallest gap is 0.048943483705 and ground residuals are below 3e-16.

The actual local conversion is applied to the two-cell helix state
in dimension 4^9=262144. It gives beta_helix=0.0625 and
0.939259535981 respectively, with errors below 5.1e-14. The positive
compensator (4) equalizes both endpoints within 2.6e-14 and the
two-cell logical bond within 3.4e-14. These are finite code and
transport tests, not CFT estimators. Results are saved in
helix_equalizer_checks.json.
