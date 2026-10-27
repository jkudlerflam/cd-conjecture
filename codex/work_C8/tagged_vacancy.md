# Mobile boundary initialization for a frustration-free core

Status: proved in the specified invariant ordered-tag sector. The
selection of that tag sector in an unrestricted ring remains open.
There is no high-channel core or new CFT in this note.

## 1. Repeated local bond

Let a q-state core of k>=2 sites have positive operators h0,a,b and

    H_core = g sum_{j=1}^{k-1} h0_{j,j+1} + a_1 + b_k,
    ker H_core = range(V),   V^dag V=I_D,   g>0.           (1)

This is frustration freeness with one-site boundary penalties, a
larger class than the boundary-free cores in vacancy_parent.md.

Enlarge the alphabet to

    C|0> direct-sum (C^3_tag tensor C^q_payload),
    tags L,B,R,   d_ext=3q+1.                             (2)

Here 0 denotes a vacancy. Let n=|0><0|, and let h0_tilde act as h0
on the payload of any two occupied sites, preserving their tags, and
vanish otherwise. Put a_L=|L><L| tensor a on the occupied sector,
b_R=|R><R| tensor b, both zero on all other letters. Use the repeated
Hermitian bond

    h_ext = U(n tensor I + I tensor n - 2n tensor n)
      - w sum_{tag,alpha} (|0;tag,alpha><tag,alpha;0| + h.c.)
      + g h0_tilde + a_L tensor I + I tensor b_R,          (3)

with U,w>0. The orientation of the last two terms matters.

On an open chain of N=k+p sites, restrict to the invariant sector
whose occupied tags, read from left to right, are L B^(k-2) R.
Vacancy hopping preserves this order. The sector is unitarily
identified with

    C^(binomial(N,p)) tensor (C^q)^(tensor k).             (4)

In this sector (3) gives

    H_ext = H_ch tensor I
       + g sum_j D_j tensor h0_j + I tensor(a_1+b_k).     (5)

H_ch has diagonal U times the number of vacant-occupied interfaces
and hopping -w for each vacancy-particle exchange. D_j projects onto
charge configurations in which occupied particles j and j+1 are
adjacent. The L particle is never at site N because k>=2, so a_L in
the left slot of (3) is counted exactly once. Likewise R is never at
site 1, so b_R in the right slot is counted exactly once. There are
no missing half-weights or charge-dependent endpoint coefficients.

## 2. Complete code and boundary maps

The charge graph is connected. Let chi_p>0 be its normalized unique
ground vector and e_p its energy. Positivity in (5) shows H_ext>=e_p.
Frustration freeness of (1) shows that

    V_p = chi_p tensor V                                  (6)

attains this bound. Conversely, any ground state must lie in the
ground space of H_ch tensor I, so it is chi_p tensor phi. Every
<chi_p|D_j|chi_p> is strictly positive, as every allowed charge
configuration has positive amplitude. Zero expectation of the
remaining positive terms in (5) forces every h0_j, a_1 and b_k to
annihilate phi. Therefore (6) is the complete ground code in this
sector. No reflection symmetry of the payload is assumed.

Extend an original left endpoint operator F to |L><L| tensor F,
zero on vacancies and other tags. The physical first site contains
the L particle whenever it is occupied. Hence its compressed action
is r_p V^dag F_1 V, where

    r_p=1-Prob_chi_p(first site vacant).

The right endpoint operator |R><R| tensor F similarly gives
r_p V^dag F_k V. Equality of the probabilities follows from the
reflection symmetry of H_ch alone. Thus boundary initialization
does not change the earlier scalar dilution factors.

For the two lengths p=1,2, a shift of the cut between cells can
preserve the occupied counts only by transferring a vacancy. All
ordered tags and payload labels stay in their respective cells.
The exact overlap therefore remains

    (V_2 tensor V_1)^dag (V_1 tensor V_2)=beta I_(D^2),
    beta=[sum_{s=0}^{k} chi_1({s}) chi_2({s,k+1})]^2>0.   (7)

Reflection of the charge state identifies the two overlap factors.
The formula concerns projected transport and allows leakage from
the code, just as the original vacancy extension did.

## 3. A quantitative gap

Write delta for the charge gap, gamma for the core gap in (1), and

    p_min=min_j <chi_p|D_j|chi_p> >0,
    M=g(k-1)||h0||+||a||+||b||.

The spectral gap above the complete code (6) obeys

    gap(H_ext) >= delta p_min gamma /
                  [ (sqrt(delta)+sqrt(M))^2+p_min gamma ]. (8)

Proof. Subtract e_p. Call the positive second part of (5) V_int,
with ||V_int||<=M. For a normalized psi orthogonal to (6), write

    psi=x+eta,  x=|chi_p><chi_p| psi,  eta perpendicular to chi_p,
    E=<psi|H_ext-e_p|psi>.

Then ||eta||^2<=E/delta, and x=chi_p tensor phi with phi perpendicular
to the core kernel. Charge compression of V_int gives

    g sum_j p_j h0_j + a_1+b_k >= p_min H_core.

Consequently <x|V_int|x> >= p_min gamma ||x||^2. On the other hand,
the triangle inequality for V_int^(1/2) gives

    sqrt(<x|V_int|x>)
       <= sqrt(<psi|V_int|psi>)+sqrt(M)||eta||
       <= sqrt(E)(1+sqrt(M/delta)).

Use ||x||^2=1-||eta||^2>=1-E/delta and rearrange to obtain (8).
This argument also works if the last lower bound is negative, since
the rearranged inequality still follows. QED.

The bound respects the small tunneling gap of the charge state. It
does not claim a uniform lower bound as k grows, or as w approaches
zero. It makes the finite-scale assumption needed for later
inter-cell perturbation estimates explicit.

## 4. Checks and limitation

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 tagged_vacancy_checks.py

The example has q=3,k=3,h0=(I-swap)/2, a=.8|2><2|,
b=1.1|2><2|, and g=.7. Its complete core is the four-dimensional
symmetric space using only colors 0 and 1. The boundary-free core
would instead have dimension ten, so the boundary penalties have
a nontrivial effect. Compressed X,Y,Z on the two retained colors
give noncommuting endpoint operators.

The script directly assembles each term of (3), including the
oriented boundary penalties, at p=1,2 and w=.3,.5 with U=1.
Sector dimensions are 108 and 270. It checks complete ground
multiplicities, the code, both endpoints, the cut overlap, and (8).
All code/port/overlap errors are below 8e-16. The smallest measured
gap is 0.00040717926, above its lower bound 0.000037792674.
Outputs are in tagged_vacancy_checks.json.

The result removes the restriction to boundary-free frustration-free
cores from the vacancy extension. It does not establish a code with
many Clifford channels, select the ordered L B...B R sector against
arbitrary tag sequences, or prove a conformal phase after coupling
cells. The number of occupied sites can be selected within the
one-L/one-R ordered sector by the separate length filter below.
