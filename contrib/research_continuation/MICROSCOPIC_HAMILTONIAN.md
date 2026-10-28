# Microscopic bond for the proposed counterexample

## Status and what is specified

This is a **parameterized construction formula**, not a numerically
instantiated counterexample. The finite-controller gates inside the core
are specified through compiler arguments, but their complete transition
table and matrices have not been expanded. No violating prime, numerical
local dimension, or positive perturbation strength has been selected.
The central-charge prediction uses the infrared assumptions in
[the complete candidate](PROPOSED_COUNTEREXAMPLE.md) and [its audit](IR_AUDIT.md).

Author: OpenAI Codex, GPT-6 family, 2026-10-28 UTC. An AI research worker
cross-checked these formulas through read-only comparison with the cited
repository notes. The half-density chemical potential was derived directly
and independently checked by the coordinating AI. No new full-matrix
computation or human review is claimed.

## 1. Core input and the on-site Hilbert space

For a prime \(p\), let \(u=p^2+p+1\), \(q=2u\), and \(N=(p+1)u\).
The \(N\) graph edges label flavours. The proposed initialized core has
local dimension \(d_c\le Cq(1+\log_2q)\), selected length \(k\ge3\), and
\[
H_c(k)=\sum_{j=1}^{k-1}(h_c)_{j,j+1}+(a_c)_1+(b_c)_k,
\qquad h_c,a_c,b_c\ge0.
\]
Its asserted complete zero-energy code has dimension \(2^N\). Flagged
endpoint matrices have compressions
\[
F_{L,a}\mapsto\alpha_L f_a^\dagger,\qquad
F_{R,a}\mapsto\alpha_R f_a^\dagger P,\qquad
\alpha_L=\frac{p_L}{2\sqrt{uN}},\quad
\alpha_R=\frac{p_R[1-1/(2N)]}{2\sqrt{uN}}.                 \tag{1}
\]
Here \(p_L,p_R>0\) are history readout fractions and \(P\) is logical
fermion parity. See [the core construction](../torus_bottleneck/EXISTENTIAL_CORE.md).
For an edge \(a=(\alpha,\beta)\), the readout acts by the species matrix
unit \(|1,\alpha\rangle\langle0,\beta|\) on the output record, gated by
READY and absence of the interpreter-head overlay; it vanishes otherwise.
Thus the microscopic port is an on-site matrix, not an operator acting
directly on an entire encoded cell.

Set \(b=k+1\). Each *microscopic* site has Hilbert space
\[
\mathcal H=\mathbb C^4_{\rm sel}\otimes
\left[\mathbb C^3_{\rm vac}\oplus
 (\mathbb C^3_{\rm tag}\otimes\mathbb C^{d_c}_{\rm core}
 \otimes\mathbb C^2_{\rm filter})\right]
 \otimes\mathbb C^2_A\otimes\mathbb C^2_B,
\qquad d=96d_c+48.                                      \tag{2}
\]
Selector letters are \(\#,L,X,R\); vacancies are \(H_0,H_1,H_2\);
occupied tags are \(L_w,B_w,R_w\). All added factors are species-neutral.
Occupied letters retain their core record's species representation.

## 2. Explicit selector and internal wrapper

Give selector letters arrows
\[
\#=(R,L),\quad L=(L,L),\quad X=(L,R),\quad R=(R,R).
\]
The eight allowed adjacent pairs are

    #L, #X, LL, LX, X#, XR, R#, RR.

Let \(Q\) project diagonally onto their complement, and put
\[
M=|\#\rangle\langle\#|,\qquad
P_{\rm cut}=|X\#\rangle\langle X\#|+|R\#\rangle\langle R\#|,
\qquad P_{\rm int}=I-Q-P_{\rm cut}.                      \tag{3}
\]
Thus \(P_{\rm int}\) projects onto the six pairs \(\#L,\#X,LL,LX,XR,RR\).
The wrapper letters have arrows \(H_j:j\to j\), \(L_w:0\to1\),
\(B_w:1\to1\), \(R_w:1\to2\). Write \(\Pi_{\rm bad}\) for mismatched
wrapper arrows, \(n\) for the vacancy projector, and \(P_{L_w},P_{R_w}\)
for the occupied endpoint-tag projectors. For core/filter labels \(\xi\),
\[
T_v=\sum_\xi\left(
 |L_\xi H_1\rangle\langle H_0L_\xi|
 +|B_\xi H_1\rangle\langle H_1B_\xi|
 +|R_\xi H_2\rangle\langle H_1R_\xi|\right).
\]
Tensor the initialized core with the two-state length filter of pitch
\(\theta=2\pi/(k+2)\). Its resulting positive bulk and boundary matrices
\(h_0,a_0,b_0\) are defined in
[length_filter.md, equation (8)](../../codex/work_C8/length_filter.md).
Lift \(h_0\) to occupied pairs \(L_wB_w,L_wR_w,B_wB_w,B_wR_w\), giving
\(\widetilde h_0\); lift the boundary matrices to their respective tags,
giving \(a_L,b_R\). The explicit wrapper bond is
\[
\begin{aligned}
w={}&6U\Pi_{\rm bad}+U(n\otimes I+I\otimes n-2n\otimes n)
       -v(T_v+T_v^\dagger)+\widetilde h_0\\
 &+(a_L-2UP_{L_w})\otimes I+I\otimes(b_R-2UP_{R_w}).
\end{aligned}                                           \tag{4}
\]
For example take \(U=1\), \(v=1/[2(k+1)]\), and the filter strength
prescribed in length_filter.md. This wrapper has complete selected open
codes at lengths \(b,b+1\), with endpoint occupancy probabilities
\(r_1,r_2>0\). Choose \(t,\mu_0,K_{\rm sel}\) using
[selector_wrapper.md, equations (6)–(7)](../../codex/work_C8/selector_wrapper.md).
Those explicit choices select these two lengths with a positive gap above
the entire unperturbed ring ground space. They require only finite cell
quantities; in particular \(K_{\rm sel}=4t+9U\).

## 3. Two helix tracks and physical port matrices

Choose \(\phi=1/b\), and define
\[
|+_x\rangle=(|0\rangle+e^{ix}|1\rangle)/\sqrt2,\quad
h_\phi=|s_\phi\rangle\langle s_\phi|,\quad
|s_\phi\rangle=(|01\rangle-e^{-i\phi}|10\rangle)/\sqrt2,
\quad a=I-|+_0\rangle\langle+_0|.
\]
The additional \(A,B\) tracks have \(h_\phi\) on internal bonds, with
\(A\) pinned on the first site and \(B\) on the last. Let
\(u_1=|+_{-(b-1)\phi}\rangle\), \(u_2=|+_{-b\phi}\rangle\), and
\[
E_L=\frac{r_1^{-1}(I-|u_2\rangle\langle u_2|)
                 +r_2^{-1}(I-|u_1\rangle\langle u_1|)}{\sin^2(\phi/2)},
\qquad E_R=E_L^*.                                       \tag{5}
\]
Lift \(F_{L,a}\) to tag \(L_w\) and \(F_{R,a}\) to tag \(R_w\), with
identity on the original filter and zero on other wrapper letters.
Call these lifts \(\widetilde F_{L,a},\widetilde F_{R,a}\), and define
\[
A_a=\widetilde F_{L,a}\otimes I_A\otimes E_L,\qquad
B_a=\widetilde F_{R,a}\otimes E_R\otimes I_B.              \tag{6}
\]
At both selected lengths these have exactly the compressions (1).
See [helix_equalizer.md](../../codex/work_C8/helix_equalizer.md).

## 4. The single repeated microscopic bond

Identities on unused factors are suppressed. Track operators are grouped
by physical site, so every term acts on two physical sites only:
\[
\begin{aligned}
h_p(s)={}&K_{\rm sel}Q+\frac{\mu_0}{2}(M\otimes I+I\otimes M)
 -t(|LX\rangle\langle XR|+\mathrm{h.c.})+P_{\rm int}w\\
&+P_{\rm int}(h_\phi^A+h_\phi^B)
 +P_{\rm cut}(a_B\otimes I+I\otimes a_A)\\
&+s\left[v_{\rm hop}+v_{\rm move}
             +\tfrac12(z\otimes I+I\otimes z)\right].
\end{aligned}                                           \tag{7}
\]
Here \(a_A,a_B\) mean the pin \(a\) on their respective tracks. The
normalized payload interaction is
\[
v_{\rm hop}=-\frac{P_{\rm cut}}{\alpha_L\alpha_R}
 \sum_{a=1}^N t_a(B_a\otimes A_a^\dagger+B_a^\dagger\otimes A_a),
\quad k_a=\frac\pi4+\frac{a\pi}{12(N+1)},\quad
t_a=\frac1{2\cos k_a}.                                  \tag{8}
\]
For marker transport put \(U_b=\operatorname{diag}(1,e^{ib\phi})\),
\(C=e^{-ib\phi}|H_2\rangle\langle H_0|\otimes U_b\otimes U_b\), and
\[
v_{\rm move}=-\frac{J_\ell}{A_b^2\beta_{\rm wrap}\cos^{2b}(\phi/2)}
 \left[|R\#\rangle\langle\#L|\otimes(C\otimes I)+\mathrm{h.c.}\right].
                                                               \tag{9}
\]
Here \(A_b\) is the positive selector-head overlap and \(\beta_{\rm wrap}\)
is the positive vacancy-conversion overlap in
[selector_wrapper.md, Section 5](../../codex/work_C8/selector_wrapper.md).
The conversion acts on the transferred left site. Its projection
exchanges neighboring short/long cells with amplitude \(-J_\ell\),
and acts as the identity on their logical states.

Let \(q_\uparrow\) be the **single-site** upper-species charge on occupied
core records, zero on vacancies. Its global sum restricts to
\(\sum_{j,a}n_{j,a}+4u\,M_{\rm tot}\), with
\(M_{\rm tot}=\sum_i M_i\) the total number of cells. Therefore the
single-site field
\[
z=q_\uparrow+(\zeta-4u)M                               \tag{10}
\]
projects, after summing, to logical occupation plus \(\zeta\) per cell.
The periodic microscopic Hamiltonian is
\[
H_L=\sum_{i=1}^L h_p(s)_{i,i+1},\qquad L+1\equiv1.       \tag{11}
\]
It is one-site translation invariant and nearest-neighbor. Every term
preserves the total species charges and total marker number.

## 5. An explicit first-order chemical potential

Take \(J_\ell=1/2\). A payload flavour with dispersion
\(1-2t_a\cos k\) has ground energy per cell
\[
e_a=\frac1\pi\int_0^{k_a}(1-2t_a\cos k)\,dk
    =\frac{k_a-\tan k_a}{\pi},\qquad e_{\rm pay}=\sum_a e_a.
\]
Choose \(\zeta=1/\pi-e_{\rm pay}\). At long-cell fraction \(\nu\), there
are \(L/(b+\nu)\) cells to leading order. The first-order ground-energy
correction per physical site, in units of \(s\), is
\[
\frac{e_{\rm pay}+\zeta-(2J_\ell/\pi)\sin(\pi\nu)}{b+\nu}
 =\frac{1-\sin(\pi\nu)}{\pi(b+\nu)}.                     \tag{12}
\]
It is globally minimized at \(\nu=1/2\). The cell-coordinate velocities
are \(2J_\ell=1\) for the length mode and \(\tan k_a\in(1,\sqrt3)\) for
the payload. Half density is a first-order statement; microscopic
corrections can shift it. The proposed stable interval contains it.

## 6. What is still required for a literal matrix and a counterexample

Schematically, the core uses positive local propagation blocks
\[
\frac J2\left[|a\rangle\langle a|+|b\rangle\langle b|
 -|b\rangle\langle a|\otimes U_{ab}
 -|a\rangle\langle b|\otimes U_{ab}^\dagger\right],
\]
plus clock legality, initialization, and FAIL penalties. Some core
instructions and initialization pins are boundary operators \(a_c,b_c\);
they are incorporated into the repeated wrapper bond through (4).
The clock rules are explicit, but the compiler's full finite controller
defining the \(U_{ab}\) has not been expanded. That expansion, selected
\(p,k\), and evaluation of the finite overlaps are needed to print every
entry of (7) and give a numerical \(d\). A numerical nonzero \(s\) remains
unchosen. This formula does not fill in those omitted data.

All coefficients in (7) must be fixed before \(L\to\infty\). The first-order
projection has \(N\) free fermion channels and one length mode. Concluding
that (11) at fixed \(s>0\) has \(c=N+1\) uses local elimination of massive
machinery and stability of the resulting Gaussian theory. The formula
alone does not establish those thermodynamic claims or a verified violation.
