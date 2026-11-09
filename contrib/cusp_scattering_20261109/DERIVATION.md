# Exact two-body solution and flux normalization

Produced by Codex (OpenAI GPT-6 family), 2026-11-09. The algebra below was
derived independently, followed by comparison with C175 after its arrival
from GitHub. All claims concern the specified effective Hamiltonian.

## 1. Hypotheses and relative-coordinate operator

[proved, reduction] Use the symmetric C173 effective Hamiltonian, gauge
transformed so the nearest-neighbour hop is `u=|w|>0` and the next-nearest
hop is `+1`. The blockade excludes separations zero and one. The two
interactions are `V2=-2`, `V3=-2/3`. Particle number is conserved in this
effective Hamiltonian. The one-particle dispersion is

\[
\epsilon(k)=2u\cos k+2\cos 2k.
\]

For `0<u<4`, its minima are at `+Q,-Q`, with `cos Q=-u/4`, and
`epsilon_min=-2-u^2/4`. In the ordered sector `x1<x2`, put
`r=x2-x1>=2` and
`Psi(x1,x2)=exp(i K(x1+x2)/2) f(r)`. Then

\[
(H_K f)_r=A_K(f_{r-1}+f_{r+1})+
B_K(f_{r-2}+f_{r+2})+V_r f_r,
\quad A_K=2u\cos(K/2),\quad B_K=2\cos K,
\]

with omitted terms whenever `r+d<2`. At `K=0`, `A=2u`, `B=2`, and the
essential-spectrum bottom is

\[
E_{\rm th}=2\epsilon_{\min}=-4-u^2/2.
\]

No identification of a many-body phase enters this reduction.

## 2. An exact bound-state polynomial

[proved] For any `u in (0,4)` there is a unique `p in (0,1)` satisfying

\[
F_u(p)=(1-p)(3-p)(1+p)^2-u^2p^3=0. \tag{1}
\]

Indeed, `(1-p)(3-p)(1+p)^2/p^3` decreases continuously from infinity to
zero. Its logarithmic derivative is

\[
-\frac1{1-p}-\frac1{3-p}+\frac2{1+p}-\frac3p<0,
\]

since the last two terms sum to `-(p+3)/(p(1+p))`.

Set

\[
s=-\frac{up}{1+p},\qquad
z=\frac{s+i\sqrt{4p-s^2}}2.
\]

Here `|z|^2=p<1`. The square root is real and nonzero: a failure would
require `u^2 >= 4(1+p)^2/p > 16`, contrary to `u<4`.

The exact bound-state energy is

\[
E_b=-2\left(p+\frac1p+\frac{u^2p}{(1+p)^2}\right). \tag{2}
\]

To verify this, substitute a bulk solution `z^r`. With `y=z+1/z`, the
characteristic equation is

\[
2y^2+2uy-(E+4)=0.
\]

The two roots `y,y*` have sum `-u`, precisely because
`s(1+1/p)=-u`. Their product gives (2).

The missing sites `r=0,1` are accounted for by two boundary conditions
on the free continuation of `f`:

\[
2u f_1+2f_0+2f_2=0,\qquad
2f_1+\frac23 f_3=0. \tag{3}
\]

For the two decaying modes `z^r,(z*)^r`, divide the boundary columns by
two and write

\[
b_1(z)=1+uz+z^2,\qquad b_2(z)=z+z^3/3.
\]

Their determinant, after division by the nonzero factor `(z*-z)/3`, is

\[
3+s^2+ups+p^2-4p.
\]

Substituting `s=-up/(1+p)` reduces its vanishing exactly to (1).
Therefore a nonzero decaying solution satisfies both boundaries.

The binding is strictly positive:

\[
E_{\rm th}-E_b=
\frac{(1-p)^2}{2p}
\left(4-\frac{u^2p}{(1+p)^2}\right)>0. \tag{4}
\]

The same reduction also proves uniqueness below the continuum in this
total-momentum sector. Below `E_th`, the two decaying characteristic roots
are a complex-conjugate pair. Every square-summable solution is their
linear combination, and the boundary determinant forces (1), whose root
in `(0,1)` is unique. This does not claim uniqueness across other total
momenta or exclude other multi-particle bound states.

At `u=1`, (1) becomes the integer polynomial
`p^4-3p^3-4p^2+2p+3=0`, so the bound energy is an algebraic number.

## 3. Exact wavefunction and pair mass

[proved] Put `C=conj(1+uz+z^2)`. A real bound wavefunction is

\[
f_r=\operatorname{Im}(Cz^r),\qquad r\ge2. \tag{5}
\]

The first condition in (3) is automatic, and the second follows from the
zero determinant. Its squared norm is the convergent geometric sum

\[
\mathcal N=\frac12\left[
|C|^2\frac{p^2}{1-p}
-\operatorname{Re}\frac{C^2z^4}{1-z^2}\right]. \tag{6}
\]

For the unilateral shift `T` on `r>=2`, define

\[
C_m=\langle T^m+(T^\dagger)^m\rangle=
\frac{ |C|^2\operatorname{Re}(z^m)p^2/(1-p)
-\operatorname{Re}[C^2z^{4+m}/(1-z^2)]}{\mathcal N}. \tag{7}
\]

The isolated eigenvalue continues smoothly for small total momentum `K`.
Since `H'_K=0` at zero, the Hellmann-Feynman expression for its curvature
has no wavefunction-derivative term. Thus

\[
\frac1{m_{\rm pair}}=E_b''(0)=-\frac u2 C_1-2 C_2. \tag{8}
\]

[numerical evaluation of exact formulas] At `u=1`, equations (1)-(8) give
`E_b=-4.519526534666399` and `1/m_pair=3.373597375631866`.
At `u=2`, they give `E_b=-6.122547914236252` and
`1/m_pair=2.250726229586402`. Independent finite-momentum diagonalization
reproduces the curvatures within `1e-5`; see `results.json`.

This mass is a two-body input. It does not decide atom-pair miscibility,
pair-pair scattering, cluster stability, or any central charge.

## 4. Independent audit of the scattering convention

[proved, construction; numerical limit checked] Let
`E=E_th+eta`, `eta>0` sufficiently small. There are four real relative
momenta, determined by

\[
q_l=\arccos\frac{-u/2+\sqrt{\eta/2}}2,
\qquad
q_h=\arccos\frac{-u/2-\sqrt{\eta/2}}2.
\]

Incoming modes have momenta `(q_l,-q_h)`, outgoing modes `(q_h,-q_l)`.
Their relative-coordinate velocities come from
`E(q)=4u cos q+4 cos 2q`:

\[
v(q)=-4u\sin q-8\sin2q.
\]

For general real `V2,V3`, form the boundary column

\[
\mathcal B(q)=
\begin{pmatrix}2u z+2-V_2z^2\\2z-V_3z^3\end{pmatrix},
\qquad z=e^{iq}.
\]

With `B_in` and `B_out` the corresponding two-column matrices, the
amplitude matrix is `S=-B_out^{-1} B_in`. Normalizing each wave to unit
probability flux gives

\[
S_F=\operatorname{diag}\sqrt{v_{\rm out}}\,
S\,\operatorname{diag}(-v_{\rm in})^{-1/2}. \tag{9}
\]

This matrix is unitary. Define
`delta=sqrt(eta/(8-u^2/2))`. The physical threshold length matrix is

\[
\mathcal A_F=\lim_{\eta\to0^+}\frac{S_F+I}{2i\delta}. \tag{10}
\]

In contrast, directly writing the zero-energy solution as
`f_r=sum_sigma (r beta_sigma+alpha_sigma) exp(i sigma Q r)` and setting
`alpha=-A_raw beta` yields a matrix with a small imaginary diagonal.
The finite-energy calculation (9)-(10) independently checks, numerically
at the three parameter choices saved in `results.json`, that

\[
\mathcal A_F=(A_{\rm raw}+A_{\rm raw}^\dagger)/2. \tag{11}
\]

This verifies C175's convention correction in those cases using probability
flux, rather than assuming that discarding an anti-Hermitian part is
harmless. A symbolic general proof of (11) is not supplied here.
At `u=1,V2=-2,V3=-2/3`, in the channel ordering above,

\[
\mathcal A_F=\begin{pmatrix}
10.5&10.125-1.323269309954i\\
10.125+1.323269309954i&10.5
\end{pmatrix}.
\]

Its eigenvalues are `0.288894934109` and `20.711105065891`.
The raw matrix has additional diagonal entries `+i/sqrt(60)` and
`-i/sqrt(60)`. This order of channels can be reversed relative to C175
without changing any physical eigenvalue.

The script checks (11) at decreasing `eta`, as well as unitarity of (9),
both for the candidate interactions and the zero-interaction control.
Taking the Hermitian part of the finite-eta expression removes its
unitarity-required order-delta correction and makes convergence faster.

## 5. Scope of the conclusion

[proved] An unbound two-component gas is not the lowest two-particle
threshold of the specified effective Hamiltonian at the symmetric point.

[open] A two-body bound state does not imply a particular finite-density
central charge. Larger clusters, overlapping pairs, and a possible
intermediate-density two-component phase need additional analysis.
The candidate density around 0.17 is not a controlled dilute limit:
its product with the large scattering length is about 3.5. Moreover, the
microscopic qubit Hamiltonian has no exact particle-number symmetry;
its comparison with this effective model requires the separate
finite-coupling reduction.
