# What a simple physical detuning can accomplish

## 1. Restriction inside the four-term family

For every real a,b,u,h in the original displayed family, the standard
boson-frame hopping is imaginary: w=i(a-b). Changing these four parameters
therefore retains the exact Sigma symmetry in the
[locking ledger](../v4_locking_operators_20261118/LEDGER.md). As long as a
balanced, symmetry-preserving two-valley description exists, its dressed
centres sum to pi. The quartet momentum stays exactly zero modulo 2pi.
A field change can move the density; it cannot detune this particular
centre-momentum resonance. A change of state that breaks the symmetry
requires a new operator ledger and is not evidence that both modes survive.

In the controlled cusp the leading separation-three source is
2Delta^3/(3V^2). Neither kappa nor the leading scaled chemical potential
enters this coefficient. For Delta,V>0 it cannot be canceled by those
controls. Delta=0 removes the generated t2 at the same time and reduces
the kinetic mechanism to one nearest-neighbour band. This is a restricted
statement about the leading source, **not** a proof that the fully matched
quartet amplitude is nonzero everywhere. Section 5 of AMPLITUDE.md exhibits
one finite state with a genuine leading cancellation.

Nor can one independently weaken V2 or V3 within this controlled family:
V2/t=-2 and V3/t=-2/3 are fixed by the same Delta,V that generate t2.
Replacing them by freely adjustable repulsions changes the effective
model without identifying an allowed microscopic control.

## 2. An exact vacuum binding check

The earlier [C175 derivation](../../codex/work_C175/PROOF.md) already found
an opposite-valley bound pair throughout 0<kappa<4 and explained why it
invalidates the unbound dilute-gas assumption. Here is a compact bound-state
equation, independently derived and evaluated as a detuning check. It is
not claimed that the existence of the bound pair is new.

At total momentum pi, write the relative-coordinate wavefunction in a
centre-of-mass gauge. For r>=2 its real hopping coefficients are 2kappa
for r to r+1 and 2 for r to r+2, with boundary potentials -2 at r=2
and -2/3 at r=3. A decaying solution is a combination of z1^r,z2^r.
Let q=z1 z2 and p=z1+z2. The reciprocal bulk polynomial gives

\[
p=-\kappa q/(1+q),\qquad
E_b=-2\left(q+q^{-1}+\kappa^2q/(1+q)^2\right).
\]

The determinant of the two boundary equations factors by z1-z2 and reduces
to

\[
\boxed{(1-q)(3-q)(1+q)^2=\kappa^2q^3,\quad 0<q<1.}\tag{1}
\]

There is exactly one such q: the left side divided by q^3 decreases
strictly from infinity to zero, since its logarithmic derivative is
-1/(1-q)-1/(3-q)+2/(1+q)-3/q<0. For 0<kappa<4 its two roots are a
complex conjugate pair with modulus sqrt(q)<1. One way to check the latter
inequality is to note that at kappa=4, q>27/50, while
(sqrt(13)-2)/3<27/50. The inequality q>(sqrt(13)-2)/3 is exactly p^2<4q.

The binding below the two-atom continuum is

\[
(-4-\kappa^2/2)-E_b
 ={(1-q)^2(3q^2+4q-3)\over2q^3}>0.
\]

It begins as kappa^4/32. At kappa=1 it is exactly specified by (1) and
numerically equals 0.0195265346664 in units t. The root and an independent
512-site relative-coordinate calculation agree to rounding error there.
At kappa=1/2 their difference is 1.9e-9, the finite-box tail.

Thus merely varying the imaginary hopping inside 0<kappa<4 does not remove
this vacuum competitor in the leading model. This does **not** prove a
finite-density gap or a specific c=1 molecular liquid: higher bound
clusters, sector selection, density effects and the number-breaking source
still matter. C175/C191's finite-density interpretations must not be
replaced by a vacuum conclusion.

## 3. One extra literal monomial breaks the resonance

Adding v X_j X_(j+1) is enough; adding both XX and YY is unnecessary for
this purpose. In the displayed basis there are then five literal monomials
and five independent real controls before removing a common scale.
After a harmless uniform z rotation to make the pair coefficient real,

\[
\Delta'=\sqrt{(a+b)^2+v^2},\qquad
w=v+i(a-b),\quad V=4u,\quad\mu=4u-2h.
\]

The physical pairing phase is included in that rotation; it is not
silently set to zero while the rest of the bond is held in another frame.
Use t=Delta'^2/(2V), r=Re(w)/t and kappa=Im(w)/t. The leading band is

\[
e(k)/t=2r\cos k+2\kappa\sin k-2\cos2k.
\]

At r=0, sin(k_A)=sin(k_B)=-kappa/4. For small r,

\[
q_4=2(k_A+k_B-\pi)=-{4\kappa r\over16-\kappa^2}+O(r^3),
\quad
|e_A-e_B|/t=4|r|\sqrt{1-\kappa^2/16}+O(|r|^3).
\]

The pi in the first expression is understood modulo 2pi; the code uses
centres near -pi and zero, so it evaluates 2(k_A+k_B+pi). At kappa=1,
the offset is sqrt(15)|r| whereas the quartet mismatch is only -4r/15.
Detuning the resonance therefore comes with a much larger valley-energy
bias, whose many-body susceptibility cannot be ignored.

The *same physical change* generates another projected pair source:

\[
H_{\rm pair}/t=\epsilon\left[{4\over3}B_3^\dagger-2r B_2^\dagger
 +\mathrm{h.c.}\right]+\cdots,
\qquad B_d^\dagger=P\sum_jb_j^\dagger b_{j+d}^\dagger P.
\]

Keeping only the shifted band would omit a physical term. This new source
need not be resonant at generic centres, but its products and density
harmonics belong in the off-symmetry matching problem. No cancellation of
the complete quartet source is claimed by comparing just one form factor.

## 4. A rational diagnostic point and its costs

For a concrete bounded test, take the coefficients fixed in L:

\[
H=\sum_j[{33\over64}XY+{31\over64}YX+4ZZ
          +{1\over512}XX]-{127\over16}\sum_jZ_j.
\]

This is a d=2, Hermitian, nearest-neighbour, one-site-invariant physical
Hamiltonian with epsilon approximately 1/16. Its leading dimensionless
parameters are kappa=0.9999962, r=0.06249976, mu-bar=-2.0000153. The
bounded calculation gives:

| Leading quantity | Value |
|---|---:|
| Valley minima | -2.8930075, -0.2569213 |
| Minimum energies divided by t | -2.3710924, -2.1290340 |
| Quartet momentum mismatch | -0.01667230 |
| Valley offset divided by t | 0.24205836 |
| Separation-two source divided by t | -0.00781249 |
| Separation-three source divided by t | 0.08333349 |

Both bare minima lie below the chosen bare chemical potential. The
momentum-pi two-body Hamiltonian is exactly independent of Re(w), so its
bound level can be compared without re-solving it: at this point it lies
0.2226623t **above** twice the lower atomic minimum. This only removes
the energetic advantage of that particular pi dimer relative to lower-valley
atoms. It does not exclude a bound state at a shifted total momentum,
other clusters, valley polarization or finite-density pairing.

Simplicity is a drawback, not a claimed success here. After removing the
constant and one common positive scale, the literal nonzero Pauli
coefficient ratio is (127/16)/(1/512)=**4064**. Without the added XX term
the same displayed four-term point has spread 508/31 approximately 16.39.
V4 itself has four terms and spread 4. The small extra coefficient is a
microscopic hierarchy, distinct from t or a small emergent velocity.
Although the chosen fractions are uncomplicated, the detuning is not a
demonstrated robust phase region or an improvement over V4's simplicity.

## 5. Why this is not a stable-region certificate

On the symmetry line the full conditional test is

\[
\min_{(r,s)\in\mathbb Z^2\setminus\{0\}}
 [a_G r^2+b_Gs^2-2h_Grs]>1.
\]

This calculation has not assigned a_G,b_G,h_G from the microscopic model.
No arbitrary Luttinger matrix is presented as such an assignment.
Consequently we cannot test that inequality for V4 or the detuned point.

Off the symmetry line, reapply the general bosonic momentum rule with the
actual dressed centres and densities:
q=n_A k_A+n_B k_B+2pi(m_A rho_A+m_B rho_B), even total charge, and the
remaining microscopic symmetries. Generic nonresonance is not a numerical
stability margin. A nonzero detuning must interrupt a relevant coupling
before strong coupling, while both compressibilities/relative stiffnesses
remain positive and both valleys remain populated. These requirements
need a finite-density matching calculation, not the vacuum binding energy
or the displayed band minima. This round identifies a minimal literal
control and tests its consequences without claiming those missing facts.
