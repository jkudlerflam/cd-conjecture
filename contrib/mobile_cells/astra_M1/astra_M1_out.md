The ring identification is exact, but factorization of the eigenstates is generally not. In the fluid regime, a length boson and a payload CFT give a useful infrared description only after checking the allowed couplings. For the toy, a reflection symmetry forbids uniform squeezed-spin dimerization even when u ≠ t. This makes a c = 2 phase plausible and stable at weak XX exchange, away from commensurate instabilities and fluid endpoints.

1. Exact ring kinematics

Take W = 2M > 0 and label sites by 0,…,L−1. A legal basis state is specified by

\[
0\le x_1<\cdots<x_W<L,\qquad
\tau\in\{X,\#\},\qquad
(s_1,\ldots,s_W).
\]

Here τ is the species of the first wall, subsequent species alternate, and s_a is its spin. Thus, as a vector space,

\[
\mathcal H_W\simeq
\mathbb C^{\binom LW}\otimes\mathbb C^2_\tau
\otimes(\mathbb C^2)^{\otimes W}.
\]

The extra two-state factor cannot be omitted. For W = 0, the legal sector consists of the all-L and all-R states.

In unwrapped coordinates the cyclic identification is

\[
(x_1,\ldots,x_W;\tau;s_1,\ldots,s_W)
\sim
(x_2,\ldots,x_W,x_1+L;\bar\tau;s_2,\ldots,s_W,s_1).
\tag{1}
\]

This identification, rather than an unconstrained product of two periodic chains, is the ring geometry.

Let C cyclically move the last wall label to the first, including \(\tau\mapsto\bar\tau\). Physical translation T acts as

\[
T|x_1,\ldots,x_W;\chi\rangle=
\begin{cases}
|x_1+1,\ldots,x_W+1;\chi\rangle,&x_W<L-1,\\
|0,x_1+1,\ldots,x_{W-1}+1;C\chi\rangle,&x_W=L-1.
\end{cases}
\tag{2}
\]

In particular, crossing the coordinate seam rotates the squeezed chain. Since W is even, \(C^W=1\).

Suppose a factorized component has \(C\chi=e^{iQ_w}\chi\). Use the convention \(T|\Psi\rangle=e^{iP}|\Psi\rangle\), with coordinate plane waves \(e^{-ikx}\). Equation (1) requires its charge wavefunction to obey

\[
\psi(x_2,\ldots,x_W,x_1+L)
=e^{-iQ_w}\psi(x_1,\ldots,x_W).
\tag{3}
\]

At u = t and J = J_z = 0, a fermionic Slater determinant therefore has

\[
e^{ik_aL}=(-1)^{W-1}e^{iQ_w},
\qquad
k_a=\frac{2\pi I_a+Q_w+\pi(W-1)}{L},
\tag{4}
\]

and

\[
P=\sum_a k_a\pmod{2\pi},\qquad
E=-2t\sum_a\cos k_a+\frac{\mu W}{2}.
\tag{5}
\]

The parity factor in (4) converts hard-core particles into fermions. Multiplying the quantization conditions verifies \(PL\in2\pi\mathbb Z\).

This is the same mechanism as Ogata–Shiba factorization. At infinite Hubbard repulsion, particle order is conserved and a cyclic spin permutation supplies the charge boundary twist. Differences here are the alternating wall species, the extra arrow-orientation register, and bosonic rather than electronic microscopic exchange statistics. At unequal t,u, or with exchange present, the kinematic identification survives but the free-charge factorization generally does not.

There is a useful cell version. Choose a marker at x and write

\[
d=(\ell_0,r_0;\ldots;\ell_{M-1},r_{M-1}),\qquad
\sum_j\ell_j=L,
\]

where r_j specifies the head position. Let S shift cell lists to the left. Then

\[
(x;d;s)\sim(x+\ell_0;Sd;Ss).
\tag{6}
\]

For a momentum eigenfunction
\(\Psi=e^{-iPx}f(d)\chi(s)\), define the payload momentum by
\(\chi(Ss)=e^{iQ}\chi(s)\). Equation (6) gives

\[
f(Sd)=e^{i(P\ell_0-Q)}f(d).
\tag{7}
\]

Using the unwrapped marker center of mass R instead of x removes the configuration-dependent phase. Under (6), \(R\mapsto R+L/M\), so

\[
f(Sd)=e^{i(PL/M-Q)}f(d).
\]

If the length factor has cyclic momentum \(q_\ell\),

\[
\boxed{\frac{PL}{M}=q_\ell+Q\pmod{2\pi}.}
\tag{8}
\]

For static cells of length b this reduces to \(bP=Q\). For coupled eigenstates, (7) is a boundary condition on the full wavefunction, potentially expanded in several payload momenta. A separately conserved payload momentum is not guaranteed.

2. Continuum theory and allowed couplings

First use the cell coordinate y. Near a b,b+1 crossing, the long-cell occupation has the expansion

\[
\eta(y)=\nu-\frac{\partial_y\phi}{\pi}
+\sum_{m\ne0}A_m e^{2im(\pi\nu y-\phi)}.
\tag{9}
\]

Our normalization is

\[
H_\ell=\frac{v_\ell^{\rm cell}}{2\pi}\int dy\,
\left[K_\ell(\partial_y\theta)^2+
K_\ell^{-1}(\partial_y\phi)^2\right],
\qquad
\dim e^{2im\phi}=m^2K_\ell.
\tag{10}
\]

Here the leading projected XX length chain has \(K_\ell=1\). Conservation of cell number, together with fixed physical length, conserves the long-cell number in this projection. It forbids charge-changing \(e^{in\theta}\) perturbations. In the general construction the payload is neutral under this U(1).

Write a payload lattice operator as

\[
h_j^p=\text{constant}+\sum_\alpha
c_\alpha e^{iq_\alpha j}\mathcal O_\alpha(y)+\cdots .
\]

All local mixed perturbations are organized by

\[
e^{i(2\pi m\nu+q_\alpha)y}
e^{-2im\phi}\mathcal O_\alpha,
\tag{11}
\]

and their derivative dressings and Hermitian conjugates. Products of payload fields are included through their OPE channels. Their dimensions and translation selection rule are

\[
\boxed{x=m^2K_\ell+\Delta_\alpha+n_\partial,\qquad
2\pi m\nu+q_\alpha=0\pmod{2\pi}.}
\tag{12}
\]

Only operators invariant under the actual internal symmetries survive. Reflection and time reversal further select the allowed real combinations.

Some important members of this classification are

| Perturbation | Dimension | Translation condition |
|---|---:|---|
| Pure length umklapp \(\cos(2m\phi+\gamma)\) | \(m^2K_\ell\) | \(m\nu\in\mathbb Z\) |
| Uniform payload perturbation \(\mathcal O_\alpha\) | \(\Delta_\alpha\) | \(q_\alpha=0\) |
| Density coupling \(\partial_y\phi\,\mathcal O_\alpha\) | \(1+\Delta_\alpha\) | \(q_\alpha=0\) |
| Current coupling \(\partial_y\theta\,\mathcal O_\alpha\) | \(1+\Delta_\alpha\) | \(q_\alpha=0\), plus time-reversal constraints |
| Density times payload stress tensor | 3 | Always |
| Charge-current times a payload current | 2 | If internal symmetries allow it |

The last line matters. A permitted marginal current coupling need not flow to zero. For two bosons it can give a mixed Gaussian theory with c = 2, without separation into the original length and payload variables.

Likewise, “the fixed-length payload chain is critical” is insufficient by itself. Length fluctuations can generate symmetry-allowed relevant payload perturbations. The payload critical point must be stable under those perturbations, or the required parameters must be tuned.

The product fixed point is stable when all allowed mixed perturbations are irrelevant or marginally irrelevant, the payload remains critical, and the length sector remains a one-component fluid. Its central charge is

\[
c=1+c_p.
\tag{13}
\]

With unequal velocities this is a product of two CFT Hamiltonians, with no common relativistic velocity.

For a payload with one site per cell, q = 0 and π give different resonances.

For q = 0, a length harmonic survives when \(m\nu\in\mathbb Z\). At reduced rational filling \(\nu=p/q\), the first pure umklapp has dimension \(q^2K_\ell\).

For a staggered payload operator,

\[
m\nu\in\mathbb Z+\frac12,\qquad
x=m^2K_\ell+\Delta_\pi.
\tag{14}
\]

At half filling, the leading mixed perturbation therefore has dimension \(K_\ell+\Delta_\pi\). Examples are

\[
\begin{array}{c|c|c}
\text{payload}&\Delta_{\rm dimer}&
x_{\rm mixed}\text{ at }\nu=\tfrac12\\ \hline
\text{XX}&1&K_\ell+1\\
\text{Heisenberg}&\tfrac12&K_\ell+\tfrac12
\end{array}
\]

At \(K_\ell=1\), the first is marginal and requires an RG calculation; the second is relevant if its coefficient is nonzero. More generally, a staggered resonance exists at reduced p/q only for even q, with smallest \(m=q/2\).

For the toy, there are two spins per cell. Write

\[
\ldots,\#_j,X_j,\#_{j+1},X_{j+1},\ldots .
\]

A spin operator at squeezed momentum π is at cell momentum 0. Thus cell translation alone would allow spin dimerization. The additional microscopic symmetry is decisive.

For all t,u, the Hamiltonian has the reflection

\[
\mathcal I:
\quad i\mapsto-i,\qquad L\leftrightarrow R,
\qquad \#_\sigma\mapsto\#_\sigma,\quad
X_\sigma\mapsto X_\sigma .
\tag{15}
\]

It preserves legality and maps each hopping term to its Hermitian conjugate with the same amplitude. It also preserves the chemical potential and exchange terms. Other symmetries include physical translation, conserved wall number, spin U(1), spin reversal, and time reversal. Spin symmetry enlarges to SU(2) when \(J_z=J\).

Arrow complementation \(L\leftrightarrow R,\ \#\leftrightarrow X\) exchanges t and u. In the legal sector it is an additional symmetry at t = u, since \(N_\#=N_X\).

Reflection (15), centered on a marker, acts on squeezed spins as site-centered reflection. Hence

\[
\epsilon_s=\sum_a(-1)^a h^s_{a,a+1}
\quad\longmapsto\quad-\epsilon_s.
\tag{16}
\]

Consequently,

\[
\boxed{\text{Uniform squeezed-spin dimerization is forbidden even for }u\ne t.}
\]

The box projection provides a direct check. A head has equal probabilities of touching its left and right markers, giving

\[
H_s^{\rm proj}
=\sum_j g(\ell_j)
\left[h^s_{\#_j,X_j}+h^s_{X_j,\#_{j+1}}\right].
\tag{17}
\]

The staggered components of the two bonds cancel at leading order. Their smooth remainder is a derivative of the dimer field.

For weak XX exchange, the spin boson has \(K_s\simeq1\). The potentially dangerous fields are

- staggered \(S^z\), dimension \(K_s\), forbidden by spin reversal;
- transverse spin fields, forbidden in the Hamiltonian by spin U(1);
- dimerization, dimension \(K_s\), forbidden by (15);
- spin umklapp, dimension \(4K_s\), irrelevant near \(K_s=1\).

A reflection-even mixed dimer perturbation can contain a reflection-odd length factor. A smooth example is

\[
(\partial_y\delta\eta)\epsilon_s,\qquad x=2+K_s,
\tag{18}
\]

which is irrelevant near XX. At commensurability, an appropriately transforming sine vertex can also multiply \(\epsilon_s\), with dimension \(m^2K_\ell+K_s\) and \(m\nu\in\mathbb Z\). Its coefficient may vanish in a particular projection.

Thus the one-spin-per-cell half-filling instability cannot be transferred directly to this toy. In the weakly coupled b,b+1 fluid, nontrivial rational fillings have denominator at least two, and these candidate commensurate perturbations are irrelevant near \(K_\ell=K_s=1\).

Exceptions include relevant commensurate locking, marginally relevant couplings, strong exchange, spontaneous spin order or dimerization, phase separation, and fluid endpoints. At an endpoint a velocity can vanish and the length transition need not have z = 1.

3. Finite rings and zero modes

A normalization issue is useful here. Put

\[
a=b+\nu,\qquad \rho_M=\frac1a.
\]

The physical cell-density boson and the long-cell boson obey

\[
K_M=\frac{K_\ell}{a^2},\qquad
v_\ell=a\,v_\ell^{\rm cell}.
\tag{19}
\]

Indeed, at fixed L,

\[
N_{\rm long}-\nu M=L-aM=-a\,\delta M,
\]

so their number-mode energies agree:

\[
\frac{(L-aM)^2}{4K_\ell}
=\frac{\delta M^2}{4K_M}.
\tag{20}
\]

If the charge is instead W = 2M, use \(K_W=4K_M\). These conventions must not be mixed.

In a stable product phase, the spectrum has the form

\[
E=e_\infty L+\frac{2\pi}{L}
\left[
v_\ell\left(x_\ell-\frac1{12}\right)
+v_p\left(x_p-\frac{c_p}{12}\right)
\right]+o(L^{-1}),
\tag{21}
\]

where

\[
x_\ell=
\frac{(M-\rho_*L)^2}{4K_M}
+K_M\left(D+\frac{\Phi}{2\pi}\right)^2
+n_R+n_L,
\tag{22}
\]

and

\[
x_p=h_a+\bar h_a+N_R^p+N_L^p.
\tag{23}
\]

Here e∞ includes the chemical potential and is evaluated at its thermodynamic optimum \(\rho_*\). Otherwise one must retain the linear chemical-potential term. The payload velocity in (21) is measured in physical distance: for a payload chain of \(N_p\) sites and squeezed velocity \(v_p^{\rm sq}\),

\[
v_p=\frac{L}{N_p}v_p^{\rm sq}.
\]

The allowed D, flux Φ, and payload sectors are glued by the ring identification. They cannot be chosen independently. In a factorized ordered-particle description, Φ is the squeezed momentum, with statistics and orientation encoded in the allowed current sectors. For free hard-core walls, a bosonic current convention gives

\[
x_{\rm charge}=
\frac{\delta W^2}{4K_W}
+K_W\left(D+\frac{Q_w}{2\pi}\right)^2+n_R+n_L.
\tag{24}
\]

Equation (4) is the corresponding exact finite-size quantization.

Near a b,b+1 crossing one can see the same gluing directly. Moving the reference marker changes its position x, so the boundary hop of the η chain acquires a phase \(e^{\pm iP}\). Its Gaussian spectrum is therefore

\[
x_\eta=
\frac{(L-aM)^2}{4K_\ell}
+K_\ell\left(d-\frac{P}{2\pi}\right)^2
+n_R+n_L,
\tag{25}
\]

in a consistent orientation convention. It must be supplemented by (8). The leading current contribution to its cyclic momentum is

\[
q_\ell=\frac{2\pi N_{\rm long}}{M}\,d+O(M^{-1}).
\]

Writing \(aP=Q+2\pi\nu d+2\pi r\) then gives

\[
d-\frac{P}{2\pi}
=\frac{bd-r-Q/(2\pi)}{a}.
\tag{26}
\]

Equations (19) and (26) recover a physical cell-current zero mode twisted by Q. They also show why the physical momentum itself appears as a boundary flux in the projected length chain.

The ground-state energy is obtained by minimizing (21) over the allowed nearby integers M, current sectors, and payload sectors. Therefore its order-1/L coefficient contains both

\[
-\frac{\pi}{6L}(v_\ell+c_pv_p)
\tag{27}
\]

and nonnegative zero-mode contributions. Rounding \(\rho_*L\), parity, and squeezed momentum can all make that coefficient depend on L.

In a stable fluid, the static-marker degeneracy is generically replaced by dispersing states. Ground-state degeneracies are then finite symmetry multiplets or finite-size level crossings; descendant degeneracies follow from the allowed CFT characters. Their detailed values require the payload operator content and gluing rules, not merely c_p.

There are two qualifications to “a finite union of conformal towers on every ring.”

First, every sufficiently large allowed ring can have a CFT description with appropriate zero modes and twists, but it need not approach one L-independent vacuum sector. At irrational density, \(M-\rho_*L\) varies quasiperiodically with L.

Second, a generic compact boson already has infinitely many primary sectors under the Virasoro algebra. A finite number of residue or twist classes is possible under additional assumptions, especially rationality of the payload CFT, but a literal finite union of Virasoro towers does not follow from the construction.

For SU(N)\(_1\), a chain whose length is not divisible by N normally still exists. It occupies a different N-ality sector, often with a nonsinglet ground multiplet. The chiral primary weights are

\[
h_r=\frac{r(N-r)}{2N},\qquad r=0,\ldots,N-1.
\tag{28}
\]

The permitted left-right pairing and momentum follow from the lattice gluing. They add an order-1/L conformal energy and a finite multiplet degeneracy. One should not automatically insert the diagonal value \(2h_r\); odd spin-1/2 chains illustrate why chiral sectors matter.

If the construction literally permits payloads only for \(M\in N\mathbb Z\), restrict the minimization to those M. A fluid can accommodate arbitrary large L by changing its lengths, paying a number-mode energy of order \(N^2/L\). A prescribed incompatible M remains excluded.

4. Central-charge estimators for the XX toy

These predictions apply at fixed J > 0 in the stable two-mode fluid. They are not assertions about every density or arbitrarily strong J.

For half-ring von Neumann entropy,

\[
S(L/2)=\frac{2}{3}\log\frac{L}{\pi a_{\rm UV}}+s_0+\delta S(L).
\tag{29}
\]

Different velocities do not change the coefficient. Near the weakly coupled XX regime, smooth lattice and irrelevant-coupling corrections are commonly of order \(L^{-2}\), with possible logarithms. There is also a useful geometric source of a logarithm. In the factorized approximation, the spin interval contains a fluctuating number \(N_A\) of walls. At half ring,

\[
S_s(N_A)=\frac13\log\left[
\frac{W}{\pi}\sin\frac{\pi N_A}{W}\right]+s_s.
\]

Expanding about \(N_A=W/2\) gives

\[
\delta S_s=-\frac{\pi^2}{6W^2}\operatorname{Var}N_A+\cdots
=O\!\left(\frac{\log L}{L^2}\right).
\tag{30}
\]

This estimate uses factorization, rather than proving the full correction coefficient. Commensurate operators near marginality can produce slower corrections. Rényi entropies can also have oscillatory branch-point corrections; their leading amplitudes need not survive in the von Neumann limit.

For finite-entanglement scaling, if both sectors acquire correlation lengths proportional to one common ξ,

\[
S=\frac{2}{6}\log\xi+s_{\rm FE}+o(1)
=\frac13\log\xi+s_{\rm FE}+o(1).
\tag{31}
\]

With two velocities and a small spin scale, the more informative expression is

\[
S\simeq\frac16\log\xi_\ell+\frac16\log\xi_s+\text{constant}.
\tag{32}
\]

A fit against only the largest transfer-matrix correlation length gives c = 2 only once the two lengths scale proportionally. Finite bond dimension can leave a substantial crossover, or favor a symmetry-broken approximation. There is no universal leading finite-entanglement correction exponent specified by this Hamiltonian alone; irrelevant operators and the MPS truncation both matter.

For the energy,

\[
E_0(L)=e_\infty L
-\frac{\pi}{6L}(v_\ell+v_s)
+\frac{2\pi}{L}\left(v_\ell x_{\ell,0}+v_s x_{s,0}\right)
+\cdots .
\tag{33}
\]

This tests c = 2 only after measuring both velocities and accounting for the ground-sector zero modes. Using one fitted velocity can give a misleading “central charge.” Dimension-three corrections can give \(L^{-2}\) contributions in nontrivial sectors; when their vacuum expectation vanishes, the first smooth vacuum correction is often \(L^{-3}\). Marginal couplings can add logarithms.

At J = 0, the spin velocity vanishes. There is no propagating spin CFT and no universal c = 2 ground-state prescription. Spins retain exponentially many zero-energy choices, subject on a ring to the cyclic-twist constraint. A polarized squeezed-spin state has only the charge contribution, c = 1. Other choices can have additional, even extensive, entanglement at no bulk energy cost.

Taking J → 0 from positive J can select the critical XX spin state, but setting J = 0 first leaves that selection undetermined. The limits are different.

[gpt_solve: backend=codex model=gpt-6-astra effort=high transcript=20261027T200718-63360]
