# Critical qubit examples: details recovered from the prior chat

## Scope and provenance

This supplements the short list in [`specs/_d2_background.md`](../../specs/_d2_background.md), rather than claiming a classification or a new proof of the conjecture. The source is the assistant's catalog in **“Prove c d conjecture”**, retrieved on 2026-10-27. The retrieved message is truncated at its final sentence. Its original model is not identified in the retrieved record. A GPT-6-based Codex assistant selected the material below, checked the indicated primary sources, and rederived the two elementary mappings. No human collaborator has independently reviewed this import and no numerical simulation was rerun.

The labels **[read; exact]**, **[read; numerical]**, and **[derived]** distinguish published exact results, published numerical evidence, and elementary algebra written out here. They do not certify a rigorous lattice-to-CFT theorem. Central charges are ordinary unitary Virasoro charges, with Ising normalized to $c=1/2$; $c_T$ is their sum. No non-unitary or non-compact example is used.

Write $X,Y,Z$ for Pauli matrices and $\mathbf S=(X,Y,Z)/2$. “Strict NN” means support on at most two adjacent physical sites. A three-site term has range two. All couplings displayed below are real. Unless explicitly qualified, the bulk is uniform and the Hilbert space is the full $(\mathbb C^2)^{\otimes L}$.

## Strict nearest-neighbour examples

The free XY/XX, axial-DM XXZ, and zero-field XYZ cases already have substantive treatments in [`codex/work_C6_d2/c6.tex`](../../codex/work_C6_d2/c6.tex). The additions most useful for the remaining interacting problem are:

| Family | Hamiltonian or parameterization | Critical behavior and evidence |
|---|---|---|
| Antiferromagnetic Ising in mixed fields | $J\sum_j Z_jZ_{j+1}-g\sum_jX_j-h\sum_jZ_j$, $J>0$ | **[read; numerical]** Ising $c=1/2$ critical curve. In this Pauli normalization it connects the transverse Ising point $(g,h)=(J,0)$ toward the classical endpoint $(0,2J)$. The classical endpoint is not an ordinary Ising CFT. The ordered phase breaks translation, so a uniform longitudinal field need not remove the transition. [Ovchinnikov et al.](https://arxiv.org/html/cond-mat/0306468) |
| Interacting Kitaev chain in spin variables | $\sum_j[-(t+\lambda)X_jX_{j+1}/2-(t-\lambda)Y_jY_{j+1}/2+VZ_jZ_{j+1}/4-(V-\mu)Z_j/2]$ | **[read; exact / numerical]** Exact $c=1$ locus below; extended floating $c=1$ regions at strong repulsion have numerically located boundaries. For generic nonzero pairing the microscopic particle symmetry is $\mathbb Z_2$, while the floating phase has emergent $U(1)$. [Chepiga–Mila](https://arxiv.org/html/2206.11754) |
| Dual ANNNI | $\sum_j[J_2\tau_j^x\tau_{j+1}^x-\Gamma\tau_j^z\tau_{j+1}^z-J_1\tau_j^x]$ | **[derived mapping; read, numerical phase]** Strict NN interacting qubits with the $c=1$ floating regime inherited from ANNNI. See mapping below and [Cea et al., Appendix C](https://arxiv.org/html/2402.11022). |

### An explicit exact $c=1$ locus beyond uniform on-site $U(1)$

For the Kitaev parameterization above, take

$$
t>0,\qquad V>2t,\qquad \mu=V,\qquad
\lambda=(V-2t)/2.
$$

**[derived]** The field vanishes and $(J_x,J_y,J_z)=(-V/4,V/4-t,V/4)$. Rotate every other spin by $\pi$ about $z$. On every bulk bond the exchanges become

$$
(J'_x,J'_y,J'_z)=(V/4,t-V/4,V/4).
$$

An axis permutation gives antiferromagnetic XXZ with

$$
J=V/4>0,\qquad \Delta=4t/V-1\in(-1,1).
$$

Thus the standard XXZ continuum result supplies $c=1$. The staggered transformation is literal on even periodic chains; other boundary choices require a twist. **[read]** This is the eight-vertex critical locus identified by [Chepiga–Mila, Eqs. (1)–(3) and following discussion](https://arxiv.org/html/2206.11754). It is an exactly located sublocus, not an exact solution of the surrounding floating-phase boundaries.

### Why the dual ANNNI model belongs in the NN list

The original variables have range two:

$$
H_A=-J_1\sum_jZ_jZ_{j+1}+J_2\sum_jZ_jZ_{j+2}
-\Gamma\sum_jX_j.
$$

**[derived]** In the bulk, set

$$
\tau_j^x=Z_jZ_{j+1},\qquad
\tau_j^z=\prod_{\ell\le j}X_\ell.
$$

Then $Z_jZ_{j+2}=\tau_j^x\tau_{j+1}^x$ and $X_j=\tau_{j-1}^z\tau_j^z$, yielding the NN Hamiltonian in the table. On a ring, parity sectors and boundary twists must be tracked; this is not a sector-free unitary identification of two periodic Hilbert spaces. Nevertheless the displayed dual Hamiltonian independently defines a uniform full-tensor-product qubit chain. The finite set of global sector choices does not introduce a bulk propagating mode.

**[read; numerical]** ANNNI's floating phase has the Luttinger-liquid description $c=1$; its boundaries include Kosterlitz–Thouless and Pokrovsky–Talapov transitions. [Cea et al.](https://arxiv.org/html/2402.11022) extract $c$ from open-chain entanglement scaling, including $N=240$. No particular decimal coupling point from the earlier chat is promoted here to an exact or independently reproduced benchmark.

## Full qubits with longer-range interactions

These examples help identify which microscopic assumption excludes an apparent counterexample. They do not satisfy strict NN support.

| Family | Uniform Hamiltonian | Continuum and status |
|---|---|---|
| Generalized cluster chain | $H_A=\sum_j[Y_jZ_{j+1}\cdots Z_{j+A-2}Y_{j+A-1}+JX_jX_{j+1}+hZ_j]$, $A\ge3$ | **[read; exact]** At $(J,h)=(\pm1,0)$, $SO(A)_1$ with $c=A/2$; at $(0,\pm1)$, $SO(A-1)_1$ with $c=(A-1)/2$. [Lahtinen–Ardonne](https://arxiv.org/html/1504.07044) |
| O'Brien–Fendley | $-\sum_j(X_j+Z_jZ_{j+1})+\lambda\sum_j(X_jZ_{j+1}Z_{j+2}+Z_jZ_{j+1}X_{j+2})$ | **[read; numerical]** An Ising segment reaches tricritical Ising, $c=7/10$, at $\lambda\approx0.428$ in this normalization. The paper uses $2\lambda_IH_I+\lambda_3H_3$, locating $\lambda_3/\lambda_I\approx0.856$. Spectral ratios and entanglement supply independent diagnostics. [O'Brien–Fendley](https://arxiv.org/html/1712.06662) |
| Three-spin Ising with two fields | $-\sum_jZ_jZ_{j+1}Z_{j+2}-h_\ell\sum_jZ_j-h_t\sum_jX_j$ | **[read; numerical]** A critical curve at nonzero longitudinal field has three-state Potts behavior, $c=4/5$, from finite-size scaling. This is a full-qubit realization using three-site interactions. Exact endpoints and a precise curve are not asserted here. [Penson–Debierre–Turban](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.37.7884) |
| Interleaved XXZ | $J\sum_j(S_j^xS_{j+R}^x+S_j^yS_{j+R}^y+\Delta S_j^zS_{j+R}^z)$, $J>0$, $-1<\Delta\le1$ | **[derived from the standard XXZ result]** For $R\mid L$, the residue classes modulo $R$ are $R$ independent chains, hence $c_T=R$. Translation by one site permutes those chains. Only two-body terms are required, but their range is $R$. |

**[read; exact]** In the cluster family, the common fermionic boundary-sector choice matters for the operator spectrum: the resulting spin CFT is $SO(N)_1$, rather than a tensor product with independently chosen Ising sectors. In particular, $\sum_j(Y_jZ_{j+1}Y_{j+2}+X_jX_{j+1})$ realizes $SU(2)_2\simeq SO(3)_1$ with $c=3/2$. [Lahtinen–Ardonne](https://arxiv.org/html/1504.07044)

## Constraints, cells, and nonconformal examples

**[read; exact critical points]** The Fendley–Sengupta–Sachdev hard-boson model is

$$
H=\sum_j[-w(d_j+d_j^\dagger)+Un_j+Vn_jn_{j+2}],
\qquad n_j\in\{0,1\},\quad n_jn_{j+1}=0.
$$

With $w>0$ and $\varphi=(1+\sqrt5)/2$, its integrable line is $w^2=UV+V^2$. It contains

| Critical theory | $V/w$ | $U/w$ |
|---|---:|---:|
| Tricritical Ising, $c=7/10$ | $-\varphi^{5/2}$ | $\varphi^{5/2}-\varphi^{-5/2}$ |
| Three-state Potts, $c=4/5$ | $+\varphi^{5/2}$ | $\varphi^{-5/2}-\varphi^{5/2}$ |

These values follow from Eqs. (4), (6), and (8) of [Fendley–Sengupta–Sachdev](https://arxiv.org/html/cond-mat/0309438). The constraint reduces the Hilbert space below $2^L$. In a full-qubit representation the permitted flip is $P_{j-1}X_jP_{j+1}$, with $P_j=1-n_j$, so it already has three-site support. Calling empty/occupied a qubit does not put these points in the strict NN class.

Other distinctions retained from the catalog:

- **[read]** The self-dual Ashkin–Teller critical line has an orbifold $c=1$ description. Its usual quantum-cell representation uses two Ising species, hence dimension four per cell; interleaving them does not preserve strict NN support on individual qubits. This is a useful warning about counting cells. [Alba–Tagliacozzo–Calabrese](https://arxiv.org/abs/1103.3166)
- **[read]** Antiferromagnetic Haldane–Shastry chains give a qubit $c=1$ example with inverse-square exchange. Their long-range coupling excludes them from the present microscopic class. [Finkel–González-López](https://arxiv.org/abs/2208.04014)
- **[derived]** Ferromagnetic Heisenberg magnons obey $E(k)-E(0)\propto1-\cos k\sim k^2$, so gaplessness at this point is not a relativistic CFT. Saturation and other quadratic band edges must likewise be distinguished from neighboring $c=1$ phases.
- **[read; numerical dynamics]** The Fredkin Hamiltonian has bulk terms $P_j^\uparrow P^s_{j+1,j+2}+P^s_{j,j+1}P_{j+2}^\downarrow$, with singlet projector $P^s$ and boundary terms selecting a Dyck state. It is a full-qubit model with three-site interactions and an exactly known ground state, but its dynamics have $z>1$. A logarithmic entropy coefficient alone is not its Virasoro central charge. [Adhikari–Beach](https://arxiv.org/abs/2011.07110)

## What this adds, and what it does not

The explicit interacting and dual examples prevent a general NN qubit argument from assuming free fermions or a microscopic uniform $U(1)$. The larger-range and constrained examples prevent incorrect counterexample claims based solely on the label “spin one-half.” The NN conformal examples collected here have $c=1/2$ or $1$, but their existence does not exhaust the general exchange-matrix-plus-field family or prove $c_T\le1$.

This is a literature and algebra supplement. There are no new numerical artifacts to reproduce. The exact algebra can be checked directly from the displayed Pauli rotations and duality identities; published numerical claims remain attributed to their original sources.
