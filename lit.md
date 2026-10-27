# c-d conjecture, literature survey (2026-10-21)

Tags. [read] = source opened. [recalled] = memory or search snippet only. [own] = my arithmetic or remark, not found in the literature, should be checked.

## 0. The statement

[read, 2403.17242v2, Eqs. (2)-(3)] $H=\sum_i h_{i,i+1}$ with one Hermitian $d^2\times d^2$ matrix $h$ and periodic BC, so one-site translation invariance is assumed. If the low-energy spectrum is "a collection of unitary CFTs", then $c_T=\sum_a c_a\le d-1$. Three-site terms are handled by blocking pairs, giving $c\le d^2-1$ (Eq. 5). The evidence is a table of known models. Only SU(N) Uimin-Lai-Sutherland (ULS, $SU(N)_1$) saturates. Listed open directions are OBC, non-unitary $c_{\rm eff}$, anyonic chains, and a 2D chiral analogue.

Their claim that Motzkin chains have "volume law" entanglement is inaccurate. The cited papers give $\tfrac12\log n$ and $\sqrt n$, and volume law needs the $q$-deformed chain (1606.07795) [read, abstracts].

## 1. Citing papers

INSPIRE (`refersto:recid:2771764`) lists 7. Semantic Scholar and OpenAlex add nothing new.

1. **Kumar, Pujari, D'Emidio, 2506.16424 (PRL 136, 2026)** [read, Suppl. "Non-unitary Violation of the c-d conjecture", Fig. 9]. The model is an SU(N) chain with fundamental and conjugate reps on alternating sublattices, continued to real $N$ as a loop model. It has $c=1-6(g-1)^2/g$ with $N=-2\cos\pi g$. For $N\in(1,2)$, $c>N-1$. They read this as support at integer $N$ and as a violation due to non-unitarity. The violation is formal, since non-integer $N$ is not a Hilbert-space dimension. See 2(e) for a unitary version.
2. **Li, Lin, McGreevy, 2509.04596 (PRL 137, 2026)** [read, pp. 4-6]. A vector-fixed-point search on 4 qubits and 4 qutrits with reconstructed NN Hamiltonians finds only Ising, XX and Heisenberg on qubits. Qutrit states with $c_\Delta>2$ all have entanglement-inferred $c<2$, "consistent with the c-d conjecture". They pitch the method as "a rich source of potential counterexamples".
3. **Prembabu, Verresen, 2509.20431** [read, Sec. II and Outlook]. Their model is the qutrit XX chain ($d=3$, NN, $\mathbb Z_3\times\mathbb Z_3$ projective on each site). It has $c=2$ in the translation-invariant (TI) case, and trimerizing gives $c=8/5$ lines that meet at a $c=2$ point. They note "Our model saturates both conjectural bounds", meaning the lower bound $c\ge d-1$ from LSM (Alavirad-Barkeshli) and the upper c-d bound.
4. **Sang, Hsieh, Zou, 2406.09555 (PRL 133, 2024)** [read, Discussion]. They assume the conjecture, saying CFTs with all $\Delta>1/2$ "may only be realized beyond spin-1/2 chains and 2-local interactions".
5. to 7. **Theses of T.-C. Lin, X. Li and S. Prembabu** [recalled, not opened]. They presumably reproduce items 2 and 3.

No citing paper proves a case or gives an in-scope counterexample.

## 2. Counterexample hunt

**(a) d = 2.**
- O'Brien-Fendley 1712.06662 [read, Eqs. 2-4, Figs. 1, 4]. They add the three-site TI term $H_3=\sum(\sigma^x_j\sigma^z_{j+1}\sigma^z_{j+2}+\sigma^z_j\sigma^z_{j+1}\sigma^x_{j+2})$. The phases are Ising, a TCI point ($7/10$), a $z=3/2$ point, and an incommensurate phase. No $c>1$ is reported.
- Rahmani-Zhu-Franz-Affleck, PRL 1504.05192 [read, abstract], find Ising and TCI. The companion PRB 1505.03966 [read, Eq. 1.7, Fig. 1] has spin form $H=t\sum\sigma^z-t\sum\sigma^x\sigma^x-g\sum\sigma^z\sigma^z-g\sum\sigma^x_j\sigma^x_{j+2}$ ($d=2$, TI, Hermitian, three-site span). Its Ising plus Luttinger phase has $c=3/2$ for $-2.86<g/t<-0.285$. This is the closest non-trivial miss, $c>d-1$ but not NN.
- [own] Trivially, $\sum_j \vec S_j\cdot\vec S_{j+2}$ is two decoupled Heisenberg chains, so $c=2$ at $d=2$. With range $r$ one gets $c=r(d-1)$. The NN restriction is essential, and at fixed $d$, $c$ grows at least linearly with range.
- Verresen-Jones-Pollmann 1709.03508 [read, Eqs. 1-3, p. 2]. Here $c$ is half the number of unit-circle zeros of $f(z)=\sum_\alpha t_\alpha z^\alpha$. [own] With $|\alpha|\le n$ (spin terms spanning $n+1$ sites), $c\le n$, saturated by $H_n+H_{-n}$. The "c = n/2" in the brief is the one-sided case $\alpha\ge 0$.
- [own] **Free-fermion case.** Take a quadratic Majorana chain with $M$ Majoranas per site and range $r$. Then $\det H(k)$ is a trigonometric polynomial of degree $\le Mr$, so it has at most $2Mr$ zeros. Also $c=\#\text{zeros}/4$, which checks on Ising, XX and VJP. So $c\le Mr/2$, and for NN with $d=2^{M/2}$, $c\le\log_2 d\le d-1$. The conjecture therefore holds for all Jordan-Wigner-free NN chains. I did not find this in the literature.

**(b) d = 3.** The largest known value is $c=2$. It occurs at the $SU(3)_1$ ULS point and in the critical region $\theta\in(\pi/4,\pi/2)$ [read, LS Suppl.]. The qutrit XX/XXYY chain also has $c=2$ for all $J\ge0$ (Alavirad-Barkeshli 1910.00589 [read, abstract]; Prembabu-Verresen [read]). I found no $c>2$.

**(c) d = 4.** $SU(4)_1$ ($c=3$) saturates, and spin-3/2 Takhtajan-Babujian gives $SU(2)_3$, $c=9/5$ [read, LS table]. $Sp(4)_2\cong SO(5)_2$ has $c=2\cdot10/5=4$ [own], and I found no $d=4$ realization. $Spin(5)_2\subset SU(5)_1$ appears at the $d=5$ ULS point (2508.13114 [read, abstract]), with $c=4=d-1$. The $d=5$ SO(5) bilinear-biquadratic chain gives $SO(5)_1$, $c=5/2$ (1011.1265 [read, abstract]). Spin-3/2 cold-atom SO(5) models are fermion chains with $d=16$ [recalled].

**(d) Integrable chains.** SU(N) fundamental gives $SU(N)_1$ and saturates. SU(2) spin $s$ gives $SU(2)_{2s}$, $c=3s/(s+1)$ [read, LS]. The SO(N) vector Reshetikhin point gives $SO(N)_1$, $c=N/2$ [read for $N=5$, recalled in general]. SU(N) rank-$k$ symmetric gives $SU(N)_k$ with $d=\binom{N+k-1}{k}$ [recalled]. All except $SU(N)_1$ sit strictly below the bound, and three-site realizations of $SU(2)_{2S}$ still use $d=2S+1$ (1301.5719 [read, abstract]).

**(e) Anyonic chains.** For the golden chain (cond-mat/0612341 [read, Eq. 1 and p. 3]), the Hamiltonian is "three-body" in the fusion-path basis. For su(2)$_k$ anyons, AFM coupling gives $c=1-6/((k+1)(k+2))$ and FM coupling gives $\mathbb Z_k$ parafermions with $c=2(k-1)/(k+2)$.

[own] Read $d$ as the quantum dimension $d_q=2\cos\frac{\pi}{k+2}$. Then AFM violates $c\le d_q-1$ for every $k\ge2$ ($k=2$: 0.5 vs 0.414; $k=3$: 0.7 vs 0.618; large $k$: margin about $(\pi^2-6)/(k+2)^2$). FM violates it more strongly ($c\to2$ while $d_q\to2$). These chains are unitary, Hermitian and TI. As a tensor-product chain they are three-site in the heights, or NN on pair labels with an NN constraint ($d'=2k$, so $d'=3$ for Fibonacci), and they satisfy the bound easily. The NN "Fibonacci" qubit term $\sum|11\rangle\langle11|$ of Li-Lin-McGreevy (Eq. 3.2) only imposes the constraint [read].

[own] This is the Temperley-Lieb algebra of Kumar et al. at $N=d_q$, where $g=(k+1)/(k+2)$ reproduces the minimal-model $c$. So their "non-unitary" violation is really a non-tensor-product effect, and $d$ must be a tensor-factor dimension. An argument that sees only the Hilbert-space growth rate would also cover these constrained NN chains.

**(f) Non-unitary and non-Hermitian chains.** Besides Kumar et al., the non-Hermitian SSH model at its exceptional point has $c=-2$ (1909.01346 [read, abstract]) and Yang-Lee has $c_{\rm eff}=2/5$ [recalled]. Non-unitary interfaces in free-fermion chains have an interface $c_{\rm eff}$ that "grows without bound" (2606.32035 [read, abstract]), but that is a defect quantity, not a bulk $c$. I found no bulk chain with integer $d$ and $c_{\rm eff}>d-1$.

**(g) Non-translation-invariant chains.**
- The random singlet phase has $c_{\rm eff}=c\ln2$ [recalled, Refael-Moore].
- Aperiodic NN XX/XXZ chains with $d=2$ under the tripling sequence give $c_{\rm eff}=\ln4/\ln3\approx1.26>1$ (Juhász-Zimborás cond-mat/0703527, Eqs. 19-21 [read]; also Iglói-Juhász-Zimborás cond-mat/0701527 [read, abstract]). These are aperiodic singlet phases, not CFTs. Any entanglement-based version of the conjecture therefore needs both TI and conformality.
- Periodic modulation with period $p$ reduces to blocking, $c\le d^p-1$.

**(h) NN chains that are not conformal.**
- The Bravyi et al. Motzkin chain ($d=3$, NN, bulk TI) has $S_{n/2}=\tfrac12\log n$ [read, abstract]. [own] Via the open-chain formula $(c/6)\log n$, this "$c_{\rm eff}$" is $3>d-1$, but $z>1$ [recalled]. Colored Motzkin has $\sqrt n$ entanglement and gap $\le n^{-2}$ (1408.1657 [read]).
- Irani 0901.1107 builds a TI NN chain with $d=21$, entanglement $\Omega(\min\{m,n\})$ and gap $1/{\rm poly}(n)$ [read, p. 2]. Gottesman-Irani 0905.2419 show TI NN ground energies are QMA$_{\rm EXP}$-complete [read, abstract].

## 3. Tools for a proof

**(a) Kapustin-Spodyneiko 1904.05491 (PRL 123, 060601)** [read, abstract and p. 2]. Equilibrium energy currents vanish for 1D lattice systems with finite-range, uniformly bounded terms and a clustering Gibbs state, so there is no flow to a CFT with $c_L\ne c_R$. TI is not needed, and finite-dimensional sites suffice (not invoked explicitly). So $c$ in the conjecture is unambiguous.

**(b) LSM lower bounds.**
- Yao-Hsieh-Oshikawa 1805.06885 [read, abstract]. A nontrivial index (Young boxes per unit cell mod N) restricts the critical theory to a class of $SU(N)_k$.
- [own] The unitary Sugawara bound $c\ge c(SU(N)_k)\ge N-1$ then gives $c\ge d-1$ for the fundamental. Together with c-d this forces $SU(N)_1$.
- Alavirad-Barkeshli 1910.00589 [read, abstract] suggest a $c$ lower bound growing with $n$ for $\mathbb Z_n\times\mathbb Z_n$ LSM.
- Lin-Shao 2207.05092 [read, abstract] give bootstrap bounds on (charged) operators of $\mathbb Z_N\times\mathbb Z_N$-LSM CFTs. These are not bounds on $c$.

**(c) Upper bounds on $c$ or on the log coefficient.** I found none in terms of $d$, the Lieb-Robinson velocity, or $\|h\|$. Area laws (Hastings, AKLV) need a gap [recalled]. Brandão-Cramer 1409.5946 [read, abstract] show that specific-heat decay implies an area law, or a subvolume law for polynomial decay.

[own] A heuristic combines $s(T)=\pi cT/(3v)\le\log d$ with $v\le v_{\rm LR}=O(\|h\|)$, giving $c\lesssim 3v_{\rm LR}\log d/(\pi T_*)$. It needs control of the crossover temperature $T_*$, and by 2(e) it cannot yield $d-1$ without tensor-product input.

**(d) Getting $c$ from lattice data.**
- Koo-Saleur hep-th/9312156 and Milsted-Vidal 1706.01436 [read, abstracts] build lattice Virasoro generators and extract conformal data.
- Zini-Wang 1706.08497 [read, abstract] give a rigorous scaling limit, proved for Ising only. Osborne-Stottmeister 2107.13834 [read, abstract] prove Koo-Saleur convergence for free fermions and WZW models.
- $E_0=e_0L-\pi cv/(6L)$ (Blöte-Cardy-Nightingale, Affleck 1986) [recalled].
- The modular commutator sees only $c_L-c_R$ (2206.00027 [read, abstract]), so by 3(a) it vanishes here. The Casini-Huerta entropic c-theorem needs Lorentz invariance [recalled].

## Status

The conjecture is open. It is believed true in its stated form, which assumes a TI, strictly NN, Hermitian chain on a tensor product of $d$-dimensional sites, flowing to unitary CFT(s), with $c_T\le d-1$. There is no proof in the literature for any class. [own] It does hold for all free-fermion NN chains, where $c\le\log_2 d$.

No counterexample in scope is known. Every reported "violation" drops one hypothesis:
- non-integer $N$ (Kumar et al.);
- $d$ read as a quantum dimension (anyonic chains, [own]);
- broken TI (aperiodic $c_{\rm eff}\approx1.26$ at $d=2$);
- longer range ($c=3/2$ or $2$ at $d=2$);
- no conformality (Motzkin).

Saturation occurs for $SU(N)_1$ ULS and the $\mathbb Z_3\times\mathbb Z_3$ qutrit XX chain. In the LSM classes (SU(N) fundamental, or $\mathbb Z_d\times\mathbb Z_d$ projective), the expected lower bound $c\ge d-1$ meets the conjectured upper bound, which would pin $c=d-1$.
