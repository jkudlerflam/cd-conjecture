# The first number-changing term in the cusp reduction

AI assistance: GPT-6 Codex, 2026-11-09. Independent derivation and implementation; no C173 or C175 code is imported.

## 1. Statement and conventions

[proved, algebra] Let `b_j` be spin lowering/raising hard-core operators, `n_j=b†_j b_j`, and let `P` project onto configurations with no adjacent occupied sites. Consider a ring with L >= 7, or an open chain, with

```
H = V D + Delta T + (Delta^2/V) W,
D = sum_j n_j n_(j+1),
T = sum_j (b†_j b†_(j+1) + h.c.),
W = sum_j (w0 b†_(j+1) b_j + h.c.) - mu0 sum_j n_j.
```

Thus the physical hopping is `w=w0 Delta^2/V`. All coefficients below use the spin occupation basis, avoiding Jordan–Wigner boundary choices. Write `R=Q D^-1 Q`, `Q=1-P`. A canonical, Hermitian Schrieffer–Wolff reduction, with no additional block-diagonal rotation, gives

```
H_eff = (Delta^2/V) K2 + (Delta^3/V^2) K3 + higher orders,
K2 = P W P - P T R T P,
K3 = P T R T R T P - P T R W P - P W R T P.
```

The usual normalization anticommutator first contributes at the next order because `P T P=0`. The third-order coefficient is exactly

```
K3 = (2/3) P sum_j (b†_j b†_(j+3) + h.c.) P
     - 2 Re(w0) P sum_j (b†_j b†_(j+2) + h.c.) P.           (1)
```

For open chains only terms with both endpoints inside the chain are summed. Projectors apply to the entire expression. In particular, the distance-three creation process requires its two interior sites and its exterior neighbours to be empty whenever those neighbours exist.

## 2. Local path proof

[proved, algebra] A pairing operation changes particle number by two. Three such operations, beginning and ending in `P`, can only create or annihilate a separated pair. A surviving creation process has support on four consecutive empty sites and is

```
0000 -> 1100 -> 1111 -> 1001,
0000 -> 0011 -> 1111 -> 1001.
```

The intermediate adjacency counts are 1 and 3. The two paths contribute `1/(1*3)+1/(1*3)=2/3`. The reversed paths give annihilation. Three creations would leave an adjacent pair, and a disconnected creation would do the same. Repeating a bond cannot create an allowed final pair. These observations exhaust the three-pair paths. External spectators do not change the denominators: any spectator adjacent to a final endpoint is forbidden by the final `P` projection.

For one pairing operation and one hopping operation, the allowed pair creation is at distance two. It has two paths,

```
000 -> 110 -> 101,
000 -> 011 -> 101.
```

The denominators are both 1 and the hops carry `w0` and `conj(w0)`. The minus sign in `-P T R W P-P W R T P` gives `-2 Re(w0)`. The diagonal chemical potential has no `P-Q` matrix elements and contributes nothing at this order. This proves (1), including cancellation of the distance-two term for imaginary hopping.

## 3. Consequence for the actual candidate

[proved, algebra] At C173's centre, `Delta=1`, `w0=i/2`, and `V=40`, the surviving coefficient is

```
gamma3 = 2 Delta^3/(3 V^2) = 1/2400,
|t2| = Delta^2/(2V) = 1/80,
gamma3/|t2| = 4 Delta/(3V) = 1/30.
```

Consequently, the full local effective Hamiltonian does not differ from its second-order, number-conserving truncation only at relative order `(Delta/V)^2`. It already contains a nonzero number-changing term at relative order `Delta/V`.

This corrects the scope of the robustness argument in [C175 PROOF section 2.2](../../codex/work_C175/PROOF.md). Its estimate may describe the first corrections to the *number-conserving coefficients*, or the effects of the number-changing perturbation after a separate off-resonant elimination. It is not an estimate on the full local remainder. A finite-size spectral error proportional to `V^-3` does not establish such a remainder estimate.

[derived, with stated low-energy assumptions] This correction does not negate C175's leading two-body binding calculation. At the symmetric point the valley momenta satisfy `k_a+k_b=pi`. The uniform distance-three pairing has total momentum zero. Its direct low-energy pair creations therefore carry

- `2 k_a` for two particles in valley a,
- `2 k_b` for two in valley b,
- `k_a+k_b=pi` for one in each valley.

For generic `0<|w|<4` none of these is zero modulo `2pi`. A direct resonant atomic pairing term is absent. For an isolated eigenstate of the number-conserving `K2`, the first-order diagonal matrix element of `K3` also vanishes because it changes N by two. Under controlled nonzero denominators the resulting energy shift can thus start at relative order `(Delta/V)^2` even though the operator perturbation is relative order `Delta/V`.

Neither momentum counting nor the finite-ring checks below supplies a bound on all thermodynamic denominators. A microscopic persistence claim requires that extra step. In particular, two inter-valley pairs have momentum `2pi`; a `Delta N=4` vertex is allowed and its coefficient and scaling dimension must be controlled. C180 already lists that vertex. If the low-energy objects become bound pairs, the stability analysis must use the pair liquid's parameters; atomic-gas values of `K_+` cannot simply be carried over. No nonzero quartet coefficient or microscopic central charge is established here.

## 4. Independent computational checks

[numerical] `check_reduction.py` builds `D,T,W` directly in the full spin basis, projects the exact perturbation matrices, and independently builds (1) by separated spin-pair flips. It makes 45 comparisons:

- open chains L=4 through 12;
- periodic rings L=7 through 12;
- `w0=i/2`, `w0=0.15+0.35i`, and `w0=1/2`;
- `mu0=1.9`.

Every matrix entry agrees to better than `1e-12`. Every third-order entry whose particle-number change is not plus or minus two vanishes to the same tolerance. Exact maximum errors and dimensions are in `results.json`.

[numerical] The script also compares the lowest 20 full-spin eigenvalues on an L=12 periodic ring to `K2/V` and `K2/V+K3/V^2`, at the candidate chemical potential `mu_eff=-2.2`:

| V | maximum error, K2 only | maximum error, K2 plus K3 |
|---:|---:|---:|
| 10 | 7.1501e-3 | 5.4208e-3 |
| 20 | 9.7752e-4 | 5.6456e-4 |
| 40 | 1.2600e-4 | 8.6147e-5 |
| 80 | 1.5888e-5 | 1.0778e-5 |

The two approximations both retain order `V^-3` errors; number-conserving fourth-order terms have not been included. Adding only `K3` is not a fourth-order approximation.

A limitation matters here. On this short ring the second-order ground state is the vacuum because the allowed momenta miss both band minima. The full ground state's particle number tends to zero as V grows (0.0254 at V=40 and 0.00639 at V=80). These data check coefficients and spectra at the candidate *couplings*, not the proposed density `n=0.17`. No thermodynamic or finite-density claim follows from them.

The C173 `blockade.py` example uses `mu=-3/V+0.5/(2V)`, which translates to `mu_eff=-11.5`; the proposed candidate uses `mu=(mu_eff+6)/(2V)=1.9/V`. The older example is therefore a vacuum-side spectral check, not a density-matched check. Our script additionally compares that chemical potential on periodic rings with imaginary hopping, explicitly as a control rather than a reproduction of C173's open-chain run and different hopping phase.
