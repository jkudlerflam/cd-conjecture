# Exponentially small thermal weight of frozen defects

## Claim and status

**[proved, using the exact local-sector hypotheses below]** For the frozen
23-letter route-A chain at the selector parameters already established in
[`route_a/README.md`](../../counterexample_campaign/route_a/README.md), the
infinite-volume free energy per physical site differs from that of the exact,
interacting, homogeneous selected-cell chain, divided by the cell length,
by at most `C T exp(-delta/T)` at sufficiently small positive temperature.
Here `C, delta > 0` are fixed independently of the ring length and temperature.
Every wrong cell length and every illegal-arrow sector is included.

Consequently, **if** the selected-cell chain has a quadratic low-temperature
free-energy term, the complete frozen microscopic chain has exactly the same
coefficient after the physical cell-length conversion. This does not establish
that quadratic term or its decomposition into Virasoro central charges.

The theorem uses trace inequalities and finite-dimensional state counting.
It does not use a spectral gap inside the selected sector, an effective
Hamiltonian truncation, a convergent perturbation series, or a CFT assumption.

## 1. Exact input from the microscopic bond

Fix one bond and parameters, with local dimension `d=23`, selected cell
length `b>=2`, and write `B=sqrt(3)|lambda|`. The following are the actual
properties proved for that bond in the linked global-sector argument:

1. The marker-position projectors and forbidden-arrow-bond projectors
   commute with the Hamiltonian. Thus the trace splits into orthogonal
   sectors of their classical positions.
2. A complete legal cell of length `ell` has ground energy

       E_ell = mu_sel - 2 cos(pi/(ell+1)).

   Its Hilbert-space dimension is at most `d^ell`. Every cell Hamiltonian
   includes its onsite terms. The conditions ensuring this lower bound at
   every length are `mu_sel, epsilon, D_x, D_M >= 0` and
   `Delta-2-|t_d|-max(|g_+|,|g_-|)>0`.
3. Adjacent complete cells couple only by the boundary port with norm at
   most `B`. The marker-controlled junction term is internal to the cell on
   its right because it acts identically on every marker state.
4. A legal end fragment has at most one token and ground energy at least
   `-2`, irrespective of its length. Its dimension is at most `d^ell`.
5. Each illegal bond costs `K`. Cutting all `r` illegal bonds produces
   `r` legal intervals. Relative to their complete-cell decomposition there
   are at most `2r` end fragments. If there are `k` wrong-length complete
   cells, there are at most `k+r` runs of selected cells and at most
   `2k+4r` port bonds touching a wrong cell or end fragment.
6. Put

       e_b = E_b/b = -a < 0,
       gamma = inf_(ell != b) (E_ell - ell e_b) > 0.

   Assume `gamma>3B` and `K>4+5B`. The existing stronger conditions
   `gamma>5B` and `K>4+7B` suffice.

The dimension bounds are intentionally crude. They sum all data and excited
states rather than retaining only the code. In particular, the proof does
not discard the growing internal degeneracy of arbitrarily long cells.

## 2. Exact finite-temperature boundary estimate

Let `Z_m^o(beta)` be the partition function of `m` selected cells with open
boundaries, retaining their entire finite Hilbert spaces and the exact
coupling. Write

    f_cell(beta) = -lim_(m -> infinity) log Z_m^o(beta)/(beta m),
    e(beta) = f_cell(beta)/b.

If self-adjoint matrices `A,C` on the same finite space obey
`||A-C||<=D`, the min-max eigenvalue inequality gives

    exp(-beta D) Tr exp(-beta C)
       <= Tr exp(-beta A)
       <= exp(beta D) Tr exp(-beta C).

Concatenating `q` open chains of `m` cells adds `q-1` bonds of norm `B`.
The same almost-additivity argument that establishes the thermodynamic
limit therefore gives, for every `m>=1`,

    |log Z_m^o(beta) + beta m f_cell(beta)| <= beta B.       (1)

For rings with at least three cells, closing the chain gives

    |log Z_m^p(beta) + beta m f_cell(beta)| <= 2 beta B.     (2)

For `m=0`, set `Z_0^o=1`, with no boundary factor. Exceptional one-cell
and two-cell periodic bond conventions are irrelevant to the limit and are
not used in (2). Open runs of one or two cells still satisfy (1). A product of empty-data
vacuum code states has energy `m E_b` and zero port expectation. The Gibbs
variational principle then gives

    e(beta) <= e_b = -a.                                  (3)

All these statements hold for arbitrary positive beta. No logical gap is
assumed.

## 3. Bound on each defect sector

First consider a legal ring with `k>=1` wrong cells. Remove all port bonds
touching a wrong cell; there are at most `2k`. Bound each resulting selected
run using (1), with at most `k` runs. Bound each wrong cell trace by its full
dimension times the exponential of its ground-energy lower bound. After
multiplication by `exp(beta e(beta)L)`, its normalized trace is at most

    product_wrong_cells w_ell(beta),

where

    w_ell(beta) = d^ell exp[-beta(E_ell - ell e_b - 3B)],
    ell != b.                                            (4)

Replacing `e(beta)` by `e_b` only enlarges this bound by (3).

For a ring with `r>=1` illegal bonds and `k` wrong complete cells, remove
port bonds touching the wrong cells or end fragments. The norm cost is at
most `(2k+4r)B`; selected-run factors contribute at most `(k+r)B`; the end
fragment energies contribute at most `4r`. Thus its normalized trace is at
most

    exp(-beta r delta_I) product_wrong_cells w_ell(beta)
       times product_end_fragments [d^ell exp(-beta a ell)],

    delta_I = K-4-5B > 0.                                 (5)

A zero-length fragment contributes one. If a legal interval has no marker,
regard its entire length as its left fragment and set its right fragment
to zero. This avoids requiring a marker in every interval.

## 4. Summing all lengths and all positions

This section gives an intentionally redundant upper count. Choose the `r`
illegal bonds and the first sites of the `k` wrong complete cells. There
are at most `binom(L,r) binom(L,k)` choices. For each wrong cell, sum its
length in (4). For each illegal interval, sum the two end-fragment lengths
independently from zero to infinity. Restrictions from total length and
non-overlap are dropped, which only raises this positive sum.

A sector description is recovered from these data together with the
selected-run offsets: wrong-cell starts and lengths give their marker
positions, each interval prefix ends at its first marker, and all remaining
complete cells have length `b`. The convention for marker-free intervals
was stated after (5). Cases in which an entire end fragment may or may not
contain a marker are combined as orthogonal subspaces of the same fragment
Hilbert space before applying its `d^ell` bound. Thus distinct such choices
are not silently discarded.

The remaining sites are selected-cell runs. Given their endpoints and a
cell length `b`, there are at most `b` possible marker offsets in each
nonempty run. There are at most `k+r` runs. Thus `b^(k+r)` is a safe upper
bound for their residual classical multiplicity. The internal states of
wrong cells and fragments have already been counted by their `d^ell`
factors. Different endpoint arrow choices are contained in those full
fragment Hilbert spaces and incur no additional state-counting assumption.

Define

    W(beta) = b sum_(ell>=1, ell!=b) w_ell(beta),
    q(beta) = d exp(-beta a),
    I(beta) = b exp(-beta delta_I)/(1-q(beta))^2.            (6)

For `q(beta)<1`, summing (4)-(5), including all redundant choices above,
is bounded by

    [(1+W(beta))(1+I(beta))]^L.                            (7)

One may include the no-defect choice in (7) even though it is accounted
separately. A completely selected ring contributes at most
`b exp(2 beta B)` to the normalized trace. A marker-free legal ring has
nonnegative energy and contributes at most
`d^L exp(beta e(beta)L) <= q(beta)^L`. These observations give the explicit
finite-volume bound

    exp(beta e(beta)L) Z_L(beta)
       <= b exp(2 beta B) + q(beta)^L
          + [(1+W(beta))(1+I(beta))]^L.                    (8)

The bound may overcount many nonexistent sector descriptions; all weights
are positive, so this is harmless. It applies equally when `b` does not
divide `L`.

## 5. Uniform exponential tail for the wrong-length sum

Set

    delta_W = gamma-3B > 0,
    C0 = 2+3B,
    kappa = a delta_W/(a+C0+delta_W) > 0.                  (9)

Since `mu_sel>=0`, every wrong-cell cost in (4) obeys

    E_ell - ell e_b -3B >= max(delta_W, a ell-C0)
                         >= kappa(ell+1).                (10)

For the second inequality, take the convex combination of the two bounds
with coefficient `kappa/a` on `a ell-C0`. Definition (9) makes its constant
term precisely `kappa`. This proves the estimate for every length, rather
than only below a cutoff.

If `d exp(-beta kappa)<=1/2`, then

    W(beta) <= b d exp(-2 beta kappa)
                   /(1-d exp(-beta kappa))
            <= 2 b d exp(-2 beta kappa).                 (11)

If also `q(beta)<=1/2`,

    I(beta) <= 4 b exp(-beta delta_I).                    (12)

All constants depend on the fixed microscopic parameters only.

## 6. Thermodynamic conclusion

Let `f_full(beta)` be the ordinary thermodynamic free energy per physical
site of the full `d`-state nearest-neighbour ring. Its existence follows
from the same bounded-boundary almost-additivity argument as (1). Taking
the logarithm of (8), dividing by `-beta L`, and sending `L` to infinity
proves

    f_full(beta) >= e(beta)
       - beta^(-1) log[(1+W(beta))(1+I(beta))].             (13)

Conversely, on rings of length `L=mb` the full partition function contains
the selected periodic sector. Equation (2) gives
`f_full(beta)<=e(beta)` in the limit. Existence of the full limit makes this
valid for all ring-length sequences. Alternatively, adding at most `b-1`
sites changes `log Z` by a bounded quantity at each fixed beta: comparing
`H_(mb+r)` to `H_(mb) tensor 1` changes at most `r+2` bond terms, of total
norm at most `(r+2)||h||`, and supplies the finite factor `d^r`.

Combining (11)-(13), with `T=1/beta`, yields

    0 <= f_cell(T)/b - f_full(T)
       <= (2bd+4b) T exp(-delta/T),

    delta = min(2 kappa, delta_I) > 0,                   (14)

for all sufficiently small positive `T`. In particular the ground-energy
density is the same, and the defect correction is smaller than every
positive power of `T`.

## 7. What this proves for the proposed counterexample

Suppose, as an additional input, that the exact selected-cell chain has

    f_cell(T) = u - (pi/6) sum_a(c_a/v_a^cell) T^2 + o(T^2).

Then (14) proves for the full frozen chain

    f_full(T) = u/b - (pi/6) sum_a(c_a/v_a^physical) T^2 + o(T^2),
    v_a^physical = b v_a^cell.

This is a rigorous transfer of the assumed thermal coefficient. It does
not identify the individual `c_a` from one unweighted number, and it does
not prove conformal finite-size spectra on incompatible rings. A family
of independently resolved source responses or a suitable scaling-limit
theorem is still needed for the complete CFT interpretation.

At C33's `n=23` point, `b=24`, the existing rational selector certificate
allows `0<|lambda|<=2^-20`, `K=10` and implies all hypotheses above.
The code in this folder computes conservative rational witnesses for
`a, delta_W, kappa, delta_I` from that certificate without running any
large Hilbert-space diagonalization. The resulting positive delta can be
extremely small; positivity at fixed nonzero parameters is what is needed.

**[research judgment]** For a bulk infrared interpretation of the original
conjecture, this removes the need to prove the separate vacancy-fluid
matching merely to eliminate frozen defects from the thermal coefficient.
The frozen d=23 chain needs `c=23>22`, and its finite band certificate
already exists. Its remaining substantive task is the microscopic
selected-cell infrared theory, including a sufficient way to identify
all 23 central-charge contributions.

Provenance: Codex GPT-6-family subagent, 2026-10-28. Derived independently
from the existing exact sector inequalities. Not human-refereed. The
coordinating agent independently checked the trace inequalities, sector
counts and exponential-tail algebra against the exact local-sector source.
No independent numerical central-charge measurement is claimed.
