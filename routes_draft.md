# Working notes for REPORT.md (Opus, 2026-10-21). Superseded by REPORT.md; kept for the derivations.

## R1. Entropic no-go (checked by hand)
Lattice block entropy S(l) of l consecutive sites in the ground state.
Subadditivity + Araki-Lieb: |S(l+1) - S(l)| <= S(site) <= log d.
Discrete Casini-Huerta function C_CH(l) := 3 l [S(l+1) - S(l)] <= 3 l log d, so C_CH(1) <= 3 log d.
If C_CH were non-increasing in l and tended to c, then c <= 3 log d.
ULS / SU(d)_1: c = d - 1. d - 1 > 3 ln d first at d = 7 (6 > 5.84). So for every d >= 7 the lattice
CH function of the SU(d) chain must increase somewhere. With the wedge-interval MI version
(C = -6 a dI/da, chain rule I(L:A) - I(L':A) = I(x:A|L') <= 2 S(x) <= 2 log d per site, so C(a=1) <= 12 log d)
the contradiction starts at d = 48. With K = 6 at d = 19.
General statement: any C(r) that is (i) a log-derivative of entanglement data, (ii) <= K log d at one lattice
spacing, (iii) monotone, (iv) -> c, gives c <= K log d, false for SU(d)_1 at large d.
Also: SU(d)_1 needs S(l*) >= ((d-1)/3) log(l*/eps) while S(l*) <= l* log d, so the conformal regime can only
begin at l* >~ (d-1)/(3 log d) up to logs. A proof has to know where conformality starts.

## R2. Thermodynamics, Lieb-Robinson, sum rules (checked by hand, standard formulas)
CFT window T <= T*: s(T) = pi c T/(3 v). Lattice: s <= log d. So c T*/v <= (3/pi) log d.
Lieb-Robinson: v <= v_LR = O(||h||). So c T* <= (3/pi) v_LR log d.
Energy current in a CFT NESS: J = (pi c/12)(T_L^2 - T_R^2); lattice ||J_E|| <= 2||h||^2 gives c T*^2 <= (24/pi)||h||^2.
f-sum for an on-site current algebra G_k: (k/2)(2 pi v/N) <= (1/2)<[S_{-q},[H,S_q]]> = O(q^2 N) gives k v <= const.
All are invariant under (c, T*) -> (lambda c, T*/lambda) up to the v-dependence: they bound c times a positive power
of the window, never c. ULS: h = P (spectrum {+1,-1}, ||h|| = 1, v_LR bounded uniformly in d) and c = d - 1,
so no bound in terms of ||h||, v_LR, or the eigenvalues of h can exist. Multiplicities (i.e. d) must enter.

## R3. Representation theory (checked: lie.py table in scratchpad, reproduced below)
G_k WZW with on-site G acting on C^d via R, c = k dim G/(k + h^v). ind(R) = embedding index of R in SU(d).
Free partons (one fermion per site, d flavours, Gutzwiller-projected) give SU(d)_1 which contains G_{ind(R)}.
GKO coset unitarity: c(G_{ind R}) <= c(SU(d)_1) = d - 1. So "IR = G_k with k <= ind(R)" implies c <= d - 1.
k_max(G,R) := max k with c(G_k) <= d - 1:
SU(N) fund: k_max = 1 = ind (conformal embedding).  SO(N) vector: 2 = ind (SO(N)_2 in SU(N)_1).
Sp(2n) fund: 1 = ind.  SU(2) spin 1: 4 = ind (SU(2)_4 in SU(3)_1).  SU(2) spin >= 3/2: no constraint (c < 3 <= d - 1).
SU(N) antisym2: N - 2 = ind (conformal).  E6 27: 6 = ind (conformal).  G2 7: k_max 3 > ind 2.  F4 26: 8 > 6.
For conformally embedded R the conjecture restricted to WZW IR is EQUIVALENT to the parton level bound k <= ind(R).

## R4. Non-integer d (checked by hand)
SU(2)_k spin-1/2 anyon chain, AFM: c = 1 - 6/((k+1)(k+2)), quantum dimension d_q = 2 cos(pi/(k+2)).
c - (d_q - 1) = 2 - 2cos(pi/(k+2)) - 6/((k+1)(k+2)) > 0 for all k >= 2 (k=2: 0.086, k=3: 0.082, large k: (pi^2-6)/(k+2)^2).
c vs log2 d_q: k = 2 equality (1/2), k >= 3 c > log2 d_q (golden: 0.7 > 0.694).
FM golden chain c = 4/5 > 0.694 too. Majorana chain with one Majorana per site: "d" = sqrt 2, c = 1/2 > sqrt2 - 1.
So "c <= d - 1" with d -> quantum dimension is false; the bound needs the integer tensor-product d.
