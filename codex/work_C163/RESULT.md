This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C163 result: the attractive liquids at intermediate K

Derivations are in PROOF.md and the log is in PROGRESS.md. Scripts in this folder are kest.py (Luttinger estimator), kest_c159.py, xpi_c159.py, c163lib.py (bond constructions), path_scan.py, refine_big.py, twisted.py and summarize.py. Copied unmodified from earlier runs: est.py, tk.py, dmrg.py, dmrg_c.py (C159/C155), vumps.py (C153), fesx.py and fesfit_c158.py (C158). Typhon inputs are in typhon/. C159's completed outputs are copied to c159_summary.txt, c159_refine/ and c159_scan/.

## Verdict

c <= 1 is not derived on the whole window, and no c > 1 bond was found. This run derives c <= 1 near the twisted curve gamma at twist of order one. It also excludes the bound-state and uniform-tilt routes to a second sector everywhere in the window. The remaining gap is one assumption (A_W), stated below, and it bites only at harmonic amplitudes comparable to the exchange.

## Mechanism (Task 1)

- [proved, C159] Every T and K_I bond has a zero-field spiral parent. A connected liquid there has zero magnetization along the spiral axis.
- [literature, derived] Below the lattice scale, the zero-field attractive XXZ parent contains only the Gaussian c = 1 tower. Its Bethe strings are massless breathers of the same boson. So the parent has no gapped mode waiting to soften. A second sector would need a lattice-scale state to come down, which takes harmonics as large as the exchange.
- [proved, gauge covariance plus Bethe] In the U(1) class every soft wavevector lies in Z alpha + Z pi. The n-string bands are single-harmonic, with minima at n alpha or n alpha + pi.
- Lemma L [derived]: suppose a composite of n host quanta formed a second boson. The conversion vertex e^{i(vartheta_2 - n vartheta_1)} has dimension n^2/(4K_1) + 1/(4K_2). It is allowed by symmetry and is relevant in attractive liquids, so it locks the two phases into one boson. Bound-state second fluids are therefore excluded.
- Proposition P [derived]: to every order in the harmonic amplitude, the infrared is one boson. On top of it sits a devil's staircase of clock-locked tongues at theta in pi Q.
- Lemma U [derived]: attraction favours a uniform Theta-odd Ising mode (a tilt along m). That mode is the boson's own density channel, the partner of the phason (C160 M3), so Ising x boson through it has codimension 2.
- The two routes still open are (a) a staggered Ising mode at momentum pi, independent of the boson, at K > 15/8, and (b) a composite with a band minimum outside Z theta + Z pi, at K > 1. [heuristic] Route (a) is energetically disfavoured, because the parent's m-coupling is ferromagnetic and a nearest-neighbour bond has no frustration. The one known example of (b) is C161's two-valley soliton band. It is preempted there by a first-order transition and is C162's subject.

## Twisted curve (Task 2)

Proposition G [derived; the asymptotic Schrieffer-Wolff expansion is the only input]. In the rotating frame, -R(phi) + eps X becomes h_F plus three parts. The first is a U(1)-class part (XXZ + DM + field along m, with c <= 1 by Bethe ansatz). The second is a set of charge-q terms at momentum q phi, which average out because their intermediate states cost at least eps^{1/2} J (out-of-plane solitons) or J(1 - cos q phi)/q (magnons and strings). The third is the resonant clock terms eps^{ceil(m/2)} cos(m vartheta) at phi = 2 pi p/m. These are relevant at large K and gap the tongues.
- The infrared is one of four things: a ferromagnet, a single liquid with K of order eps^{-1/2}, a Z_m clock-locked phase, or their KT, Pokrovsky-Talapov and z = 2 boundaries.
- Hence c <= 1 on a neighbourhood of every compact piece of gamma with phi in (0, pi).
- At the ends, phi -> 0 and phi = pi (the staggered image), the problem reduces to C160's blow-up. So all of gamma is covered, away from the soliton-onset layers.

## Numerics (Task 3)

- C159's job 17597033 completed on 2026-11-08 [numerical]. On rings, kest.py gives K_S(16) = 1.9 to 7.5 at 16 attractive points (parent Delta from -0.64 to -1.18, harmonic fraction 0.1 to 0.58). That is inside the Ising x boson window. These points read c_prof(128) = 0.98 to 1.03, with flat slopes (L = 128, which is 10 to 30 healing lengths).
- One point meets C159's pre-fixed c > 1 criterion. dist_T_50 has rising slopes 1.20 and then 1.34, but K_S = 0.65 and it sits within 0.2 of XZ. It is a repulsive liquid near the KT point in the XZ crossover window, outside the attractive window.
- hi_T_40 also has a rising slope, 1.21 at L = 128, and K_S of about 0.45. It is outside the window too.
- ising_TFM_74 has c_two = 1.76. It is a Theta-broken ferromagnetic cat with a domain wall, since S(1) = ln 2 and S_mid is flat at L = 32 and 64.
- Typhon job 17600999 [pending; PENDING (Priority) at the last check, then ssh to Typhon dropped]. Directory /scratch/lustre/jkudlerflam/c163/.
  - big/ has DMRG at L = 192 and 256 for hi_T_40, dist_T_50, hi_T_37, hi_K_1 and two displacements of hi_T_40.
  - fes/T40.json and fes/T50.json are two-site VUMPS runs with chi up to 192.
  - twisted/ holds 16 bonds near gamma at L = 64 and 128.
  - paths/ holds 40 paths into the window for T, K_I and N. These use entanglement-peak-first targeting at L = 48, then DMRG at L = 64 and 128 plus K_S(16).
  - The job ends by writing summary.txt. The criteria are in PROOF 4.2, fixed in advance.

## Narrowest remaining assumption (Task 4)

(A_W) [assumed] Take an attractive half-floating liquid of T or K_I (K > 1), with harmonics comparable to the exchange. Then no lattice-scale excitation outside the boson's charge lattice reaches zero energy while the liquid survives. Concretely, (a) holds at K > 15/8 and (b) at K > 1.

It bites for parents with -1 < Delta < -0.67 (case a) or -1 < Delta < 0 (case b), at harmonic fraction above roughly 0.3, away from the neighbourhoods of h_F and gamma. Case (b) also bites at soliton onsets (C162). A_W replaces C159's A6_TK and is narrower than it. In N, Ising x floating boson has codimension 1 at every K, so N stays under C154's (A6).
