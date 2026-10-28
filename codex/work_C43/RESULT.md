# C43. One-vacancy completion of C39

[proved, verdict] The literal dimension-47 completion fails. C39's 46-letter object is an initialized open-cell core, not a periodic chain with a spontaneous length selector. Its repeated bond leaves all-left and all-right data strings at zero internal energy. Species-independent vacancy hopping preserves an exponentially large ground space built from these strings. Neither the weak ports nor C39's readout acts on them. This is a finite microscopic obstruction, before the infrared question.

[assumed, inherited fallback; proved dimension] The final retained construction is C39's mobile wrapper with d0=96*46+48=4464. Its c=n+1 conclusion for every fixed finite n remains conditional on the inherited wrapper and infrared inputs. This task establishes neither an unconditional fixed-d counterexample nor every-large-ring uniqueness for that wrapper. It does not rule out a cheaper vacancy completion after a separate autonomous selector is constructed.

[proved, files] The complete derivations are in [PROOFS.md](PROOFS.md), the literal implementation in [microscopic.py](microscopic.py), and the executable checks in [checks.py](checks.py). [checks.json](checks.json) contains the numerical results, and [PROGRESS.md](PROGRESS.md) is the local handoff.

## 1. Exact obstruction and what dimension 47 omits

[proved] C39's literal repeated core bond penalizes an incoming A endpoint and an outgoing B endpoint. Its external left and right boundary operators are necessary for the initialized A, data, B cell. Repeating that bond on a ring does not insert those external boundaries or select n. This is explicit in C39 RESULT and PROOFS Section 7, which count the selector and wrapper separately.

[proved] Include both C39 bus ports with coefficient lambda and the nonnegative weak readout. The occupied bond stays positive when

    3 |lambda| (1+sqrt(eta)) < K,

where K is the grammar penalty. The ports act only in the B,A block, whose unperturbed energy is at least K; each of the three channel products per bus has norm at most one. Outside that block the positive core decomposition is unchanged.

[proved] For an occupied squeezed ring with N>=2, the complete internal zero space is all binary L words and all binary R words. Endpoints are excluded by positive grammar costs, excited tokens by positive mass, and a vacuum data token would introduce a one-way grammar transition around the ring. Thus there are no other zero sectors.

[proved] For 2<=N<L, the exact vacancy-lift ground energy and multiplicity are

    E0(L,N) = -2t sin(pi N/L)/sin(pi/L) - mu N,
    g_N = (2/N) sum_(r=0)^(N-1) 2^gcd(N,r)
        >= 2^(N+1)/N.

[proved] Positivity bounds the full energy below by the zero-twist hard-core charge sea. Tensoring that sea with any cyclic-invariant L or R word attains the bound. Saturation forces the same internal kernel. Counting cyclic-invariant binary words gives the displayed multiplicity. At interior occupied density this degeneracy grows exponentially, including after grand-canonical minimization over N. The desired selected period and unique ground state are absent.

## 2. Squeezing and Ogata-Shiba gluing

[proved] For ordered occupied positions x_1<...<x_N and their internal labels, the exact Hamiltonian is

    H = H_hop[T] - mu N + J sum_j A_j(x) h_(j,j+1),

where A_j tests physical adjacency and the hopping seam rotates the internal labels by T. In an internal cyclic sector T=e^(iQ), the free charge orbitals obey

    k_I=[2pi I+Q+pi(N-1)]/L.

[proved] The last term is the hard-core fermion seam sign, independent of payload parity. At J=0 the charge sector factorizes within each internal twist. At J>0 the adjacency indicators couple position and internal dynamics, so exact additive factorization fails.

[proved] Within a fixed cyclic sector and a chosen charge sea, the first-order compression is J alpha_(L,N) H_s, with H_s=sum_j h_j and

    alpha_(L,N)=N/L-sin^2(pi N/L)/[NL sin^2(pi/L)],
    alpha(rho)=rho-sin^2(pi rho)/(pi^2 rho)>0.

[proved] The proof uses P_Q h_j P_Q=P_Q H_s P_Q/N, not equality of individual rooted adjacency expectations. Cross-twist blocks and degenerate charge seas must be retained. A charge ground-state projection is not a uniform thermodynamic elimination.

[proved conditional on an independently supplied selected period b] If N=bM and cell translation has phase P, the b offset branches give

    Q_r=(P+2pi r)/b,
    W_c=j+Q_r/(2pi)=(k+P/(2pi))/b,
    r=0,...,b-1,    j,k integers.

[proved conditional on that finite cell structure] The period becomes charge twist sectors with fractional winding spacing 1/b, coupled to payload momentum. Low-energy occupied-number changes are multiples of b. Logical fermions retain C27's common seam (-1)^(F-1), where F is their total number, rather than independent flavor seams. These identifications preserve the oscillator count.

[heuristic, standard bosonization] The physical marker harmonic has wavevector 2pi rho/b and dimension 1/b^2 at the free charge reference. Its correlation decays as distance^(-2/b^2), so a selected squeezed crystal becomes fluctuating physical density order. It does not retain b independent frozen physical vacua.

[proved conditional on positive dressed defect cost] Rounding N to the nearest multiple of b costs O(b^2/L), whereas a wrong-period defect costs a fixed positive amount. Thus compatible occupied sectors exist on every sufficiently large physical ring, including rings with incompatible L modulo b. This removes the crystallographic obstruction only under the stated selector and gap assumptions.

[proved, limitation] Gluing does not establish every-ring uniqueness. Charge-sector crossings and equal twist minima can occur. The complex C39 bond has no supplied positivity proof within each interacting sector. For the literal 47-state lift, uniqueness is affirmatively false by the multiplicity formula above.

## 3. Mixed vertices, margins and the all-n schedule

[assumed, continuum dictionary] For an already selected period b, let z be integer payload chiral counts, sum z=0 by C39's exact U(1), and q_pay=p.z. The free payload dimension is |z|^2/2. With occupied counting field y=rho x-phi_c/pi and cell coordinate u=y/b, a local cell harmonic m has

    w=(m+q_pay/(2pi))/b,
    k_phys=2pi rho w mod 2pi,
    Delta_ref=w^2+|z|^2/2.

[proved conditional on this dictionary] At physical umklapp rho w=ell, ell integer, the dimension is ell^2/rho^2+|z|^2/2. At rho=1/2, nonzero ell costs at least four. At ell=0, the resonance condition is the original payload momentum condition. C39's distinct roots and unordered-pair-sum margins exclude every nonforward neutral bilinear and quartic. All remaining non-Gaussian resonances have reference dimension at least three. Fractional marker harmonics cannot be declared relevant uniform perturbations while omitting their physical momentum.

[proved conditional on the dictionary] With delta_int=min(delta2,delta4), every potentially dangerous non-Gaussian vertex below reference dimension three has physical mismatch at least

    min(delta_int/(2b), pi*(2-sqrt(3))) > 0.

[proved conditional on the Gaussian description] A full scaling metric G>=4I/5 in the stated reference field coordinates makes reference-dimension-three vertices have dimension at least 12/5. Forward density and current couplings remain Gaussian. Exact occupied U(1) forbids charge-dual vertices; exact payload U(1) forbids pairing and includes exact parity. No individual flavor U(1) is assumed.

[proved, first-order compatibility] Multiplying the complete occupied bond by J alpha preserves its relative coefficients and Fermi momenta. A core gap of order J alpha and a port of order J alpha lambda produce the intended logical scale J alpha lambda^2. The field must be lifted with the same occupied-pair gate. C39's lacunary choices can therefore be held fixed before choosing the new small ratio J/t.

[assumed, remaining matching input] A successful autonomous occupied model still needs local low-frequency matching that retains the charge boson, preserves the period-defect gap and compactification, and gives corrections tending to zero relative to the slowest fixed-n payload stiffness as J/t tends to zero. Finite-ring denominators do not prove this. The finite free-charge momentum-restoration bound in PROOFS Section 5 supports the separation of nonforward channels but does not establish the interacting matching by itself.

[heuristic, conditional all-n inference] Fix n and C39's nonzero lacunary parameters, reduce g before lambda to suppress its normalized non-Gaussian remainder, stabilize an autonomous selector, and take J/t last. Under the matching and Gaussian conditions, the theory has n+1 bosons. There is no uniform-in-n coupling window. The literal 47-state model cannot use this inference because it has the wrong ground sector.

## 4. Microscopic ED

[numerical, complete Hilbert spaces] Both n=1 and n=2 were tested on the full three-site microscopic chain, including vacancies, all charges, all markers and all illegal words. Every matrix component was diagonalized with multiplicity.

| n | L | Full dimension | Ground energy | Ground multiplicity | Next energy |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | 3 | 103823 | -2.34 | 6 | -2.17 |
| 2 | 3 | 103823 | -2.34 | 6 | -2.17 |

[numerical, parameters] These use t=1, J=.07, occupied mu=.17, a=4, K=10, D=.2, w=.8, v=.65, g=.4, q=.8, eta=.6, lambda=.03 and readout coefficient .03 lambda^2. They are diagnostic parameters, not the lacunary certificate points. The occupied charge commutator vanishes entrywise and its minimum eigenvalue is zero. The ground energies and multiplicities also follow from the exact formulas above.

[numerical, larger invariant microscopic sectors] Additional tests retain every initialized-cell excitation and data state, all cyclic offsets and one moving vacancy. These sectors are closed under the literal physical bond, but are not its unrestricted global ground sectors. For these checks only, readout and occupied mu are zero.

| n | (L,N) | Dimension | Physical vs squeezed entry error | Complete J=0 spectrum error | First-order matrix error |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | (4,3) | 648 | 0 | 8.9e-16 | 1.8e-15 |
| 2 | (5,4) | 3200 | 0 | 1.7e-15 | 3.6e-15 |

[numerical] The interaction's internal-translation commutator norms are 203.548 and 398.427, confirming finite-J coupling between twist factors. Low-energy ED at J=.02,.01 has residual below 3.2e-11 and nonzero order-J^2 corrections. The tests ran in about 3.3 seconds in this environment. They check microscopic identities and the failure, not a central charge or the 4464-state wrapper.

[proved, reproduction] Run `python3 -B checks.py` from work_C43. Only numpy, scipy and the standard library are used. No network, cluster job, subagent, or outside-folder write was used. Project STATE.md was left untouched under the explicit folder restriction.

## 5. Final accounting

[proved] The finite vacancy kinematics and the literal 47-state failure are established. The selected-period harmonic classification is established conditional on the stated bosonization dictionary. No smaller autonomous selector has been constructed, and there is no no-go theorem here against constructing one.

[assumed, inherited construction; proved alphabet count] The fallback onsite space is C^4 tensor [(C^3 tensor C^2 tensor C^46) direct-sum C^3] tensor C^2 tensor C^2, of dimension 4464. It already has the mobile length mode, so another vacancy is unnecessary. Its C39 all-n payload margins and C8/C32 wrapper properties are inherited with their stated qualifications.

[heuristic, final claim] The surviving claim is the conditional fixed-d0=4464 family with c=n+1 for every finite n. The proposed d=47 reduction and its unique-ground-state claim fail for the explicit microscopic reason above. Thermodynamic infrared transfer and every-large-ring uniqueness for the fallback remain unproved here.
