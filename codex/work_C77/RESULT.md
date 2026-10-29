# C77. What limits the size of a counterexample

2026-10-28. Derivations and conventions are in [PROOF.md](PROOF.md).

[proved conditional statements; heuristic overall assessment] The strongest
supported conclusion is that ordinary q-component liquids on q+1 states
cannot beat c=d-1. None of the proposed general counting arguments extends
this to all nearest-neighbor chains without an additional completeness or
crossover assumption. The smallest local dimension left open is already
d=2. There is no established counterexample at that dimension in the
supplied results, and this work does not produce one.

## (i) Hard-core species do not give an extra hole mode

[assumed complete liquid; proved count] For q particle species and one hole,

    n_hole + sum_a n_a = 1.

The hole fluctuation is minus the total particle-density fluctuation.
With one regular sea per particle species, there are q canonical bosons.
Spin-charge separation reorganizes them as

    c_total = c_charge + c_spin = 1 + (q-1) = q = d-1.

It does not add spin modes to q already-counted species modes. At the
SU(q)-symmetric fixed point, the spin sector is SU(q)_1, whose Sugawara
central charge is (q^2-1)/(q+1)=q-1. A charge gap, spin gaps, or depleted
species lower the ordinary liquid count. Fixing total particle numbers
does not remove finite-momentum bosons. Strict unit filling does remove
the local charge fluctuation.

[proved bare-band count; assumed continuum applicability] A scalar
nearest-neighbor dispersion is a shifted cosine, with at most one regular
occupied arc and two Fermi points. The two endpoints give one c=1 sector,
not two. Complex hopping shifts the arc. A quadratic band edge has zero
velocity and is not an extra conformal mode. Applying this reasoning to
a strongly projected model requires identifying its complete low-energy
fields first. Projected species operators are not independent canonical
fermions, and the unprojected q-flavor Fock site has dimension 2^q.

[proved correction to the literal model class] Hopping into a hole plus
density interactions alone preserves the ordered color word of the
particles, up to cyclic rotation on a ring. Particles cannot pass.
Therefore this literal model does not automatically have a propagating
SU(q)_1 spin sector. For species-independent hopping and interactions,
an open chain at fixed N has

    H = H_spinless-charge tensor I_(q^N)

in ordered coordinates. Its ordinary liquid has one dynamical charge
boson and frozen color degeneracy. The latter is not q-1 finite-velocity
CFTs. The q=2 infinite-U example is explicit in
[Kumar](https://arxiv.org/html/0812.2977v1).

[assumed established infrared input; proved membership] If particle-color
exchange is included, q is attainable on the same q+1 states. The
antiferromagnetic permutation bond

    h=J sum_(a,b=0)^q X^{ab} tensor X^{ba},  J>0,

has the standard SU(q+1)_1 limit with c=q. It includes hole hopping,
density terms, and occupied-occupied color exchange. The last terms
are absent from the literal hopping-plus-density model. SU(N) saturation
is recorded in [Latorre and Sierra](https://arxiv.org/html/2403.17242v2).
The spin-charge algebra is derived in
[Assaraf et al., Eqs. (19)-(26)](https://arxiv.org/pdf/cond-mat/9903057),
but their finite-U Hubbard model is not itself a q+1-state realization.

[proved limitation] Consequently, “maximum c=q” is justified for the
complete q-mode liquid description and is saturated in the exchange
extension. It is not a theorem established here for every strongly
interacting choice of species-dependent hopping and density couplings.
Word conservation alone also does not prove a global c<=1 theorem for
that entire literal subclass. A universal maximum for that unrestricted
subclass would require further dynamical input. No excess mechanism is
supplied by ordinary spin-charge or Fermi-sea counting.

## (ii) General chains

The status of each argument is as follows. Every implication to an upper
bound on c retains its stated hypothesis.

| Argument | Valid conclusion | Missing ingredient for a general small-d exclusion |
| --- | --- | --- |
| Commuting additive on-site densities | At most d-1 independent charges modulo identity | All CFT sectors must be Gaussian and resolved by these microscopic charges |
| Fermi-point counting | For r canonical free range-one orbitals, c<=r; full Fock dimension is 2^r | A general interacting chain has no bounded-degree Bloch determinant controlling all collective sectors |
| Lieb-Robinson speed | Locally visible propagating modes have speed at most v_LR | A speed ceiling does not bound the number of modes or the conformal temperature window |
| Entropy per site | s(T)<=log d | The conformal slope can persist only below an arbitrarily small crossover temperature |

[proved conditional current statement] r surviving abelian currents give
c=r+c_neutral with c_neutral>=0. Thus their existence gives a lower bound
on c. An upper bound requires proving c_neutral=0 and excluding emergent
unresolved currents. Ising already illustrates a critical sector without
a continuous on-site charge. Nonabelian rank is also insufficient.

[proved conditional thermal bound] There is a normalization issue in the
question. With the usual nonchiral c, a sector with equal left/right
speed v contributes pi c T/(3v) to entropy. The coefficient pi/6 is
correct when summing separately over chiral branches. In either convention,
if at temperature T_0 the conformal entropy has relative error at most
epsilon<1, the speed ceiling and entropy bound give

    c_total <= 3 v_LR log d / [pi (1-epsilon) T_0].

This controls c_total T_0/v_LR. It does not control c_total without a
positive lower bound on T_0/v_LR. Such a bound is absent. Small velocities
make the entropy constraint stronger at a fixed temperature, but need
not leave that temperature within the conformal regime.
[Nachtergaele and Sims](https://arxiv.org/abs/math-ph/0506030) supply the
standard finite propagation-speed input, not the missing thermal scale.

[proved logical countertest] For arbitrary formal c>0 and v<=v_LR,

    tau=3 v log d/(pi c),
    s_trial(T)=log d (T/tau)/sqrt[1+(T/tau)^2]

has the required conformal slope, positive heat capacity, and entropy
bounded by log d at every temperature. Its crossover shrinks as 1/c.
This is a test of the inequalities' logical strength, not a constructed
lattice model. It shows why no numerical bound on c follows from those
inequalities alone. Block-entropy counting similarly retains an unknown
crossover length.

## (iii) Smallest dimension left open and a candidate

[proved conclusion about the available arguments] d=1 is trivial. d=2
is not excluded globally by any argument above. C63, C67, and C68
explicitly leave interacting qubit strata unresolved. The background's
statement that nothing exceeds c=1 at d=2 must be read as absence of an
established example. It cannot be promoted to an exclusion theorem.
There is no basis here for naming d=3 as the first logically open case.

[heuristic search candidate] Use the short Hermitian qubit bond

    h=(1+gamma) X tensor X + (1-gamma) Y tensor Y
       +u Z tensor Z +D(X tensor Y-Y tensor X)
       +(b/2)(Z tensor I+I tensor Z).

For gamma D b !=0, this is an interacting parity-symmetric family without
a continuous additive on-site spin charge. It is the C67 unresolved
family, not a new high-c discovery. A concrete calibrated starting point
is (gamma,u,D,b)=(0,0.4,0.5,0), C68's chiral XXZ liquid. Continue into the
interacting finite-parameter region and look for an independent neutral
Ising transition while the existing boson retains positive speed.
The proposed target is c=3/2>1. The weak neighborhood itself has c=1,
and no location of an additional critical sector is known.

[assumed tests; heuristic reachability] This target must survive every
local invariant mixed perturbation. If a boson vertex V has dimension
x_V, sigma_Ising V has dimension 1/8+x_V and epsilon_Ising V has dimension
1+x_V. Symmetry and momentum must be checked on each product. Parity
alone need not protect them. Marginal energy-current couplings also
need control. No parameter-count argument proves that the requisite
matching coefficients can be tuned to zero.

A convincing candidate would need both Ising and boson spectral towers,
nonzero limiting velocities, and a converged entropy coefficient.
C63 already exhibits an apparent entropy coefficient 3/2 at a z=2,
extensively degenerate boundary. C67's finite-size candidates through
L=20 do not establish c>1. Neither an extra correlation peak nor a large
short-chain entropy fit is enough.

## Checks and completion scope

[numerical] `python3 checks.py` checks the frozen-word spectrum at
q=2,L=6,N=3, a 160-dimensional matrix. It agrees with eight copies of
the 20-dimensional charge spectrum to 1.20e-14. The candidate bond passes
Hermiticity and parity checks, and its uniform charge-commutator map has
full rank at the algebra test point (0.17,0.4,0.5,0.23). Formal entropy
functions pass the stated inequalities for c=1,1.5,10,100. These are finite
algebra and formula checks, not central-charge measurements.

[completion] Files are RESULT.md, PROOF.md, PROGRESS.md, STATE.md,
checks.py, and checks.json. All writes stayed in work_C77. The parent
project memory was not edited because the task explicitly prohibits
writes outside this folder. Local STATE.md provides the handoff.
No counterexample, general hard-core maximum theorem, large computation,
or independent referee review is claimed.
