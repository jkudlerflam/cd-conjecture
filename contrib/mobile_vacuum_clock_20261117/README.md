# Mobile vacuum-clock control with moderate coefficients

Date: 2026-11-17. Author and AI assistance: OpenAI Codex, GPT-6 family.
No independent agent or human review is claimed.

## Result

**[Proved, awaiting independent review.]** The literal d=22 model has a
positive, volume-independent compressed vacuum-clock pin cost, retaining
the terminal selector term, its actual contact halves, marker energies,
surrounding interacting cells and the gapless vacancy dynamics.

A simple sufficient allowance is

    0<q<=1/4,  nf<=1/64,  0<J/t<=1/80,
    t/2<=s<=3t/4.

Here s is the positive length tilt used in the previous mobile result.
The theorem covers that ensemble and its canonical minimizing
vacancy-number blocks. It does not assume equivalence at every
prescribed density. No charged-mass factor D or growing cell-length
factor multiplies J/t in this clock estimate.

For the normalized frustration-free clock
u proportional to (1,q,...,q^n), pin a selected cell to the vacuum space
orthogonal to u. Its exact tilted ground energy satisfies

    F_clock-pin-F >= J [d_n/257-tau_n] > 0,
    d_n=1+q^2-2q cos(pi/(n+1)),
    tau_n=q^(2n+2)/(1+q^2+...+q^(2n)).

The uniform lower bound is 15J/69904. At q=1/4 it improves to
2J/4369 for all n>=1. These are pin costs, not a global spectral gap
inside the critical payload or charge fluid.

**[Proved excitation estimate.]** A second variational comparison bounds
the vacuum-clock excitation weight in the full interacting ground state,
including nonzero ports. Under the same schedule and the already
sufficient frozen allowance

    B=sqrt(2+2 eta)|lambda|
      <=(1-q^2)^2 q^(2n+4)/10,

this weight is less than 1/2, uniformly in volume and n, in the rooted
tilted model. The sharper parameter-dependent formula is in
[PROOF §6](PROOF.md). Charged-rail occupation is a separate quantity;
this bound does not by itself show that most cells lie in the full code.

For physical cyclic rings the proof retains one global vanishing error,
with the explicit occupation formula in PROOF §7. It is not charged
per cell or per defect. These statements hold in the exactly invariant
legal material sector at finite K; illegal competitors remain unexcluded.

## What makes the estimate work

The contact-weighted frustration-free edges all annihilate the same u.
Their common projector is therefore **exactly invariant** before the
terminal term, even with fluctuating contacts and interacting surroundings.
The marker contact energy commutes with it.

The physical q^2 terminal term breaks that invariance. Its clock leakage
has norm at most J q^2 sqrt(p_n(1-p_n)), with
p_n=q^(2n)/(1+...+q^(2n)). At n=1, q=1/4 this is J/68. Sending J/t
to zero does not remove the leakage relative to J. No other fixed
clock-only rank-one projector is invariant under all the actual contacts.

The proof uses a quantum channel that resolves **coherent clock
gradients**, rather than measuring a classical head position. It removes
a vacancy at the corresponding edge and retains the actual joint
contact/gradient expectation. A Cauchy--Schwarz estimate controls the
kinetic change. The affected material terms have norm bounded by

    c_n+q^2+|mu_sel|/2+nf,
    c_n=1+q^2+2q cos(pi/(n+1)).

Charged masses vanish on the active vacuum branches and are unchanged
elsewhere. This is why D drops out. Replacing a pinned cell by u removes
its edge energy and adds at most J tau_n in terminal energy. All unchanged
bulk terms cancel exactly. There is no contact-average substitution,
vanishing charge-gap denominator, or extensive O(J^2 L/t) error.

## Sufficient bounds versus small-chain observations

The checks use q=1/4, two selected cells, and half occupancy. At
J/t=1/80, lambda=f=0:

| n | Proved pin cost divided by J | Observed finite-ring pin cost divided by J | Observed clock-excitation weight |
| --- | ---: | ---: | ---: |
| 1 | 0.000458 | 0.3193 | 0.000174 |
| 2 | 0.002933 | 0.2475 | 0.0000146 |

The finite tests retain positive pin costs through J/t=1. No instability
was found in that tested range. This does **not** establish a uniform
thermodynamic threshold: finite-dimensional variational pin positivity
alone is much weaker. The proof is conservative by large factors at
these low-n points. We have not responded to that looseness by adding
penalties or requiring a smaller q.

The parameter-dependent sufficient J/t allowances at q=1/4, f=0 are
about 0.0553 for n=1 and 0.0281 for n=2; 1/80 is a simple common
allowance valid at every n. The theorem's physical-ring error disappears
only after the thermodynamic limit at fixed nonzero parameters.

## Simplicity audit

[PARAMETERS.md](PARAMETERS.md) lists the independent controls, actual
physical scales and term counts, and distinguishes construction
requirements from conservative proof allowances.

Two changes were checked without adding states or interactions:

- **Delta=4 exactly** preserves grammar, selector and the finite-order
  two-charge mechanism while removing the tiny +E_n calibration. It
  changes one coefficient; it does not remove a physical term.
- **D=0** genuinely removes the extra express/terminal onsite projector
  and preserves grammar, the selector and this clock theorem. It loses
  the existing large-D band proof, so it is not a certified simpler
  replacement for the proposed high-c family.

Deleting the q^2 terminal term would restore an exact clock code but
destroy finite-length selection. Material identity and coboundary
shifts are generally physical after the vacancy lift; they cannot be
used to make its term count look smaller.

This architecture still has many processes and delicate band scales.
The new result removes an unnecessary mass hierarchy from **one proof**;
it does not make the whole construction a simple moderate-ratio example.

## Checks and reproduction

**[Numerical corroboration; not interval certificates.]**
[checks.py](checks.py) independently assembles the contact-weighted
clock sources and compares them with the literal bond at n=1,2. It
checks exact common-kernel invariance, the terminal leakage identity,
and the coherent-gradient channel with noncommuting edge operators.
It also retains the finite pin scans, full-payload/nonzero-port checks,
and literal Delta=4 and D=0 variants, including their Gaussian pair
coefficient checks. Outputs and versions are in [checks.json](checks.json).
The largest full material space has dimension 27,440, the largest saved
eigenpair residual is below 7.3e-12, and the recorded run took 1.4 seconds.

An initial added full-payload test passed letter tuples to an integer
basis assembler; that representation mismatch was corrected before the
final run. It did not change an inequality or a reported result.

From the repository root, with NumPy, SciPy and mpmath installed:

```sh
OPENBLAS_NUM_THREADS=1 python3 -B contrib/mobile_vacuum_clock_20261117/checks.py
```

The script writes only its adjacent small JSON. The previous frozen and
mobile contributions are imported as literal-bond cross-checks and are
left unchanged. No external physics premise is needed for these
elementary variational and finite-matrix arguments.

## What remains

The bare clock projector is not invariant in the true selector model.
The positive pin cost and controlled local excitation weight do not
yet construct a dressed product code or control its dynamical influence
on arbitrarily small payload energy scales. That relative dynamical
matching is the remaining vacuum-clock task.

Wrong material lengths, illegal sectors and the complete interacting
multichannel CFT remain separate. The all-n clock estimate at q=1/4
is not an all-n band certificate. Proposed total nonchiral Virasoro
central charge c=n+1 remains conditional; no microscopic c is measured
or established here.
