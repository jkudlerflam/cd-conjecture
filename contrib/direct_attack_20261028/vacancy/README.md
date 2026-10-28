# A gap-free, spatially weighted ultraviolet normal form

Status: **[proved analytic lemma; full vacancy matching remains open]**.

This proves a controlled Gaussian-filter variant of the polynomial-locality task in
[`uniform_matching`, Section 4](../../counterexample_audit_20261028/uniform_matching/README.md).
The strongest construction uses a single exponential of a truncated,
strictly local generator. An earlier finite-depth circuit proof is also
retained. The useful point is a finite power of the frequency cutoff,
not an optimal exponent.

For a bounded finite-range chain `H=t H0+J V`, let `s=J/t` and let
`0<delta<1` be a frequency cutoff in units of `t`. There is a local unitary
change of basis such that

    U H U* = t H0 + J V_low,local + R,
    ||R||_1 <= C t s^2 delta^-10
                  + C J delta^-5 exp(-c/delta^2).                 (1)

Here `||Phi||_1=sup_y sum_(X contains y)(1+diam X)||Phi(X)||` is a
specified first-spatial-moment interaction norm. The constants depend only
on the original range and bounded local strengths, not on volume or a
spectral gap. `V_low,local` is an explicitly truncated Gaussian filter of
`V` under the `H0` dynamics. Its support radius is `O(delta^-2)`; the
nonlinear remainder has exponentially summable spatial tails. The estimate
requires `C|s|delta^-3<=1/2`, with `C` fixed by the local interaction bounds.

Taking `delta=s^a` with `0<a<1/10` gives `||R||_1/J -> 0`. The error is
small in precisely the spatial norm needed to charge a finite number of
interfaces for common-bulk-subtracted interactions. This result controls
the ultraviolet error only. A small energy-transfer frequency does not
imply a small operator on high-energy bath states. The low-frequency
ground-state comparison, and hence the `d=24` counterexample, remain
unproved.

The filter has Gaussian energy-transfer tails; it is not a strict spectral
cutoff. The single-exponential construction is translation covariant and
avoids the phase dependence of the older color circuit. Hamiltonian
interaction locality alone still does not establish cancellation of
connected ground-state energies. The strongest proof is in
[SINGLE_EXPONENTIAL.md](SINGLE_EXPONENTIAL.md). The original, all-`s`
finite-depth proof with exponent 18 remains in
[UV_NORMAL_FORM.md](UV_NORMAL_FORM.md). Both are valid for
bounded interactions on finite- or infinite-dimensional on-site spaces,
including the phase-model shifts used in the vacancy representation.
It assumes the standard finite-range Lieb-Robinson estimate, whose
applicability is separated from the elementary circuit proof. No gapped
spectral-flow theorem is used.

Author: an OpenAI GPT-6-family Codex subagent in the direct-attack campaign,
2026-10-28. AI-assisted; no human review. No numerical CFT measurement or
resolution of the conjecture is claimed. The coordinating agent and a
separate alternatives agent independently checked the original filter
signs, causal-cone count and weighted-norm exponent. The coordinating
agent then proposed the single-exponential improvement, independently
verified by the vacancy worker using explicit connected-commutator
counting, and subsequently reviewed by the coordinating agent.

Two further finite-cell RG contributions are recorded in
[RG_EXTENSION.md](RG_EXTENSION.md) and
[DEGREE_REPARAMETERIZATION.md](DEGREE_REPARAMETERIZATION.md). They prove
bandwidth-normalized initialization and Taylor bounds, then recover the
published endpoint degree powers by sacrificing fixed-cell smallness.
The latter uses actual analytic perturbative order and finite microscopic
degree; a generic small coefficient norm would not suffice. The complete
multichannel induction and central-charge identification remain open.

`check_normal_form.py` is a compact finite-matrix sanity check of the exact
first derivative, boundary/filter error identity, and quadratic Taylor
scaling. It does not test the asymptotic locality estimate or the vacancy
chain. Run it with Python, numpy and scipy; it writes `checks.json` and took
about two seconds here. The first-derivative identity error is `5.2e-15`
and the remainder divided by `s^2` converges from `3.152656` to `3.152780`
over `s=0.02,0.01,0.005,0.0025`.

The environment used was:

```sh
PYTHONPATH=/private/tmp/cd-campaign-python \
  /Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  -B contrib/direct_attack_20261028/vacancy/check_normal_form.py
```
