# An entropy-sensitive contact bound for the mobile construction

Date: 2026-11-09. AI-assisted derivation by OpenAI Codex (GPT-6 family),
independently reviewed by the coordinating Codex agent (AI only).
No human review.

**Status: [proved, checked by a separate AI agent; conditional only on the
stated Hamiltonian hypotheses].**
This strengthens the legal-space contact and pin bounds of
[C92](../../codex/work_C92/PROOF.md), without adding assumptions. It does
not establish finite-penalty matching, wrong-length selection, a dressed
clock gap, or a central charge.

## Statement

Use the periodic, independently labelled gap representation of C92, with
at least three material sites:

\[
 H=-t\sum_j(E_j^+E_{j+1}^-+E_j^-E_{j+1}^+)
   +s\sum_j g_j+J\sum_j A_jh_j,
 \qquad A_j=|0\rangle\langle0|_j.
\]

Here each gap ranges over the nonnegative integers; the shifts have norm
one; the material internal space is finite at fixed volume; the material
operators may fail to commute with each other, but act independently of
the gap coordinates and satisfy \(\|h_j\|\le M\). Take \(t>0\),
\(J\ge0\), \(s\ge s_{\min}>JM\), and write
\(c=s_{\min}-JM\). The Gibbs state includes every gap number, with
positive length tilt \(s\). It is not a state at fixed vacancy number.

For the contact expectation \(a=\langle A_j\rangle\) at any site,

\[
 c(1-a)\le4t\sqrt{a(1-a)}+T h_2(a),\qquad
 h_2(a)=-a\log a-(1-a)\log(1-a).                 \tag{1}
\]

The logarithms are natural. The endpoint values of \(h_2\) are defined
by continuity. In particular, the same source and temperature window
used by C92,

\[
 s\in[t/2,3t/4],\qquad JM\le t/4,\qquad
 0<T\le t/(16\log2),                            \tag{2}
\]

implies the uniform bound

\[
 \boxed{\langle A_j\rangle>1/264.}              \tag{3}
\]

For zero-temperature limits of these Gibbs states, a stronger statement
holds:

\[
 \langle A_j\rangle\ge {c^2\over c^2+16t^2},
 \quad\text{hence}\quad\langle A_j\rangle\ge1/257
 \quad\text{under (2)'s source assumptions}.     \tag{4}
\]

The bounds are independent of volume, cell lengths, pin pattern and
internal excitation gaps. For \(t=0\), (1) remains valid with its kinetic
term omitted; (2) is intended only for \(t>0\).

## Proof of (1)

Retain C92's vacancy-removal channel at site \(j\):

\[
 \Phi(\rho)=P\rho P+S^*\rho S,\qquad
 P=A_j,\quad S=E_j^+.
\]

The two Kraus operators are \(P,S^*\), and
\(P+SS^*=1\), so the channel is trace preserving. Its Stinespring
environment is two dimensional. Its diagonal probabilities in the Kraus
basis are exactly

\[
 \operatorname{Tr}(P\rho P)=a,\qquad
 \operatorname{Tr}(S^*\rho S)=1-a.
\]

Off-diagonal environment entries need not vanish. Dephasing can only
increase its entropy, so \(S(\rho_E)\le h_2(a)\). Subadditivity after
the Stinespring isometry gives

\[
 S(\rho)-S(\Phi(\rho))\le S(\rho_E)\le h_2(a). \tag{5}
\]

This replaces C92's state-independent upper bound \(\log2\) on the
same entropy loss. No assumption about a classical material source or a
factorized Gibbs state is made.

C92's operator compression identities and Cauchy--Schwarz estimate give
the unchanged energy inequality

\[
 \operatorname{Tr}H[\Phi(\rho)-\rho]
 \le-c(1-a)+4t\sqrt{a(1-a)}.                    \tag{6}
\]

For clarity, the kinetic constant comes from the two hopping bonds
meeting \(j\), whose off-diagonal block satisfies
\(\|PH(1-P)\|\le2t\). Removing a vacancy lowers the tilt by at least
\(s\) on the nonempty subspace and changes the material contact energy
by at most \(JM\).

The Gibbs variational principle and (5)--(6) imply

\[
 0\le F_H(\Phi(\rho))-F_H(\rho)
 \le-c(1-a)+4t\sqrt{a(1-a)}+Th_2(a),
\]

which is (1). Finite volume and positive tilt make the gap heat trace and
its first moment finite. The vacancy-removal channel never increases the
gap number, so the same total-gap cutoff argument as C92 justifies the
identities and entropy bound before removing the cutoff.

## Exact rational verification of (3)

Suppose \(a\le1/264\). Under (2), \(c\ge t/4\). On this interval,
each of \(-\tfrac14(1-a)\), \(4\sqrt{a(1-a)}\), and \(h_2(a)\)
increases with \(a\). Therefore the right-hand side of the final
inequality above, divided by \(t\), is at most

\[
 -{263\over4\cdot264}
 +{4\sqrt{263}\over264}
 +{h_2(1/264)\over16\log2}.                    \tag{7}
\]

The elementary estimates
\(-(1-a)\log(1-a)\le a\), \(\log264<6\) and
\(\log2>2/3\) give

\[
 {h_2(1/264)\over16\log2}
 <{21\over32\cdot264}.
\]

Thus (7) is strictly less than

\[
 {128\sqrt{263}-2083\over32\cdot264}<0,
\]

because
\(2083^2-128^2\cdot263=29897>0\). This contradicts the variational
inequality. The logarithm estimates can themselves be verified without
floating point: \(e^6>\sum_{r=0}^{7}6^r/r!>264\), and
\(\log2=2\operatorname{artanh}(1/3)>2/3\).

For (4), take \(T\to0\) in (1), using \(h_2(a)\le\log2\).
If \(a<1\), division by \(\sqrt{1-a}\) and squaring give the stated
bound. The case \(a=1\) already satisfies it.

## Consequence for the d=13 pin comparison

Import precisely the legal, infinite-grammar-penalty pin replacement of
C92 equations (10)--(17), including its nonzero-port allowance and its
sum over payload charge. For a length-15 pinned cell,

\[
 \delta_{\rm pin}
 =J\left[2a_*\cos(\pi/16)-4\sqrt2|\lambda|\right].
\]

For the C88 port \(\lambda=2^{-608}\), use \(a_*=1/264\). Since
\(\cos(\pi/16)>49/50\),

\[
 {2\cos(\pi/16)\over264}-{1\over135}
 >{49\over6600}-{1\over135}={1\over59400}
 >6\,2^{-608}>4\sqrt2\,2^{-608}.
\]

Consequently

\[
 \boxed{F_{\rm pinned}-F_{\rm unpinned}
 \ge k\left[{J\over135}-T\log24576\right].}     \tag{8}
\]

Here the unpinned sector has the same cell lengths and the same retained
gap ensemble; all other cell lengths and quantum material states remain
arbitrary. Summing over the length-15 pin flags yields

\[
 0\le p_{\rm all\ pin\ flags}(T)-p_{\rm no\ pins}(T)
 \le{T\over15}\log\left[1+24576e^{-J/(135T)}\right]
 =o(T^2).                                       \tag{9}
\]

This is over thirty times C92's conservative activation-energy lower
bound \(J/4096\). The exact number is a bound, not a measured pin gap.

All ensemble and projection qualifications survive: (8)--(9) do not
bound a contact in a fixed-vacancy-number state, do not eliminate illegal
words at finite \(K\), and do not compare wrong-length sectors. A global
seam remains subextensive, as in C92; it cannot be charged per pin.

## Reproduction

Run `python3 -B contrib/route_audit_20261109/check_contact_bound.py` from
the repository root. It uses only the standard library and checks every
rational inequality used above, the zero-temperature constant, and the
rounded displayed thresholds. Its output is saved in `contact_checks.json`.

The check verifies arithmetic in the analytic proof; it is not a
thermodynamic simulation or a central-charge measurement.

`check_entropy_channel.py` additionally checks the entropy step on 72
finite complex density matrices, including pure states and entangled
gap--material states, and seven diagonal examples saturating the binary
entropy bound. It requires NumPy and saves `entropy_checks.json`.
Saturation for \(\rho=\operatorname{diag}(a,1-a)\) on the gap levels
0 and 1 also shows analytically that \(h_2(a)\) is the sharp universal
entropy-loss bound for this channel at fixed \(a\).
