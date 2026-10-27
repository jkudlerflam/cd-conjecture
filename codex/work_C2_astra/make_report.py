import json,numpy as np
r=json.load(open('results_extended.json'))
text=r'''# Ground-state Koo-Saleur test

All calculations were closed book. The complete sparse sweep covers N=8,10,...,22.
See ks_test.py, results_extended.json, free_ising_check.py, stress_tracking.json,
validation.json, mode_one_results.json, and tables_fits.txt.

## Conventions and analytic result

The Ising and O'Brien-Fendley densities are site-centered,
\[
h_j=-\tfrac12(X_{j-1}X_j+X_jX_{j+1})-Z_j
+\lambda(X_{j-1}X_jZ_{j+1}+Z_{j-1}X_jX_{j+1}).
\]
For spin chains the displayed forward-bond density is used. An additional table
uses the literal Ising summand. This distinction changes both current variance
and the limiting answer.

Write q=4 pi/N and h(q)=sum_j exp(-iqj)h_j. Continuity gives
\[
[H,h(q)]=i(1-e^{-iq})J(q),\qquad
J_{\rm tot}=i\sum_j\sum_{r>0}r[h_j,h_{j+r}].
\]
Only r=1,2 can contribute for these densities. For a fixed high-energy cutoff,
\[
A_{-p}=\lim_{N\to\infty}\frac1N\sum_{a\in\mathrm{UV}}
\frac{|\langle a|J_{\rm tot}|I\rangle|^2}{(E_a-E_0)^p}.
\]
A_0 is the high-energy part of Var(J_tot)/N. Smoothness of the high-energy
spectral measure at q=0 yields
\[
m_0^{\rm UV}=\frac{16\pi^2 A_{-2}}{N}+o(N^{-1}),\qquad
m_1^{\rm UV}=\frac{16\pi^2 A_{-1}}{N}+o(N^{-1}),
\quad c_2\sim\frac{128\pi^2 A_{-2}^3}{A_{-1}^2 N}.
\]
The CFT contributions are O(N^-2) and O(N^-3). The variance alone does not
fix the coefficient: its inverse-energy moments are required. Current
conservation removes these leading UV terms, but need not remove an O(N^-3)
UV contribution to m_1.

The correct low-energy normalization is c_MV=8 w_T/Delta_T^2. Indeed
w_T=(c/2)alpha^2 and Delta_T=2alpha. The target momentum is relative to the
ground momentum, P_T=P_0+4 pi/N. The code selects the largest overlap among
the lowest ten states in that sector.

## Numerical tables

All ground states are unique. Both parities and every momentum were checked
for Ising models. Every momentum in Sz=0 was checked for SU(2) models; every
integer-spin multiplet has an Sz=0 representative. Complex conjugation relates
k and -k. Spin-chain ground momentum is pi for N=2 mod 4 and zero otherwise;
Ising ground momentum is zero. Translation eigenvector residual is below
4.4e-15 and ground eigenvector residual below 1.6e-11. Integrable current
variances are roundoff, below 1e-22. N=8 unrestricted Hilbert-space calculations
agree with the symmetry reduction, and centered Ising results agree with an
independent free-fermion calculation.

'''
for label in ['TFIM','Heisenberg','OF0.2','OF0.35','J1J2']:
 text+='### '+label+'\n\n|N|c2|m0 N²|m1 N³|Var(J)/N|cMV|\n|---:|---:|---:|---:|---:|---:|\n'
 for x in r:
  if x['label']==label:
   a=x['centered'];text+=f"|{x['N']}|{a['c2']:.6f}|{a['m0N2']:.6f}|{a['m1N3']:.6f}|{a['JvarN']:.8f}|{a['cMV']:.6f}|\n"
 text+='\n'
text+='### Literal Ising density\n\n|N|c2|m0 N²|m1 N³|Var(J)/N|cMV|\n|---:|---:|---:|---:|---:|---:|\n'
for x in r:
 if x['label']=='TFIM':
  a=x['literal'];text+=f"|{x['N']}|{a['c2']:.6f}|{a['m0N2']:.6f}|{a['m1N3']:.6f}|{a['JvarN']:.8f}|{a['cMV']:.6f}|\n"
text+=r'''
## Fits

Fit the current moments and variance to C+b/N^2+d/N^4 on N=14,...,22.
The fitted coefficients in m0~a/N, m1~d/N, c2~K/N are

|Model|Var(J)/N limit|a|d|K|
|---|---:|---:|---:|---:|
|OF0.2|0.01403394|0.03803650|0.28566191|0.00539495|
|OF0.35|0.03776324|0.11627186|0.80095630|0.01960181|
|J1J2|0.00743247|0.21006892|0.47547904|0.32802929|

These are extrapolations, not rigorous error bounds. Varying the lower fit
cutoff from 12 to 16 gives K ranges 0.0053914-0.0054034,
0.019546-0.019643, and 0.327868-0.328111. Full q=0 current resolvents were
used to estimate the UV moments; the low-energy current contribution is
assumed to vanish in their thermodynamic limits. Accessible sizes are not
yet in the K/N regime. Fitting c2 directly to K/N would be misleading.

Fits c_MV=C+b/N^2+d/N^4 on N=14,...,22 give

|Model|C|b|d|
|---|---:|---:|---:|
|TFIM|0.4999967|-7.39907|38.6445|
|Heisenberg|1.0022824|23.0368|338.4935|
|OF0.2|0.4979573|9.19727|-37.0018|
|OF0.35, requested largest-overlap rule|0.4073702|79.8833|-6295.2636|
|J1J2|1.0004728|33.2736|681.6547|

The OF0.35 fit is not a central-charge extrapolation. At N=20,22,24 the
largest overlap selects the second state: its c_MV is 0.569018,0.544475,
0.520935. The lowest state instead gives 0.588676,0.576694,0.566675. Its
C+b/N^2+d/N^4 fit on 14,...,24 gives C=0.510962; changing the lower cutoff
from 12 to 16 changes C from 0.516254 to 0.506819. This supports convergence
of the stress-tensor estimate to 1/2 but does not establish that the
largest-overlap prescription has reached its asymptotic state selection.
N=24 was used only for this supplementary sector calculation, not a full
uniqueness sweep.

Direct c2 fits to C+b/N^2+d/N^4 on 14,...,22 give apparent intercepts
0.153083,0.988241,0.0845593,0.181506,0.363814 in the table order. These
intercepts are not trustworthy asymptotic estimates: the generic expansion
has an eventual 1/N leading term, and the exactly solvable centered Ising
control already disproves the naive fitted intercept.

## Exact Ising correction to the premise

The literal density has J_tot=2 sum_j X_j Y_{j+1}, which is not conserved.
It has Var(J_tot)/N=4, A_-2=1/4, and A_-1 -> 8/(3 pi). Thus
\[
m_0\sim4\pi^2/N,\quad m_1\sim128\pi/(3N),\quad
c_2\sim9\pi^4/(32N).
\]
Its largest-overlap low-energy prescription also chooses an energy-family
state rather than the stress tensor; it tends to 8/9.

For the centered density J_tot=sum_j(X_jY_{j+1}-Y_jX_{j+1}) is conserved.
Nevertheless, its high-energy pair form factor has magnitude
2(1-cos(q/2)) |cos(k/2-q/4)|. With antiperiodic momenta
k=(2l+1)pi/N, the Majorana trace gives a factor 1/2 in the sum over pairs.
Away from the low-energy crossing region, |F|^2=q^4 cos^2(k/2)/16+O(q^5).
Consequently
\[
m_0^{\rm UV}=Nq^4/64+o(Nq^4),\qquad
m_1^{\rm UV}=Nq^4/(6\pi)+o(Nq^4).
\]
Here epsilon(k)=4 sin(k/2), and the second expression follows from
(1/2)(q^4/16) times the average of 8 sin(k/2) cos^2(k/2), which is
(1/2)(q^4/16) times 16/(3pi).

Adding the CFT terms, with v=2 and c=1/2,
\[
m_0=4\pi^2/N^2+4\pi^4/N^3+O(N^{-4}),\qquad
m_1=(224\pi^3/3)/N^3+o(N^{-3}),
\]
so
\[
c_2\longrightarrow9/98=0.0918367347,\qquad
c_2=9/98+27\pi^2/(98N)+O(N^{-2}).
\]
The independent free-fermion sum gives c2=0.0918639284 at N=100000.
In contrast c_MV=(1/2)[cos(pi/N)/cos(pi/(2N))]^4 -> 1/2.

For Heisenberg, finite-size c2 is close to one but current conservation alone
cannot prove that limit. A supplementary n=1 calculation, whose CFT vacuum
contribution vanishes, gives m1(n=1)N^3=0.594301,0.628170,0.656830 at
N=20,22,24. It exhibits residual UV weight. Small-system fits cannot settle
the limiting correction, particularly in the presence of marginal logarithms.

## Conclusion

Under the smooth UV spectral-limit assumptions, the non-integrable models
have c2 -> 0, not their central charges. The numerical current moments support
this mechanism, but the accessible sizes do not demonstrate the final 1/N
law directly. Low-energy projection removes the UV sum and gives estimates
consistent with the expected central charges. The largest-overlap rule needs
care near the OF tricritical crossover. Current conservation alone is not a
sufficient condition for the unfiltered ratio to converge to c.
'''
open('ks_report.md','w').write(text)
