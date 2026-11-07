# C142 progress log

Run by a Claude Opus agent (2026-11-07), not by Codex. J=1 throughout.

## Step 0
Read the spec, _Cd_target.md, _bound_background.md, C140 (RESULT, PROOF, REFEREE), C139 and C138 (RESULT, PROOF), C134 and C136 (RESULT, C136 REFEREE). Copied the solvers bethe_scan.py, check_c129.py, continuum_rep.py (repulsive) and bethe_ferro.py, continuum.py (attractive) from C140.

## Step 1, strip (L*)
Found the loss in C140 Lemma 3.3. The tail integral of r0 was bounded by (2/b) r0(Q), while the sharp bound is (pi/(2b)) r0(Q) (the function atan(t)(1+t^2)/t is increasing on (0,1], with value pi/2 at t=1). This gives L r0(x_n) >= b/(2(pi-gamma)) >= b/(2 pi) = 0.159 b, above the needed 0.1324 b at every gamma. Theorem 3.4 then holds with margin theta >= 1 - 2(2-ln2)/pi = 0.168 uniformly, and C140 Theorem 4.1 extends to all 0 <= Delta < 1.
Also: the negative mass beyond x_n is at most (pi-2gamma)/(2 pi L) per side (C140 used the two-sided total 1/L), so the sawtooth is sup|Cum R| <= (1-gamma/pi)/L <= 1/L, not 3/(2L).
Plan: redo C139 Section 8.2 line by line with sawtooth s/L, s <= 1; constants script; numerical checks of each new inequality.

## Step 2, strip (L*) assembled
constants_rep.py: every constant of the redone C139 8.2 as a function of gamma. g = L||G|| <= 45/4, K_large <= 815 analytic (grid sup 811.3, attained as gamma -> 0), so C_rep = 5*815 + pi^2/4 <= 4078 on [0,1]. No L_cut beyond L >= 29 is needed (empty-core case handled by the first-order bound).
check_repulsive.py: 160 cases, every new inequality holds (check_repulsive_summary.json). ed_check.py: Bethe = ED to 2.3e-14 at Delta 0.99, 0.999 (L <= 12); ED defect max 4.000 at Delta = 1, L = 2.

## Step 3, strip (D)
Proved (D) and (D1) up to a factor log(1/eta), uniformly in filling: sigma_L - r_Q = int calR_Q(x,y) dR_L(y), the oscillation bound |int f dR_L| <= Var(f)/L (z_L monotone on the attractive side), and Var_y calR_Q(x,.) <= 2(1+q) calR(0) uniformly in Q by comparison with the full-line kernel through the outside channel, with calR(0) <= (log(1/eta)+pi/2)/(pi eta). So eta L |sigma - r| <= (4/pi)(log(1/eta)+pi/2) everywhere on [-Q_L,Q_L].
Not closed: removing the log needs second-order control on the shell C_0 <= eta L r <= C log(1/eta); (W1)-(W3) not addressed beyond a remark. check_attractive.py runs in parallel per anisotropy.

## Step 4, final
check_attractive_<eps>.json: 95 cases, all pass; true eta L sup|sigma - r_Q| <= 0.317 against proved C_D(eta) = 0.67 to 9.97. PROOF.md (Parts I to III) and RESULT.md written. Outcome: strip (L*) closed with C = 4078 on [0,1]; strip (D) open, reduced to removing log(1/eta) from (D) on the shell C_0 <= eta L r <= C log(1/eta), plus C138's (W1)-(W3).
