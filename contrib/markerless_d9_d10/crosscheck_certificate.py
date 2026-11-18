"""Cross-check of the integer certificate (certificate/, C84's pipeline with the marker vertex deleted) against the
independent arbitrary-precision scalar formula scalar9.py for the same member:
a = 4, D = 1000, phi_A = pi/(2n), psi_A = pi, psi_B = pi/2, t_d = beta 2^-eta_power, field = beta^2 2^-field_power
(units lambda^2 alpha^4 g^2).  The certificate works in the gauge where the hopping is real, so the two sets of
Fermi roots must agree after one common shift (the phase of the hopping)."""
import mpmath as mp, json, sys, numpy as np
from scalar9 import blocks
from analyze_mp import roots_mp
n, cert_json, dps = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
cert = json.load(open(cert_json)); eta_power, field_power = -cert['eta_power'], -cert['field_power']
mp.mp.dps = dps
a, D = mp.mpf(4), mp.mpf(cert['D'])
d0, d1 = mp.mpf(1), a + D
for k in range(2, n+2): d0, d1 = d1, (a+D)*d1 - d0
beta = 1/d1; td = beta/mp.mpf(2)**eta_power
T, muL, muR = blocks(n, a=a, D=D, td=td, phiA=mp.pi/(2*n), psiA=mp.pi, psiB=mp.pi/2)
mu = muL + muR
th = mp.arg(sum(T[i,i] for i in range(n))); Th = T*mp.expj(-th)
scT = max(abs(Th[i,j]) for i in range(n) for j in range(n))
M = mu.copy(); f = beta**2/mp.mpf(2)**field_power
for i in range(n): M[i,i] += f
zs = roots_mp(M, T); off = max(abs(abs(z)-1) for z in zs)
mine = np.sort(np.array([float(mp.arg(z*mp.expj(th))) for z in zs]))        # roots in the real-hopping gauge
certroots = np.sort(np.array(cert['roots']))
err = min(np.abs(np.sort(np.mod(mine + s + np.pi, 2*np.pi) - np.pi) - certroots).max() for s in (0.0, np.pi))
rep = dict(n=n, D=int(cert['D']), eta_power=-eta_power, field_power=-field_power, dps=dps,
    log10_beta_mine=float(mp.log10(beta)), log10_beta_certificate=cert['log10_beta'],
    hopping_phase_over_pi=float(th/mp.pi), imaginary_part_of_gauged_hopping=float(max(abs(mp.im(Th[i,j])) for i in range(n) for j in range(n))/scT),
    real_part_of_onsite_block_over_T=float(max(abs(mp.re(mu[i,j])) for i in range(n) for j in range(n))/scT),
    my_roots_distance_from_unit_circle=float(off), max_root_difference=float(err),
    certificate_delta2=cert['delta2'], certificate_delta4=cert['delta4'])
print(json.dumps(rep)); json.dump(rep, open(f'crosscheck_certificate_n{n}.json','w'), indent=1)
