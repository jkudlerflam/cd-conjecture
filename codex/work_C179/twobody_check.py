"""Cross-checks of the phi = pi/2 inter-valley bound state.
(a) ringk.py (bitmask momentum ED, hard-core bosons) at N = 2, L = 40 to 60, K = pi (L even).
(b) C173's own twobody.py at phi = pi/2 (its K sign convention differs; scan all K).
(c) fine K profile near the intra-valley momenta 2k_a, 2k_b and near K = pi, and the bound-state size."""
import numpy as np, sys
sys.path.insert(0, '../work_C173')
from ringk import reps, build_k
import twobody as tb173
H2 = None
exec(open("twobody_pi2.py").read().split("ks = np.linspace")[0])  # defines H2 without running the scan
ks = np.linspace(-np.pi, np.pi, 40001)
eps = lambda k: 2 * np.sin(k) - 2 * np.cos(2 * k)
emin = eps(ks).min()
for L in (40, 50, 60):
    H, R, P = build_k(L, 2, np.pi, 1j, data=reps(L, 2))
    print(f"(a) ringk L={L} K=pi: E_2 - 2eps_min = {np.linalg.eigvalsh(H.toarray())[0]-2*emin:+.5f}")
e173 = min(tb173.two_body(K, 1j, L=400) for K in np.linspace(-np.pi, np.pi, 97))
print(f"(b) C173 twobody.py, phi=pi/2 |w|=1 L=400: min_K E_2 - 2eps_min = {e173-2*emin:+.5f}")
s = np.arcsin(0.25)
for K in (-2 * s, 2 * s, np.pi - 0.05, np.pi, np.pi + 0.05):
    c = (eps(ks) + eps(K - ks)).min(); e = H2(K, 1j, 800)
    print(f"(c) K/pi={K/np.pi:+.4f}: E_2 - continuum(K) = {e-c:+.6f}")
m = 1 / 7.5
Eb = 0.01953
print(f"    bound-state size estimate 1/sqrt(m_valley * E_b) = {1/np.sqrt(m*Eb):.1f} sites (m_valley = 1/eps'' = 1/7.5)")
