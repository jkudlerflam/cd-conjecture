# Lattice test of the continuum conclusion: is c_A(L) = (6L/pi)(e0 L - E0(L)) nonincreasing in L,
# and is Theta(q) = -L E0(L) convex in q = L^2?  Data: exact Bethe ground energies of the SU(N) Sym^2 chains
# (contrib/casimir_survey/astra_bethe/bethe_exact_energies.json), normalization h = sum_j H_{k-j} P_j, ||h|| = H_k/2.
import json, sys, numpy as np
from scipy.special import digamma
path = sys.argv[1] if len(sys.argv) > 1 else "../casimir_survey/astra_bethe/bethe_exact_energies.json"
r = json.load(open(path))
for f in r['families']:
    N, k = f['N'], f['k']; Hk = sum(1/j for j in range(1, k+1)); hn = Hk/2
    e0 = Hk - sum(digamma(1+(j-1)/N) - digamma(j/N) for j in range(1, k+1))/N
    A_inf = k*(N*N-1)/(N+k)*np.pi/N
    L = np.array([x['L'] for x in f['data']], float); E = np.array([x['E'] for x in f['data']])
    cA = 6*L*(e0*L - E)/np.pi
    Th = -L*E; q = L**2
    slope = np.diff(Th)/np.diff(q)                     # discrete convexity in q: slopes nondecreasing
    conv_viol = np.where(np.diff(slope) < -1e-12)[0]
    mono_viol = np.where(np.diff(cA) > 1e-12)[0]
    print(f"N={N} k={k} d={N*(N+1)//2}  A_inf/||h||={A_inf/hn:.3f}  L={int(L[0])}..{int(L[-1])}")
    print("   c_A(L)/||h|| first five:", np.round(cA[:5]/hn, 4), " last:", round(cA[-1]/hn, 4))
    print("   c_A increases at L =", [int(L[i]) for i in mono_viol][:8], "  Theta(q) convexity fails after L =", [int(L[i+1]) for i in conv_viol][:8])
