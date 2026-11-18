"""Item 4, step 1: evaluate T, mu^L, mu^R for a member and store them at full precision.
usage: item4_blocks.py n a D eta dps impl [matched]
   impl = author (scalar9.blocks) | referee (readme_scalar.blocks_mp)
   eta may be written 2^-35;  "matched" uses phi_A = pi/(2n) (the certificate's choice) instead of pi/(2(n+1)).
"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.dirname(HERE))
import mpmath as mp
n, a, D, eta, dps, impl = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
mp.mp.dps = dps
matched = len(sys.argv) > 7 and sys.argv[7] == "matched"
a_, D_ = mp.mpf(a), mp.mpf(D)
eta_ = mp.mpf(2) ** int(eta[2:]) if eta.startswith("2^") else mp.mpf(eta)
phiA = mp.pi / (2 * n) if matched else mp.pi / (2 * (n + 1))
d0, d1 = mp.mpf(1), a_ + D_            # d_0, d_1
for k in range(1, n + 1):              # d_{k+1} = (a+D) d_k - d_{k-1}, up to d_{n+1}
    d0, d1 = d1, (a_ + D_) * d1 - d0
beta = 1 / d1
td = eta_ * beta
t0 = time.time()
if impl == "author":
    import scalar9
    T, muL, muR = scalar9.blocks(n, a=a_, D=D_, td=td, phiA=phiA)
else:
    import readme_scalar as rs
    T, muL, muR = rs.blocks_mp(n, a_, D_, td, mp, phiA=phiA)
def dump(M):
    return [[[mp.nstr(M[i, j].real, dps), mp.nstr(M[i, j].imag, dps)] for j in range(M.cols)] for i in range(M.rows)]
out = dict(n=n, a=a, D=D, eta=eta, dps=dps, impl=impl, phiA="pi/(2n)" if matched else "pi/(2(n+1))", beta=mp.nstr(beta, dps), log10_beta=mp.nstr(mp.log10(beta), 15), log10_td=mp.nstr(mp.log10(td), 15),
           seconds=round(time.time() - t0, 1), T=dump(T), muL=dump(muL), muR=dump(muR))
fn = os.path.join(HERE, f"blocks_n{n}_a{a}_D{D}_eta{eta}_dps{dps}_{impl}{'_matched' if matched else ''}.json")
json.dump(out, open(fn, "w"))
print(fn, out["seconds"], "s", "log10 beta", out["log10_beta"], "log10 td", out["log10_td"])
