"""XX at m = infinity: exact log Z~(L, beta) from C125's crossed product, and the smallest
majorant parameters. Writes data/xx_exact.npz and prints a summary table."""
import os, numpy as np
from majorant import xx_logZt, vmin_needed, F

HERE = os.path.dirname(os.path.abspath(__file__))
M = 0.5
Ls = np.array([2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 1024, 2048])
betas = np.geomspace(0.1, 8192, 60)


def table():
    g = np.zeros((len(Ls), len(betas)))
    for j, b in enumerate(betas):
        g[0, j] = xx_logZt([2.], b)[0]
        g[1:, j] = xx_logZt(Ls[1:], b)
    return g


if __name__ == '__main__':
    path = os.path.join(HERE, 'data', 'xx_exact.npz')
    if os.path.exists(path):
        g = np.load(path)['g']
    else:
        g = table(); np.savez(path, g=g, Ls=Ls, betas=betas)
    LL, BB = np.meshgrid(Ls, betas, indexing='ij')
    print('max over grid of (L/beta) log Z~ :', (LL / BB * g).max(), ' CFT value pi A/6 =', np.pi / 6)
    print('(L/beta) log Z~ at beta = 8192:', dict(zip(Ls[:6], np.round((Ls / betas[-1] * g[:, -1])[:6], 4))))
    for S in [2, 4]:
        for K in [0., np.log(2), 1., 2., 3.]:
            V = vmin_needed(g, LL, BB, S, K)
            i = np.unravel_index(np.argmax(V), V.shape)
            print('S=%d K=%.3f  V\'_min=%.4f (V\'/M=%.3f) attained at L=%d beta=%.3g' % (S, K, V.max(), V.max() / M, Ls[i[0]], betas[i[1]]))
    for K in [0., np.log(2), 1., 2.]:
        kap = ((LL / BB) * (g - K)).max() / M
        print('kappa_eff(K=%.3f) = sup (L/beta)(log Z~ - K)/M = %.4f ; (6/pi) kappa M = %.4f vs A = 1' % (K, kap, 6 / np.pi * kap * M))
