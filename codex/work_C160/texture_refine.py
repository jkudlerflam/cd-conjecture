"""C160: refine the flagged cases of texture_periodic.py (off-harmonic weight > 1e-2).
(a) refine q1 by maximizing the weight on its harmonics (narrow windows); (b) rerun with 16 starts and N = 600;
(c) local wavelength: spacing of successive crossings of the dominant Fourier component, to see defects."""
import numpy as np, texture_periodic as tp
rng = np.random.default_rng(7)
def harm_weight(n, q, width_fac=1.0):
    N = len(n); b = n[N // 4: 3 * N // 4]; M = len(b)
    j = np.arange(M); win = np.hanning(M)
    bb = (b - b.mean(0)) * win[:, None]
    tot = np.sum(np.abs(np.fft.fft(bb, axis=0)) ** 2) / M      # Parseval normalization over all k
    w = 0.0
    for m in range(1, 40):
        ph = np.exp(-1j * m * q * j)
        w += 2 * np.sum(np.abs(ph @ bb) ** 2) / M if (m * q) % (2 * np.pi) > 1e-9 else 0
    return w / tot
for trial in range(16):
    eps = 0.12
    D = rng.normal(size=3); D /= np.linalg.norm(D)
    S = rng.normal(size=(3, 3)); S = (S + S.T) / 2; S -= np.trace(S) / 3 * np.eye(3)
    b = rng.normal(size=3); sc = rng.uniform(0.1, 0.8)
    if trial not in (2, 3, 5, 6, 11):
        continue
    J = -np.eye(3) + eps * tp.A(D) + sc * eps ** 2 * S; B = sc * eps ** 2 * b
    for N, st in ((400, 6), (600, 16)):
        tp.rng = np.random.default_rng(100 + trial)
        E, n = tp.ground(J, B, N, eps, starts=st)
        kind, q1, off = tp.periodicity(n)
        qs = np.linspace(q1 - 0.01, q1 + 0.01, 401)
        ws = [harm_weight(n, q) for q in qs]
        qb = qs[int(np.argmax(ws))]
        # local wavelength from zero crossings of the bulk component along the dominant direction
        bulk = n[N // 4: 3 * N // 4]; c = bulk - bulk.mean(0)
        u = np.linalg.svd(c, full_matrices=False)[2][0]
        s = c @ u; zc = np.where(np.diff(np.sign(s)) > 0)[0]
        lam = np.diff(zc)
        print(f'trial {trial} N={N} starts={st}: E/N={E/N:+.7f} q1={qb:.5f} weight on harmonics={max(ws):.4f} '
              f'local wavelengths {lam.min() if len(lam) else 0}-{lam.max() if len(lam) else 0} (2pi/q1={2*np.pi/qb:.2f})')
