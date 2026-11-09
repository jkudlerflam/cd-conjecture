"""C167: digest of two-site VUMPS runs (fesx.py json): per chi, energy, S on both bonds, xi, |staggered moment|,
uniform moment, leading transfer eigenvalues (modulus, phase/pi per two-site cell), and c_pair = 6 dS / dln xi.
usage: python3 ana_path.py file1.json [file2.json ...]"""
import sys, json, numpy as np
for fn in sys.argv[1:]:
    d = json.load(open(fn)); prev = None
    print('==', fn)
    for r in d['results']:
        mg = np.array(r['mag']); st = (mg[:, 0] - mg[:, 1]) / 2; un = (mg[:, 0] + mg[:, 1]) / 2
        S = r['S']; xi = r['xi']; lam = np.array(r['lam']); z = lam[:, 0] + 1j * lam[:, 1]
        cp = f" cpair={6 * (np.mean(S) - prev[0]) / np.log(xi / prev[1]):5.2f}" if prev and xi > prev[1] * 1.01 else ''
        prev = (np.mean(S), xi)
        print(f" chi={r['chi']:3d} e={r['e']:.9f} S={S[0]:.4f},{S[1]:.4f} xi={xi:7.2f} stag={np.linalg.norm(st):.4f} "
              f"unif={np.round(un, 3).tolist()} lam2={np.abs(z[1]):.4f}@{np.angle(z[1]) / np.pi:+.3f} "
              f"lam3={np.abs(z[2]):.4f}@{np.angle(z[2]) / np.pi:+.3f}{cp}")
