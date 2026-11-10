# E1 on C176 effective-model runs cut to a central 48-site window (calibration of the window length).
import json, numpy as np, sys
import c193 as C
for tag in sys.argv[1:]:
    rec = json.loads(open(f'../work_C176/typhon_out/out/{tag}.jsonl').read().strip().split('\n')[-1])
    m = rec['meas']; L = rec['cfg']['L']
    G = np.array(m['Gb_re']) + 1j*np.array(m['Gb_im']); W = G.shape[0]
    lo_full = m.get('gwin', [L//2 - W//2])[0] if isinstance(m.get('gwin'), list) else L//2 - W//2
    a = W//2 - 24; Gc = G[a:a+48, a:a+48]
    res, rmax = C.E1(Gc, L, lo_full + a, C.ranges_E1(48))
    print(tag, 'W_full', W, 'gwin', m.get('gwin'), ' '.join(f"k={p['k']:+.3f} eta(R3,R4)={p['eta'][2]:.4f},{p['eta'][3]:.4f}" for p in res))
