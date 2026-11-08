"""C163: K_S(12), K_S(16) at the C159 stage-B points (Claude Opus agent, 2026-11-08). Output kest_c159.out."""
import json, glob, numpy as np
from kest import k_est
names = ['hi_TFM_59', 'hi_TFM_60', 'hi_TFM_62', 'hi_TFM_63', 'hi_TFM_64', 'hi_TFM_65', 'hi_TFM_66', 'dist_TFM_77', 'dist_TFM_78',
         'hi_T_34', 'hi_T_36', 'hi_T_37', 'hi_T_39', 'ising_T_43', 'ising_T_46', 'ising_T_41', 'hi_T_40', 'dist_T_50', 'hi_T_35', 'hi_T_38',
         'dist_T_49', 'dist_T_52', 'hi_K_1', 'dist_K_18', 'dist_K_19']
for n in names:
    d = json.load(open(f'c159_refine/{n}.json'))
    v = d['scan']['vec']
    K12 = k_est(v, 12)[0]; K16, info = k_est(v, 16)
    L = d['dmrg']['L']; sm = [L[k]['Smid'] for k in ('32', '64', '128')]
    sl = [6 * (sm[1] - sm[0]) / np.log(2), 6 * (sm[2] - sm[1]) / np.log(2)]
    print(f"{n:12s} parentDelta={d['scan']['Delta']:+.2f} harm={d['scan']['harm']:.2f}  K_S(12)={K12:.2f} K_S(16)={K16:.2f}  "
          f"|mag_m|={abs(np.dot(info['m'], info['mag'])):.3f}  DMRG c_prof(128)={L['128']['c_prof']:.3f} slopes={sl[0]:.2f},{sl[1]:.2f}", flush=True)
