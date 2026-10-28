import json,numpy as np,sys

def fit_open(row,L,margin):
    dat=[r for r in row['entropy'] if margin<=r['cut']<=L-margin]
    x=np.log(2*L/np.pi*np.sin(np.pi*np.array([r['cut'] for r in dat])/L))
    return dict(margin=margin,c_vn=float(6*np.polyfit(x,[r['S'] for r in dat],1)[0]),c_renyi2=float(8*np.polyfit(x,[r['S2'] for r in dat],1)[0]))

files=[('calibration_xx_L24.json',24),('calibration_xx_L48.json',48),('calibration_su3_L24.json',24),('calibration_su3_L48.json',48),('dmrg_candidate_L48_d09_mu16.json',48)]
for f,L in files:
    try:
        with open(f) as z:rows=json.load(z)
    except FileNotFoundError:continue
    row=rows[-1]
    half=row['entropy'][L//2-1]
    print(f, 'sweep',row['sweep'],'N',row['N'],'half',half)
    print([fit_open(row,L,m) for m in [L//8,L//6,L//4]])
