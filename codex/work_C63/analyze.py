import json
from pathlib import Path
import numpy as np
from ed_staggered import sector,low,entropy_profile


def summarize(rows):
    groups={}
    for row in rows: groups.setdefault((row['u'],row['b']),{})[row['L']]=row
    table=[]
    for (u,b),rr in sorted(groups.items()):
        r=[rr[L] for L in sorted(rr)];tail=[x for x in r if x['L']>=10]
        L=np.array([x['L'] for x in tail],float);ener=np.array([x['E0']/x['L'] for x in tail]);vv=np.array([x['v_phonon'] for x in tail])
        design=np.array([np.ones(len(L)),L**-2,L**-4]).T
        ep=np.linalg.lstsq(design,ener,rcond=None)[0];vp=np.linalg.lstsq(design,vv,rcond=None)[0]
        tailS=[x for x in r if x['L']>=12]
        csize=3*np.polyfit(np.log([x['L'] for x in tailS]),[x['S'][-1] for x in tailS],1)[0]
        last=r[-1];Lc=last['L'];vel=last['v_phonon'];g1=last['charge_gaps']['1'];g2=last['charge_gaps']['2']
        table.append(dict(u=u,b=b,c_chord=last['c_chord'],c_size=float(csize),v_extrap=float(vp[0]),casimir_weight=float(-6*ep[1]/np.pi),c_E_assuming_one_velocity=float(-6*ep[1]/(np.pi*vp[0])),K_spectral=float(np.pi*vel/(2*Lc*g1)),charge_ratio=g2/g1,neutral_x_same_T2=last['neutral_same_T2']*Lc/(2*np.pi*vel),max_residual=max(x['residual'] for x in r)))
    return table

if __name__=='__main__':
    rows=[json.loads(s) for s in Path('scan.jsonl').read_text().splitlines()]
    if Path('extended.jsonl').exists():rows += [json.loads(s) for s in Path('extended.jsonl').read_text().splitlines()]
    table=summarize(rows);Path('summary.json').write_text(json.dumps(table,indent=2)+'\n')
    print('u b c_chord c_size cE(one v) v(inf fit) K(Lmax) gap2/gap1 x_neutral')
    for r in table: print(' '.join(f'{r[k]:.6f}' for k in ['u','b','c_chord','c_size','c_E_assuming_one_velocity','v_extrap','K_spectral','charge_ratio','neutral_x_same_T2']))
