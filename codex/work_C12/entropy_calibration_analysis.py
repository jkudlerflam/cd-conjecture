import json,numpy as np

def smooth_fit(row,L,period,alpha=1,margin=0.15):
    key='S' if alpha==1 else 'S2';data=row['entropy']
    ls=np.array([r['cut'] for r in data],float);ss=np.array([r[key] for r in data])
    x=np.convolve(ls,np.ones(period)/period,mode='valid')
    y=np.convolve(ss,np.ones(period)/period,mode='valid')
    use=(x>=margin*L)&(x<=(1-margin)*L)
    xx=np.log(2*L/np.pi*np.sin(np.pi*x[use]/L))
    coef=np.polyfit(xx,y[use],1)
    residual=np.sqrt(np.mean((np.polyval(coef,xx)-y[use])**2))
    return dict(c=float(12*coef[0]/(1+1/alpha)),rms=float(residual),margin=margin)

if __name__=='__main__':
    out={}
    for kind,period in [('xx',2),('su3_refined',3)]:
        rows={}
        for L in [24,48]:
            f='calibration_'+kind+'_L'+str(L)+'.json'
            try:r=json.load(open(f))[-1]
            except FileNotFoundError:continue
            rows[L]=r
            out[f]=dict(sweep=r['sweep'],chi=r['chi'],E=r['E'],half=r['entropy'][L//2-1],fits={str(alpha):[smooth_fit(r,L,period,alpha,margin) for margin in [.125,.1667,.25]] for alpha in [1,2]})
        if len(rows)==2:
            out[kind+'_size_scaling']={key:float(factor*(rows[48]['entropy'][23][key]-rows[24]['entropy'][11][key])/np.log(2)) for key,factor in [('S',6),('S2',8)]}
    print(json.dumps(out,indent=2))
    with open('entropy_calibration_analysis.json','w') as f:json.dump(out,f,indent=2)
