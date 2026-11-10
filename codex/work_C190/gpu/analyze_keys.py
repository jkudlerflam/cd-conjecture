import json,sys,math,os
D=sys.argv[1]
def load(t):
    f=os.path.join(D,t+'.json')
    return json.load(open(f)) if os.path.exists(f) else None
def conv(d):
    h=d['hist']; chi=h[-1]['chi']; fin=[x for x in h if x['chi']==chi]; prev=[x for x in h if x['chi']<chi and x.get('Sbar6') is not None]
    a,b=fin[-1],fin[-2] if len(fin)>1 else fin[-1]
    c1=a['res_max']; c2E=abs(a['E']-b['E']); c2S=max(abs(a[k]-b[k]) for k in ('Sbar1','Sbar6','Sbar12')); c2n=abs(a['n_mid']-b['n_mid'])
    c3=a['trunc']; c4=abs(a['Sbar6']-prev[-1]['Sbar6']) if prev else float('nan')
    return dict(chi=chi,ns=len(fin),nsw=len(h),E=a.get('E_phys',a['E']),n=a['n_mid'],S=(a['Sbar1'],a['Sbar6'],a['Sbar12']),c1=c1,c2E=c2E,c2S=c2S,c2n=c2n,c3=c3,c4=c4,
                ok=(c1<=1e-8 and c2E<=1e-8 and c2S<=1e-3 and c2n<=5e-4 and c3<=1e-7 and c4<=0.01), okc13=(c1<=1e-8 and c2E<=1e-8 and c2S<=1e-3 and c2n<=5e-4 and c3<=1e-7))
import glob
KEYS=[(k.split(':')[0],float(k.split(':')[1])) for k in sys.argv[2:]] or [('rat_n316',0.1957),('sym40_n18',0.1326)]
# usage: analyze_keys.py <folder> [key:target_n_mid ...]. Every done run of a (key, L, P) is listed, and the
# lowest-energy run among those passing C1 to C4 is used (among those passing C1 to C3 if none passes C4).
for key,tgt in KEYS:
    print('==',key,'target n_mid',tgt)
    R={}
    for L in (64,128,256):
        for P in (0,1):
          cands=[]
          for f in sorted(glob.glob(os.path.join(D,f'{key}_L{L}_P{P}_*.json'))):
            t=os.path.basename(f)[:-5]; d=load(t)
            if d is None or d.get('status')!='done': print(t,'missing/not done'); continue
            c=conv(d); c['tag']=t; cands.append(c)
            print(f"{t}: chi={c['chi']} fin_sweeps={c['ns']} tot={c['nsw']} E={c['E']:.10f} n_mid={c['n']:.4f} S=({c['S'][0]:.4f},{c['S'][1]:.4f},{c['S'][2]:.4f}) res={c['c1']:.1e} dE={c['c2E']:.1e} dS={c['c2S']:.1e} dn={c['c2n']:.1e} trunc={c['c3']:.1e} C4step={c['c4']:.4f} C1-3={'ok' if c['okc13'] else 'FAIL'} C4={'ok' if c['c4']<=0.01 else 'FAIL'}")
          pool=[c for c in cands if c['ok']] or [c for c in cands if c['okc13']] or cands
          if pool:
            R[(L,P)]=min(pool,key=lambda c:c['E'])
            if len(cands)>1: print(f"   -> used for L={L} P={P}: {R[(L,P)]['tag']}")
    g={}
    for L in (64,128,256):
        if (L,0) in R and (L,1) in R:
            P=min((0,1),key=lambda p:R[(L,p)]['E']); g[L]=R[(L,P)]; print(f" L={L}: ground sector P={P}, E0-E1={R[(L,0)]['E']-R[(L,1)]['E']:+.3e}, n_mid P0/P1 {R[(L,0)]['n']:.4f}/{R[(L,1)]['n']:.4f}")
    for (a,b) in ((64,128),(128,256)):
        if a in g and b in g:
            print(f" ground-branch c_two({a}->{b}) W=1,6,12:", ' '.join(f"{6*(g[b]['S'][i]-g[a]['S'][i])/math.log(b/a):.3f}" for i in range(3)))
        for P in (0,1):
            if (a,P) in R and (b,P) in R:
                print(f"   sector P={P} c_two({a}->{b}):", ' '.join(f"{6*(R[(b,P)]['S'][i]-R[(a,P)]['S'][i])/math.log(b/a):.3f}" for i in range(3)))
    if 256 in g: print(f" C6: n_mid(256) ground {g[256]['n']:.4f} vs target {tgt} -> {'ok' if abs(g[256]['n']-tgt)<=0.006 else 'FAIL'} (both sectors: {R[(256,0)]['n']:.4f}, {R[(256,1)]['n']:.4f})")
