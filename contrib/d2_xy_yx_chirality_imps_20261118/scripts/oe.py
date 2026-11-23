import glob,re,numpy as np
R={}
for fn in glob.glob("out_obc/tg_*_rot_L*_p*_chi*.npz"):
    m=re.match(r"out_obc/tg_d([\d.]+)_h([\d.]+)_rot_L(\d+)_p(\d)_chi(\d+).npz",fn)
    d,h,L,p,c=m.groups(); R.setdefault((d,h,int(L)),{}).setdefault(int(p),{})[int(c)]=float(np.load(fn)["E"])
for (d,h) in sorted(set((k[0],k[1]) for k in R)):
    row=[]
    for L in sorted(L for (dd,hh,L) in R if (dd,hh)==(d,h)):
        e=R[(d,h,L)]
        if 0 in e and 1 in e:
            cm=max(set(e[0])&set(e[1]))
            prev=sorted(set(e[0])&set(e[1]))
            gap=e[1][cm]-e[0][cm]
            dg=(gap-(e[1][prev[-2]]-e[0][prev[-2]])) if len(prev)>1 else np.nan
            row.append("L=%d chi=%d dE=%+.2e L*dE=%+.4f (chi-chg %+.1e)"%(L,cm,gap,L*gap,dg))
    if row: print("d=%s h=%s  "%(d,h)+" | ".join(row))
