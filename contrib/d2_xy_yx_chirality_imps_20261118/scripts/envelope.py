import json,numpy as np,sys
for fn in sys.argv[1:]:
    d=json.load(open(fn))['meas']; L=d['L']; R=d['Gb_rows']
    out=[]
    for r in [2,4,8,16,24,32,48,64,96,128,160]:
        i=L//2-r//2; j=i+r
        if str(i) not in R: continue
        re,im=R[str(i)]; g=abs(complex(re[j-i-1],im[j-i-1]))
        # envelope: max |G| over r-2..r+2 (oscillating)
        env=max(abs(complex(R[str(L//2-rr//2)][0][rr-1],R[str(L//2-rr//2)][1][rr-1])) for rr in range(max(1,r-3),r+4) if str(L//2-rr//2) in R)
        out.append((r,env))
    print(fn, " ".join("%d:%.2e"%x for x in out))
    D=d['Dpair']; k='3_%d'%(L//2-3); v=np.array(D[k]); a=np.hypot(v[:,0],v[:,1])
    print("   pair3 env", " ".join("%d:%.2e"%(r,a[max(0,r-3):r+3].max()) for r in [4,8,16,32,64,96,128,160] if r+3<len(a)))
