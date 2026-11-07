"""Finite-slab increments, retaining the full thermal history contractions."""
import json
import time
from checks import ROOT
import xx_kernels as xx
import ed_kernels as ed


def main():
    start=time.time()
    out=[]
    for L in [128,8192]:
        for a in [.25,1.,4.]:
            rows=[xx.xx_case(L,a,N) for N in range(1,7)]
            for i,row in enumerate(rows):
                row['model']='XX_free'
                if i: row['increment']=row['G']-rows[i-1]['G']
                out.append(row)
    for model,d,delta,M,L in [('XX',2,0.,.5,6),('XXZ_m09',2,-.9,.5,6),
            ('XXZ_p05',2,.5,.5,6),('Heisenberg',2,1.,.5,6),('SU3',3,1.,1.,3)]:
        small,big=ed.sectors(L,d,delta),ed.sectors(2*L,d,delta)
        for a in [.25,1.,4.]:
            rows=[ed.one_case(small,big,a,N,M,model) for N in range(1,9)]
            for i,row in enumerate(rows):
                if i: row['increment']=row['G']-rows[i-1]['G']
                out.append(row)
        print(model,time.time()-start,flush=True)
    (ROOT/'interval_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('complete',len(out),time.time()-start,flush=True)


if __name__=='__main__': main()
