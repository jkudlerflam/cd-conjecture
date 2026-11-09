"""Compare the canonical SW source with all-orders finite-ring diagonalization."""
from pathlib import Path
import sys, json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'codex/work_C186'))
import sw_general as sw

f=np.sqrt(.3)/2
rows=[]
for L in (7,8):
    for V in (8.,16.,32.,64.,128.):
        w=(-f*f+1j*1.035/2)/V
        mu=(3-2*f*f-2.35/2)/V
        He,pos,gap,smin=sw.des_cloizeaux(L,V,w,1.,-f,f,mu,0.,True)
        source=He[pos[4],pos[0]]
        rows.append(dict(L=L,V=V,source_re=float(source.real),source_im=float(source.imag),
                         im_times_V5=float(source.imag*V**5),leading_coefficient=2*f*1.035/2,
                         projection_min_sv=float(smin),band_gap=float(gap)))
Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')
