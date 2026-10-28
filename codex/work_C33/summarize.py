"""Assemble display scales from already computed certificates and matrices."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,json,math
from structured import visibility
n=23
cert=json.load(open('root_certificate_n23.json'))
bands=json.load(open('full_loop_bands_n23.json'))
a=np.load('full_loop_n23.npz',allow_pickle=True)
F=a['F'];bits=int(a['bits']);dc=int(a['detC']);gp=cert['finite_g_power']
logunit=math.log10(2*visibility(n))-2*math.log10(dc)+2*gp*math.log10(2)
maxF=max(abs(int(x)) for x in F.ravel())
logF=math.log10(maxF)-(2*bits+1)*math.log10(2)
r=dict(n=n,log10_beta=-math.log10(dc),log10_td=-math.log10(dc)-int(a['eta_power'])*math.log10(2),
       log10_g=gp*math.log10(2),log10_E_unit=logunit,
       log10_T_entry_scale_leading=logunit+logF,
       log10_finite_g_hermitian_margin_lower=logunit+cert['log10_finite_g_hermitian_margin_normalized'],
       log10_common_field=logunit+cert['field_normalized_power']*math.log10(2),
       log10_min_speed_leading=logunit+bands['log10_min_speed_normalized'],
       log10_higher_conversion_remainder_upper=math.log10(visibility(n))+4*gp*math.log10(2)-math.log10(4))
print(json.dumps(r,indent=2));json.dump(r,open('scales_n23.json','w'),indent=2)
tree=json.load(open('tree_results.json'))
print('TREE D1000:')
for n in [2,4,8,12,16,23,32,48,64]:
 t=next(x for x in tree if x['n']==n and x['D']==1000)
 print(n,t['log10_T_scale_g2'],t['equilibrated_min'])
