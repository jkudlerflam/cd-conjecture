"""Extra calibration: finite XX x XX, momentum error and filter envelopes."""
import csv
import gzip
import json
import numpy as np
from common import BASE,HERE,write_json
from audit_saved import historical,raw_matrix,filter_matrix,eigenvalues,signed_fit,forward_fit


def main():
    out=HERE/'analysis';A=historical();rows=[]
    for name in ('T1_L256','T1_L384','T3_L256','T3_L384'):
        m=json.loads(gzip.decompress((BASE/'meas/quartet'/f'{name}.json.gz').read_bytes()))
        kappa=A.fit_kappa(m['P']['2'],m['refs'][0])
        for key in ('2_5','3_6'):
            for dk in (-.002,-.001,0,.001,.002):
                for degree in (0,1,2):
                    for width in (40,48,56,64):
                        vals=[]
                        for ref in m['refs']:
                            r,C=raw_matrix(m['Q'][key],ref)
                            rc,K,cond=filter_matrix(r,kappa+dk*np.pi,width,2,degree)
                            vals.append((r,rc,K,np.einsum('ar,rij->aij',K,C),cond))
                        n=min(len(v[1]) for v in vals);F=np.mean([v[3][:n] for v in vals],axis=0)
                        rc=vals[0][1][:n];K=vals[0][2][:n];r=vals[0][0]
                        _,g=eigenvalues(F)
                        fit=signed_fit(rc,g,32,56);forward=forward_fit(r,rc,K,g,32,56)
                        rows.append({'dataset':name,'operator':key,'dkappa_pi':dk,'degree':degree,'width':width,
                                     'condition':max(v[4] for v in vals),'signed_fit':fit,'forward_fit':forward})
    finite=[]
    for L in (128,256,384,1024):
        sites=np.arange(1,L+1);modes=np.arange(1,L//2+1)
        U=np.sqrt(2/(L+1))*np.sin(np.pi*np.outer(sites,modes)/(L+1))
        C=U@U.T
        refs=[L//2-17,L//2-8,L//2]
        for width in (32,48,64):
            data=[]
            for i in refs:
                r=np.arange(10,L-i-1)
                j=i+r
                G=(C[i,j]*C[i+1,j+1]-C[i,j+1]*C[i+1,j])**2
                rc,K,cond=filter_matrix(r,.2077*np.pi,width,2)
                if not len(K):continue
                data.append((r,rc,K,(K@G).real))
            if not data:continue
            n=min(len(v[1]) for v in data);r,rc,K,_=data[0];rc=rc[:n];K=K[:n]
            g=np.mean([v[3][:n] for v in data],axis=0)
            finite.append({'L':L,'width':width,'refs':refs,'signed_fit':signed_fit(rc,g,32,56),
                           'forward_infinite_power_fit':forward_fit(r,rc,K,g,32,56)})
    write_json(out/'calibration_sensitivity.json',{'momentum_and_polynomial_filters':rows,
        'exact_finite_XX_times_XX':finite,'finite_control_note':'Even OBC half-filled free chains, tA=1,tB=0.7, d=4. Exact Fermi projector and Wick determinant squared. Known asymptotic p=4; pure infinite power fitting need not remove boundaries.'})
    print('calibration variants',len(rows),'finite controls',len(finite))


if __name__=='__main__':main()
