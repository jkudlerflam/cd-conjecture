"""Finite loop, exact g^2 matrix entries enclosed by integer arithmetic.
All species have the same scalar loop flux after the tree gauge. Physical
loop phases are (pi, pi/2, pi/2), express neutral phase pi/(2n).
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np, time,json,math
from scipy.linalg import eig
from fractions import Fraction
from rational_bands import integer_ldl,congruence
from loop_scan import margins

def compute(n=23,D=1000,bits=512,eta_power=136,max_seconds=175):
    start=time.monotonic();m=n+1;S=2*n+1;scale=1<<bits
    d0,d1=1,D+4
    for k in range(2,m+1):d0,d1=d1,(D+4)*d1-d0
    detC=d1;loopden=detC*(1<<eta_power)
    steps=detC.bit_length()+bits+20
    assert loopden>100*steps
    den=np.array([D+4]*m+[4]*n,dtype=object);den=den[:,None]+den[None,:]
    def K(X):
        Z=np.zeros_like(X);Z[0,...]-=X[-1,...];Z[-1,...]+=X[0,...]
        Z[:,0,...]-=X[:,-1,...];Z[:,-1,...]+=X[:,0,...]
        return Z
    def solve(Br,Bi):
        Xr=np.zeros(Br.shape,dtype=object);Xi=Xr.copy()
        div=den if Br.ndim==2 else den[:,:,None]
        for k in range(steps):
            Yr=Br.copy();Yi=Bi.copy()
            for Y,X in [(Yr,Xr),(Yi,Xi)]:
                Y[1:,...]+=X[:-1,...];Y[:-1,...]+=X[1:,...]
                Y[:,1:,...]+=X[:,:-1,...];Y[:,:-1,...]+=X[:,1:,...]
            Yr+=K(Xi)//loopden;Yi-=K(Xr)//loopden
            Xr=Yr//div;Xi=Yi//div
            if k%30==0 and time.monotonic()-start>max_seconds:raise TimeoutError((n,k,Br.shape))
        return Xr,Xi
    Br=np.zeros((S,S),dtype=object);Br[0,m]=detC*scale
    Xr,Xi=solve(Br,np.zeros_like(Br))
    Ur=np.zeros((S,S,2*n),dtype=object);Ui=Ur.copy()
    for j in range(n):
        Ur[m+j,:,j]=Xr[m+j,:];Ui[m+j,:,j]=Xi[m+j,:]
        Ur[:,m+j,n+j]=Xr[:,m+j];Ui[:,m+j,n+j]=Xi[:,m+j]
    Yr,Yi=solve(Ur,Ui)
    F=np.zeros((n,n),dtype=object);M=F.copy();EF=F.copy();EM=F.copy()
    for i in range(n):
        for j in range(n):
            f=mi=ef=em=0
            for s in range(S):
                xr,xi=Xr[m+i,s],Xi[m+i,s]
                yr,yi=Yr[m+i,s,n+j],Yi[m+i,s,n+j]
                f+=xr*yr+xi*yi
                ef+=10*(abs(xr)+abs(xi)+abs(yr)+abs(yi))+100
                for r,c,ind in [(m+i,s,j),(s,m+i,n+j)]:
                    xr,xi=Xr[r,c],Xi[r,c];yr,yi=Yr[r,c,ind],Yi[r,c,ind]
                    mi+=xr*yi-xi*yr
                    em+=10*(abs(xr)+abs(xi)+abs(yr)+abs(yi))+100
            F[i,j]=4*f;M[i,j]=-2*mi;EF[i,j]=4*ef;EM[i,j]=2*em
    M=M-M.T;EM=EM+EM.T
    # F and M share denominator 2 scale^2. Exact hopping = -F, onsite = i M
    # in units 2 alpha^4 beta^2, with the stated entry radii.
    np.savez(f'full_loop_n{n}.npz',F=F,M=M,EF=EF,EM=EM,bits=bits,D=D,detC=np.array(detC,dtype=object),eta_power=eta_power)
    print(json.dumps(dict(stage='finite_loop_integer',n=n,seconds=time.monotonic()-start)),flush=True)
    return F,M,EF,EM

def analyze(n=23):
    a=np.load(f'full_loop_n{n}.npz',allow_pickle=True);F,M,EF,EM=[a[k] for k in ['F','M','EF','EM']]
    bits=int(a['bits']);scale=1<<bits
    t=np.load(f'derivative_n{n}.npz',allow_pickle=True);N=t['F_numerator']+t['F_numerator'].T
    R,ds=integer_ldl(N)
    treebits=int(t['bits']);extra=1<<(2*(bits-treebits));scale=1<<treebits
    T=-congruence(F,R,ds,extra=extra)/2;mu=1j*congruence(M,R,ds,extra=extra)/2
    cert=json.load(open(f'certificate_n{n}.json'))
    # Rational dyadic field lying below 0.3 times the certified lower bound.
    field_power=math.ceil(-cert['log10_lambda_min_lower']/math.log10(2)+2)
    field=2.**(-field_power)*congruence(np.eye(n,dtype=object)*(scale*scale),R,ds)
    B=mu+field;I=np.eye(n);Z=np.zeros((n,n))
    zs=eig(np.block([[-B,-T.T],[I,Z]]),np.block([[T,Z],[Z,I]]),right=False)
    roots=np.sort(np.angle(zs));d2,d4,wi=margins(roots)
    centres=(roots[:n]+roots[n:][::-1])/2
    # Rigorous row bounds in the dyadically scaled congruence. Its diagonal
    # tree metric lies in [1,4], so they also bound the whitened matrix norm.
    powers=[(int(d).bit_length()-1)//2 for d in ds]
    E=abs(R)@(2*EF+EM)@abs(R).T
    rowbounds=[sum((Fraction(int(E[i,j]),extra*(1<<(powers[i]+powers[j]+1))) for j in range(n)),Fraction(0)) for i in range(n)]
    rad=max(rowbounds)
    W=congruence(np.eye(n,dtype=object)*(scale*scale),R,ds)
    velocities=[];residuals=[];slopes=[]
    for k in roots:
        h=B+T*np.exp(1j*k)+T.T*np.exp(-1j*k)
        w,U=np.linalg.eigh(h);idx=np.argmin(abs(w));u=U[:,idx]
        hp=1j*T*np.exp(1j*k)-1j*T.T*np.exp(-1j*k)
        slope=np.vdot(u,hp@u).real
        slopes.append(abs(slope));residuals.append(abs(w[idx]))
        velocities.append(abs(slope)/np.vdot(u,W@u).real)
    rec=dict(n=n,status='numerical_roots_with_certified_entry_errors',D=int(a['D']),eta_power=-int(a['eta_power']),
        field_power=-field_power,roots=roots.tolist(),centres=centres.tolist(),delta2=d2,delta4=d4,witness=wi,
        unit_circle_error=float(max(abs(abs(zs)-1))),whitened_entry_error_bound=float(rad),
        max_whitened_root_residual=float(max(residuals)),min_whitened_root_slope=float(min(slopes)),
        log10_min_speed_normalized=float(math.log10(min(velocities))),
        log10_beta=-math.log10(int(a['detC'])))
    print(json.dumps(rec),flush=True);json.dump(rec,open(f'full_loop_bands_n{n}.json','w'),indent=2)
    np.savez(f'full_loop_whitened_n{n}.npz',T=T,mu=mu,field=field,R=R,ds=np.array(ds,dtype=object),field_power=field_power)
    return rec
if __name__=='__main__':
    compute();analyze()
