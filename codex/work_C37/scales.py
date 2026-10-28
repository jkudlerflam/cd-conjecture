"""Numerical display of proved bounds, using logarithms only."""
import math,json
out=[]
for n in range(1,24):
    d=(2*n+1)**2;N=n*(2*n+1);B=32*n
    logQ=1000+9*d*(2000+math.log10(2014))
    logH=n*(logQ+math.log10(B))
    logdelta=-N*(N-1)*(math.log10(2)+logH)
    logv=-n*logQ+(2*n-1)*logdelta-(n-1)*math.log10(B)
    logtau=-n*logQ+2*n*(logdelta-math.log10(64))-(n+4)*math.log10(2)-math.log10(n)-(n-1)*math.log10(B+1)
    logg=(logtau-math.log10(256*n))/2
    a4=(2/(n+2))**2*math.sin(math.pi/(n+2))**4
    p0,p1=1,1004
    for j in range(2,n+2):p0,p1=p1,1004*p1-p0
    logbeta=-math.log10(p1)
    # The loop contribution is negligible for the display, but retained in the proof bound.
    logwide=math.log10(a4*n/16)+2*logbeta
    out.append(dict(n=n,log10_denominator_bound=logQ,log10_momentum_margin=logdelta,
        log10_normalized_speed_bound=logv,log10_analytic_tolerance=logtau,
        log10_conversion_g=logg,log10_alpha4=math.log10(a4),
        log10_tree_hopping_upper_before_g2_lambda2=logwide,
        log10_beta=logbeta))
json.dump(out,open('scales.json','w'),indent=2)
for n in [1,2,3,4,8,12,16,23]:print(json.dumps(out[n-1]))
