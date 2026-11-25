"""Reproducible scientific figure; inputs are committed numerical results."""
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/quartet-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import HERE


def main():
    out=HERE/'analysis'
    synth=json.loads((out/'synthetic_controls.json').read_text())['results']
    target=json.loads((out/'calibration_sensitivity.json').read_text())['momentum_and_polynomial_filters']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    a,b=axes
    rr=[x for x in synth if x['control']=='exact_XX_times_XX_quartet_d4']
    a.plot([x['width'] for x in rr],[x['signed_fit']['log_p'] for x in rr],'o-',color='#B5472B',label='Fit filtered output to raw power')
    a.plot([x['width'] for x in rr],[x['forward_power_fit']['p'] for x in rr],'s-',color='#176B73',label='Forward-filtered power model')
    a.axhline(4,color='#777777',linestyle='--',label='Known asymptotic p = 4')
    a.set(title='Exact two-boson quartet control (d = 4)',xlabel='Demodulation width (sites)',ylabel='Fitted decay exponent p')
    a.legend(fontsize=8,loc='upper left');a.grid(alpha=.18)
    for op,color,marker in [('2_5','#176B73','o'),('3_6','#9653A1','s')]:
        for degree,style in [(0,'-'),(2,'--')]:
            rows=[x for x in target if x['dataset']=='T1_L384' and x['operator']==op and x['dkappa_pi']==0 and x['degree']==degree and x['forward_fit']]
            b.plot([x['width'] for x in rows],[x['forward_fit']['p'] for x in rows],marker+style,color=color,
                   label=f"Q({op.replace('_',',')}), envelope degree {degree}")
    b.axhline(4,color='#777777',linestyle=':')
    b.set(title='Target L = 384, chi = 1024',xlabel='Demodulation width (sites)',ylabel='Conditional forward-fit exponent p')
    b.legend(fontsize=8,loc='best');b.grid(alpha=.18)
    fig.suptitle('Quartet extraction: filtering bias is large; operator dependence remains',fontsize=13)
    fig.savefig(out/'filter_calibration.png',dpi=180)
    fig.savefig(out/'filter_calibration.svg')


if __name__=='__main__':main()
