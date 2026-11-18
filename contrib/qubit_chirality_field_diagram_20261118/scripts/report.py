"""Small summaries and standalone scientific figures from saved results."""
import csv
import hashlib
import json
import os
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/qubit-field-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main():
    batches=[json.loads((HERE/'results'/f).read_text()) for f in ('coarse.json','refine.json','larger.json')]
    records=[r for b in batches for r in b['records']]
    summary=[]
    for r in records:
        es=sorted(e for s in r['sectors'] for e in s['energies'])
        summary.append(dict(L=r['L'],h=r['h'],energy_per_site=r['e'],M=r['magnetization'],
                            parity=r['parity'],k=r['k'],gap1=es[1]-es[0],gap2=es[2]-es[0],
                            opposite_parity_gap=r['opposite_parity_gap'],
                            Neel2=r['observables']['Z']['pi_over_L'],
                            dimer2=r['observables']['bond']['pi_over_L'],
                            B2=r['observables']['B']['pi_over_L'],
                            chirality2=r['observables']['kappa']['pi_over_L'],
                            Z_peak_q_pi=r['observables']['Z']['peak_q_over_pi']))
    with (HERE/'results'/'ed_summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]),lineterminator='\n');w.writeheader();w.writerows(summary)
    roots=[]
    for L in (12,16,20):
        rows={r['h']:r for r in records if r['L']==L}
        g0,g1=(rows[h]['opposite_parity_gap'] for h in (1.66,1.68))
        roots.append(dict(L=L,linear_secant_root=1.68-g1*.02/(g1-g0)))
    xs=1/np.array([r['L'] for r in roots])**2;ys=np.array([r['linear_secant_root'] for r in roots])
    rootfits=dict(all_three_intercept=float(np.polyfit(xs,ys,1)[1]),
                  L12_L16_intercept=float(np.polyfit(xs[:2],ys[:2],1)[1]),
                  L16_L20_intercept=float(np.polyfit(xs[1:],ys[1:],1)[1]))
    operator=[]
    for filename in ('operator_levels.json','operator_levels20.json'):
        for r in json.loads((HERE/'results'/filename).read_text())['records']:
            b=next(b for b in r['blocks'] if b['shift']==r['L']//2)
            weights=b['operators']['kappa']['low_level_weights']
            index=next(i for i,w in enumerate(weights) if w>1e-10)
            operator.append(dict(L=r['L'],h=r['h'],gap=b['gaps'][index],
                                 L_gap=r['L']*b['gaps'][index],weight=weights[index],
                                 total_weight=b['operators']['kappa']['total_weight']))
    out=dict(first_crossing_secants=roots,conditional_inverse_L_squared_fit=rootfits,
             caveat='Not thermodynamic error bars; assumes the first odd-parity instability and 1/L^2 corrections.',
             chirality_levels=operator,
             total_ED_batch_seconds=sum(b['runtime_seconds'] for b in batches),
             completed_ED_points=len(records),larger_batch_complete=batches[-1]['complete'])
    (HERE/'results'/'interpretation.json').write_text(json.dumps(out,indent=2)+'\n')

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':170,'font.family':'DejaVu Sans',
                         'svg.hashsalt':'qubit-chirality-field-diagram-20261118'})
    fig=plt.figure(figsize=(12,9.5),layout='constrained')
    gs=fig.add_gridspec(3,1,height_ratios=[1.3,1.5,2.7])
    ax=fig.add_subplot(gs[0]);ax.set_xlim(-6,6);ax.set_ylim(-.4,1.7)
    colors=dict(neel='#78a9cf',high='#82bba1',unknown='#dddddd',critical='#d37f37',two='#a88ac0')
    def band(a,b,color,hatch=None,axis=ax):
        axis.add_patch(Rectangle((a,0),b-a,.55,facecolor=color,edgecolor='#666666',lw=.7,hatch=hatch))
    band(-6,-1.9,colors['high']);band(-1.9,-1.66,colors['unknown'],'////')
    band(-1.66,1.66,colors['neel']);band(1.66,1.9,colors['unknown'],'////');band(1.9,6,colors['high'])
    ax.text(0,.27,'Néel regime\nfinite-size evidence',ha='center',va='center')
    for x in (-4,4):ax.text(x,.27,'Gapped high-field phase\nnumerical support; proved asymptotically',ha='center',va='center',fontsize=9)
    ax.text(0,1.28,'Whole field axis: exact reflection h ↔ −h',ha='center',fontsize=13,weight='bold')
    ax.annotate('unresolved zones near |h| ≈ 1.7–1.9',xy=(1.78,.55),xytext=(0,.94),ha='center',arrowprops=dict(arrowstyle='-',color='#555'))
    ax.annotate('',xy=(-1.78,.55),xytext=(-1.2,.94),arrowprops=dict(arrowstyle='-',color='#555'))
    ax.set_yticks([]);ax.set_xlabel('h  (Pauli normalization, full field on every site)')
    ax.set_xticks([-6,-3,-1.9,0,1.9,3,6]);ax.spines['left'].set_visible(False)
    ax.text(-5.9,-.24,'← −∞',ha='left');ax.text(5.9,-.24,'+∞ →',ha='right')

    bx=fig.add_subplot(gs[1]);bx.set_xlim(1.7365,1.760);bx.set_ylim(-.25,1.7)
    band(1.7365,1.760,'#eeeeee',axis=bx)
    band(1.739,1.741,colors['two'],'///',bx)
    band(1.746,1.7472,colors['critical'],'///',bx)
    for h in (1.738,1.749,1.750,1.752,1.755,1.758):bx.scatter(h,.27,s=38,color='#196c92',zorder=4)
    bx.annotate('c ≈ 1 slope\nκ not detected',xy=(1.738,.55),xytext=(1.7379,1.0),ha='center',arrowprops=dict(arrowstyle='-',color='#666'),fontsize=9)
    bx.annotate('c ≈ 2 slopes\nphase unresolved',xy=(1.740,.55),xytext=(1.7413,1.40),ha='center',arrowprops=dict(arrowstyle='-',color='#666'),fontsize=9)
    bx.annotate('Boson × Ising candidate\ninherited hI = 1.7466 ± 0.0006',xy=(1.7466,.55),xytext=(1.7475,1.0),ha='center',arrowprops=dict(arrowstyle='-',color='#666'),fontsize=9)
    bx.text(1.7548,.98,'Flow toward c ≈ 1\nwith chirality candidate',ha='center',fontsize=9)
    bx.set_yticks([]);bx.set_xlabel('|h|: supplied finite-entanglement window (candidate labels)')
    bx.spines['left'].set_visible(False);bx.set_xticks([1.738,1.741,1.7466,1.750,1.755,1.758])
    cx=fig.add_subplot(gs[2]);audit=json.loads((HERE/'results'/'ladder_audit.json').read_text())
    selected={}
    for r in audit['target']['runs']:
        if r['branch']=='_pin' or r['h'] in (1.739,1.7405) and r['branch']=='_pinB':selected[r['h']]=r
    for key,label,mark in [('xi_chi96','χ ≥ 96','o'),('xi_chi128','χ ≥ 128','s')]:
        points=[(h,r['fits'][key]['c_slope']) for h,r in sorted(selected.items()) if r['fits'][key]]
        cx.plot(*np.array(points).T,marker=mark,ms=4,lw=1,label=label)
    for c in (1,1.5,2):cx.axhline(c,color='#bbbbbb',lw=.8,ls='--')
    cx.axvspan(1.746,1.7472,color=colors['critical'],alpha=.16)
    cx.axvspan(1.739,1.741,color=colors['two'],alpha=.15)
    cx.set(xlabel='h',ylabel='6 dS / d log ξ (fit slope)',xlim=(1.7365,1.760),ylim=(.9,2.5))
    cx.legend(loc='upper right');cx.grid(axis='y',alpha=.15)
    fig.suptitle('Fixed qubit bond: (9/16) XY + (7/16) YX + ZZ, field −h Z',fontsize=16)
    fig.supxlabel('Hatching marks unresolved intervals, not confidence bounds. No target Virasoro c is established.\n'
                  'Néel edge: search band 1.66–1.70 (conditional extrapolation ≈1.668). Upper edge: unlocated beyond 1.758.',fontsize=9)
    fig.savefig(HERE/'phase_diagram.png');fig.savefig(HERE/'phase_diagram.svg',metadata={'Date':None});plt.close(fig)
    svg=HERE/'phase_diagram.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')

    fig,axs=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for L,color in zip((12,16,20),('#3884b8','#d17b29','#688847')):
        rows=sorted((r for r in summary if r['L']==L),key=lambda r:r['h'])
        hs=[r['h'] for r in rows]
        axs[0,0].plot(hs,[r['M'] for r in rows],'o-',label=f'L={L}',color=color,ms=4)
        axs[0,1].plot(hs,[r['Neel2'] for r in rows],'o-',color=color,ms=4)
        axs[1,0].plot(hs,[r['opposite_parity_gap'] for r in rows],'o-',color=color,ms=4)
        target=next(r for r in records if r['L']==L and r['h']==1.7466)
        rr=np.arange(3,L//2+1);cc=np.array(target['observables']['kappa']['connected_re'])
        axs[1,1].plot(rr,(-1.)**rr*cc[rr],'o-',color=color,ms=4)
    axs[0,0].set(xlim=(1.60,1.92),ylim=(-.02,.86),xlabel='h',ylabel='Magnetization per site')
    axs[0,0].legend();axs[0,1].set(xlim=(1.60,1.92),xlabel='h',ylabel='S_Z(π) / L')
    axs[1,0].set(xlim=(1.60,1.92),ylim=(-.005,.24),xlabel='h',ylabel='Opposite-parity gap')
    axs[1,1].set(xlabel='Separation r',ylabel='(−1)^r connected ⟨κ(0)κ(r)⟩',title='h=1.7466: oscillations retained')
    for ax in axs.flat:ax.grid(alpha=.2)
    fig.suptitle('Independent rings: strong size and sector effects near the proposed transition',fontsize=13)
    fig.savefig(HERE/'ed_diagnostics.png');plt.close(fig)
    manifest={}
    for p in sorted((HERE/'scripts').glob('*.py')):
        manifest[str(p.relative_to(HERE))]=hashlib.sha256(p.read_bytes()).hexdigest()
    (HERE/'results'/'script_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
