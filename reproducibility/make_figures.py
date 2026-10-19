"""Regenerate manuscript vector figures from the saved experiment results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
data=json.loads((ROOT/'reproducibility/results/results.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight')
    plt.close(fig)

fig,ax=plt.subplots(figsize=(7.2,3.1));ax.set_xlim(-.7,6.7);ax.set_ylim(-2.15,1.3);ax.axis('off')
for j,label in enumerate(['0','1',r'$\cdots$',r'$T-1$']):
    x=2*j
    ax.add_patch(Circle((x,0),.23,facecolor='#e7eef5',edgecolor='#234d6b'))
    ax.text(x,0,label,ha='center',va='center')
    for start,end in [((x,1),(x,.25)),((x,-.25),(x,-1))]:
        ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=12,color='#234d6b'))
    ax.text(x,1.13,r'$s$',ha='center');ax.text(x,-1.15,r'$z$',ha='center')
    ax.text(x+.12,.62,r'$a_t$',ha='left');ax.text(x+.12,-.65,r'$d_t$',ha='left')
for j in range(3):
    ax.add_patch(FancyArrowPatch((2*j+.25,0),(2*j+1.75,0),arrowstyle='-|>',mutation_scale=12,color='#234d6b'))
    ax.text(2*j+1,.12,r'$E$',ha='center')
ax.add_patch(FancyArrowPatch((6,-.25),(0,-.25),connectionstyle='arc3,rad=-.5',arrowstyle='-|>',mutation_scale=12,color='#b05b32'))
ax.text(3,-1.92,r'Cyclic carry-over: capacity $E$',ha='center',color='#934820',bbox={'facecolor':'white','edgecolor':'none','pad':2})
save(fig,'cyclic_network')

fig,axs=plt.subplots(1,3,figsize=(8.3,2.8),layout='constrained')
cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#d8d8d8')
for ax,eps in zip(axs,[0,.01,.05]):
    rows=[r for r in data['grid'] if r['eps']==eps]
    pc=sorted(set(r['Pc'] for r in rows));pd=sorted(set(r['Pd'] for r in rows))
    Z=np.full((len(pd),len(pc)),np.nan)
    for r in rows:
        if np.isfinite(r['E']):Z[pd.index(r['Pd']),pc.index(r['Pc'])]=r['E']
    im=ax.pcolormesh(pc,pd,Z,shading='nearest',cmap=cmap,vmin=0,vmax=12)
    ax.set(xlabel='Charge rating (kW)',title=rf'$\varepsilon={eps:g}$')
axs[0].set_ylabel('Discharge rating (kW)')
fig.colorbar(im,ax=axs,label='Minimum usable energy (kWh)',shrink=.85)
save(fig,'capacity_frontiers')

fig,axs=plt.subplots(1,2,figsize=(7.5,2.8),layout='constrained')
for ax,row in zip(axs,data['timing']):
    g=np.array(row['pv'])+row['wind'];t=np.arange(9)
    ax.stairs(g,t,color='#234d6b',lw=2,label='Generation')
    ax.axhline(2,color='#b05b32',ls='--',label='Demand')
    ax.set(xlabel='Interval boundary (h)',ylabel='Power (kW)',ylim=(-.2,4.6),xticks=range(9),title=f"{row['name'].capitalize()}: $E_{{min}}={row['E']:g}$ kWh")
axs[0].legend(loc='upper right',frameon=False,fontsize=8)
save(fig,'chronology')

fig,ax=plt.subplots(figsize=(5.8,2.9),layout='constrained')
E=np.linspace(0,2.2,221);u=np.maximum(4-2*E,0)
ax.plot(E,u,color='#234d6b',label='Minimum unmet energy')
ax.axhline(1,color='#b05b32',ls='--',label='Shared allowance (1 kWh)')
ax.plot([1.5],[1],marker='o',color='#b05b32');ax.axvline(1.5,color='#b05b32',ls=':',lw=1)
ax.set(xlabel='Usable capacity (kWh)',ylabel='Unmet energy (kWh)',xlim=(0,2.2),ylim=(0,4.2))
ax.legend(frameon=False)
save(fig,'shared_shortfall')
print('Wrote four vector figures and PNG previews.')
