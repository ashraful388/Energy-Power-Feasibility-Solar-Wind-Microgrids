"""Generate statistical summaries, figures and results prose from executed runs."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from annual_study import OUT
R=OUT.parents[1];F=R/'figures'
a=json.loads((OUT/'annual_results.json').read_text());b=json.loads((OUT/'extended_checks.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
rng=np.random.default_rng(17881)
def stats(x):
 x=np.asarray(x);boot=np.median(rng.choice(x,(2000,len(x)),replace=True),axis=1)
 return dict(median=float(np.median(x)),q25=float(np.quantile(x,.25)),q75=float(np.quantile(x,.75)),ci=np.quantile(boot,[.025,.975]).tolist())
mixstats=[]
for alpha in sorted(set(r['alpha'] for r in a['mix'])):
 rows=[r for r in a['mix'] if r['alpha']==alpha]
 mixstats.append(dict(alpha=alpha,finite=int(sum(np.isfinite(r['E']) for r in rows)),E=stats([r['E'] for r in rows]),normalized=stats([r['normalized_E'] for r in rows]),yield_ratio=stats([r['yield_ratio'] for r in rows])))
paired=[]
for seed in a['config']['seeds']:
 rows=[r for r in a['mix'] if r['seed']==seed];e0=rows[0]['E'];em=next(r['E'] for r in rows if abs(r['alpha']-.7)<1e-8);paired.append(100*(e0-em)/e0)
timing={name:stats(values) for name,values in [('scan',[r['scan_seconds'] for r in a['annual_checks']]),('matrix',[r['seconds'] for r in b['matrix']]),('LP',[r['lp_seconds'] for r in a['annual_checks']])]}
summary=dict(mix=mixstats,paired_reduction_percent=stats(paired),timing=timing,max_annual_error=max(abs(r['E']-r['analytical_E']) for r in a['annual_checks']),max_matrix_error=max(abs(x['E']-y['analytical_E']) for x,y in zip(b['matrix'],a['annual_checks'])),max_dispatch_residual=max(r['max_residual'] for r in a['annual_checks']))
(OUT/'statistics.json').write_text(json.dumps(summary,indent=2))
def save(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(8,3.0),layout='constrained')
for ax,key,title in zip(axs,['E','normalized'],['Fixed installed capacity: 600 kW','Equal generated energy: 1.4 × demand']):
 x=[r['alpha'] for r in mixstats];med=[r[key]['median']/1000 for r in mixstats];lo=[r[key]['q25']/1000 for r in mixstats];hi=[r[key]['q75']/1000 for r in mixstats]
 ax.fill_between(x,lo,hi,alpha=.18,color='#234d6b',label='Interquartile range');ax.plot(x,med,'o-',color='#234d6b',ms=3,label='Median (30 seeds)');ax.set(xlabel='PV share',ylabel='Usable energy (MWh)',title=title,xlim=(0,1))
axs[0].legend(frameon=False,fontsize=8)
save(fig,'annual_mix')

fig,axs=plt.subplots(1,3,figsize=(9,3),layout='constrained')
pieces=b['breakpoints']['pieces'];x=np.linspace(b['breakpoints']['charging_floor'],9,300);y=np.maximum(15-4*x,5)
axs[0].plot(x,y,color='#234d6b');axs[0].axvline(2.5,color='#b05b32',ls='--');axs[0].set(xlabel='Charge rating (kW)',ylabel='Usable energy (kWh)',title='Exact 12-hour breakpoint')
rows=[r for r in a['frontier'] if np.isfinite(r['E'])];axs[1].plot([r['Pc'] for r in rows],[r['E']/1000 for r in rows],label='Base trace',color='#234d6b')
rows=[r for r in a['robust_grid'] if np.isfinite(r['E'])];axs[1].plot([r['Pc'] for r in rows],[r['E']/1000 for r in rows],'o-',label='14-scenario bound',color='#b05b32',ms=3);axs[1].axvline(a['charging_floor'],ls=':',color='#555');axs[1].set(xlabel='Charge rating (kW)',ylabel='Usable energy (MWh)',title='Annual capacity boundary');axs[1].legend(frameon=False,fontsize=7)
rows=a['reliability'];axs[2].plot([100*r['eps'] for r in rows],[r['E']/1000 for r in rows],'o-',color='#234d6b');axs[2].set(xlabel='Unserved-energy allowance (%)',ylabel='Usable energy (MWh)',title='Annual reliability relaxation')
save(fig,'annual_frontiers')

prof=np.loadtxt(OUT/'annual_base_profiles.csv',delimiter=',',skiprows=1);fig,axs=plt.subplots(3,1,figsize=(8,4.5),layout='constrained')
daily=prof[:,[1,2,4,5]].reshape(365,24,4).mean(axis=1)
axs[0].plot(daily[:,0]*300,label='PV');axs[0].plot(daily[:,1]*300,label='Wind');axs[0].plot(daily[:,2],color='black',lw=.8,label='Demand');axs[0].set(ylabel='Daily mean (kW)',xlabel='Day');axs[0].legend(ncol=3,frameon=False,fontsize=8)
i=slice(6858,7129);axs[1].plot(prof[i,0],prof[i,5],lw=.7,label='Generation');axs[1].plot(prof[i,0],prof[i,4],lw=.7,label='Demand');axs[1].set(xlabel='Hour',ylabel='Power (kW)',title='Base-case critical interval');axs[1].legend(frameon=False,fontsize=8)
dispatch=np.loadtxt(OUT/'annual_base_dispatch.csv',delimiter=',',skiprows=1);axs[2].plot(dispatch[i,0],dispatch[i,3]/1000,color='#234d6b');axs[2].set(xlabel='Hour',ylabel='Inventory (MWh)')
save(fig,'annual_profiles')

base=a['annual_checks'][0]['E'];best=a['optimum'];m7=mixstats[7];normbest=min(mixstats,key=lambda row:row['normalized']['median']);reduction=summary['paired_reduction_percent'];err=summary['max_annual_error']
text=r'''\subsection{Annual Verification, Critical Intervals, and Renewable Mix}
\label{subsec:annual-results}

The annual study retained every generated case, including the deliberately recharge-infeasible adequacy test. All 30 base-mix annual boundaries agreed with the independent linear programme, with a maximum absolute discrepancy of ERR~kWh. The restricted Ren difference-matrix calculation also agreed in all 30 cases. All 15 fixed-capacity checks correctly classified designs immediately below, on, and above the boundary. The 15-minute annual case gave 7238.117~kWh by both the interval scan and linear programming. Its capacity difference from the hourly case reflects different random forcing as well as resolution and is not attributed solely to discretisation.

For the base seed at $\alpha=0.5$, the capacity is BASE~kWh, the charging-power floor is 58.815~kW, and the discharge-power floor is 140.931~kW. The limiting interval begins at hour 6858 and spans 271 hours (Figure~\ref{fig:annual-profiles}). Its duration exceeds an individual low-generation episode because partial recharge between deficits does not necessarily restore the previous inventory. The 80\% usable state-of-charge window corresponds to NAMEPLATE~kWh nameplate capacity. A three-day autonomy rule based on mean demand and discharge efficiency assigns 7578.947~kWh and is infeasible for this trace; the half-, one-, and two-day rules are also insufficient.

\begin{figure}[H]\centering
\includegraphics[width=\textwidth]{figures/annual_profiles.pdf}
\caption{Synthetic annual daily means and the base-case limiting interval, with the corresponding inventory from the independent optimal dispatch. Daily means are displayed only for visualisation; sizing retains hourly chronology.}\label{fig:annual-profiles}
\end{figure}

At fixed $K=600$~kW, the refined base-seed minimum occurs at $\alpha=ALPHA$, with BEST~kWh. Pure wind and pure photovoltaic designs require 39101.515 and 75813.880~kWh, respectively. The optimum is specific to the generator and seed; it is not a recommended equipment mix for a site. Across 30 seeds, the lowest median on the prespecified 0.1-spaced grid occurs at $\alpha=0.7$, with median MEDIAN~kWh, interquartile range QLOW--QHIGH~kWh, and a 95\% bootstrap interval of CILOW--CIHIGH~kWh for the median. The paired reduction relative to wind-only sizing has median REDUCTION\%, with bootstrap interval RLOW--RHIGH\%. All 330 fixed-capacity mix cases were finite. Under the separate equal-generated-energy control, the lowest grid median occurs at $\alpha=NORMALPHA$, demonstrating that the installed-capacity and generated-energy controls answer different questions (Figure~\ref{fig:annual-mix}).

\begin{figure}[H]\centering
\includegraphics[width=\textwidth]{figures/annual_mix.pdf}
\caption{Resource-mix comparisons over 30 seeds. The left panel preserves installed renewable capacity; the right rescales generated energy and therefore does not preserve installed capacity. Shading gives the interquartile range, not a real-weather confidence band.}\label{fig:annual-mix}
\end{figure}

For an exact breakpoint illustration, the 12-hour residual sequence, in kW, is
\[(-2,-1,-2,8,3,-5,2,9,9,-5,4,8).\]
Efficiencies are unity and $P_{\mathrm d}=6$~kW. The charge floor is $13/6$~kW. Above that floor,
\begin{equation}
 E_{\min}(P_{\mathrm c})=\max\{15-4P_{\mathrm c},5\}\quad\mathrm{kWh}.
\end{equation}
The breakpoint at 2.5~kW is not an individual surplus value. A nine-interval wrapping certificate controls below it; a three-interval deficit controls above it. Complete interval enumeration and linear programming verify both branches. Additional annual convexity checks found no violation in 86 randomly selected feasible pairs; 14 sampled pairs with an infeasible endpoint were not used in this finite-value inequality check. The convexity guarantee itself follows from Proposition~\ref{prop:mix-convex}, not from these numerical checks.

\begin{figure}[H]\centering
\includegraphics[width=\textwidth]{figures/annual_frontiers.pdf}
\caption{Exact small-instance breakpoint, annual base and finite-scenario capacity requirements, and annual reliability relaxation. The dotted line in the middle panel marks the base-case charging floor; scenario feasibility can impose a larger floor.}\label{fig:annual-frontiers}
\end{figure}

\subsection{Annual Reliability, Nonidealities, and Ablations}

Table~\ref{tab:annual-reliability} reports the capacity and shortage statistics for the annual base trace. Allowing a small energy shortage substantially changes the capacity requirement, but it does not bound outage duration. Reported shortage durations are those of the solver's returned capacity-optimal dispatch and are not claimed to be duration-optimal. Across the 30 zero-shortfall annual dispatches, power-balance residuals were below RESIDUAL~kW and state-equation residuals were below RESIDUAL~kWh; no violation exceeded the stated tolerance.

\begin{table}[H]
\caption{Annual base-trace reliability sensitivity at fixed power ratings.}\label{tab:annual-reliability}
\begin{tabularx}{\textwidth}{Xrrrr}
\toprule Allowance (\%) & $E$ (kWh) & ENS (kWh) & Hours & Max. run (h)\\\midrule
RELIABILITYROWS
\bottomrule\end{tabularx}\end{table}

At unit charging and discharging efficiencies, the annual requirement is 7172.595~kWh, compared with 8146.632~kWh at 0.95 each. All 16 independently varied efficiency pairs agree with the reference linear programme. Changing the usable state-of-charge window changes nameplate capacity by the reciprocal window width, while leaving the usable-energy theorem unchanged. Removing the charging cap at the baseline 200~kW has no material capacity effect. At 75~kW, however, the correct capacity is 11760.873~kWh; omitting that cap underestimates it as 8146.632~kWh. At 40~kW, omission changes an infeasible design into a falsely finite one. Likewise, ignoring a 100~kW discharge rating gives a falsely finite design because the actual requirement reaches 140.931~kW.

Removing terminal cyclicity did not change the optimum in this particular annual base case; no universal effect is inferred. Imposing empty initial storage instead made the trace infeasible. Hourly retention factors of 0.9999 and 0.999 increased the LP-only capacity to 8275.220 and 9553.278~kWh. These results quantify sensitivity to a changed model and are not validations of the no-leakage cut formula under leakage. Rescaling annual generation to 0.95 times demand caused recharge infeasibility. Ratios of 1.05, 1.2, and 1.4 produced capacities of 31962.500, 15438.536, and 10070.352~kWh; annual energy surplus alone therefore did not determine the required inventory.

\subsection{Annual Scenario and Published-Method Comparisons}

Every one of the 14 declared annual scenarios was feasible at the baseline power ratings, with matching interval and LP capacities. Repeated 72-hour joint low-generation events increased the requirement to 12847.189~kWh. The most persistent wind forcing, with autoregressive coefficient 0.99, controlled the shared design at 16532.318~kWh through a 567-hour interval beginning at hour 1025. This exceeds the base-case requirement but supplies no guarantee beyond the specified scenarios or for a shared causal dispatch. The controlling scenario and interval for each tested charge rating are supplied in the machine-readable results.

The restricted Ren cumulative matrix, the interval calculation, and the reference LP returned the same zero-shortfall annual capacity within tolerance. Paired median times across 30 seeds were SCANTIME~ms for the interval scan, MATRIXTIME~ms for the difference matrix, and LPTIME~ms for the LP including dispatch extraction. These measurements support a speed difference between these implementations on this generator, not a claim that the cumulative principle is new or that the full published method has been outperformed. Peak Python-traced allocations are provided for four horizons; they do not measure total native-solver resident memory. Capacity differences between the implementations are numerical rounding, rather than meaningful optimality gaps.
'''
def texsci(value):
 m,e=f'{value:.2e}'.split('e');return '$'+f'{float(m):g}'+r'\times10^{'+str(int(e))+'}$'
repl={'ERR':texsci(err),'BASE':f'{base:.3f}','NAMEPLATE':f'{base/.8:.3f}','ALPHA':f"{best['alpha']:.6f}",'BEST':f"{best['E']:.3f}",'MEDIAN':f"{m7['E']['median']:.3f}",'QLOW':f"{m7['E']['q25']:.3f}",'QHIGH':f"{m7['E']['q75']:.3f}",'CILOW':f"{m7['E']['ci'][0]:.3f}",'CIHIGH':f"{m7['E']['ci'][1]:.3f}",'REDUCTION':f"{reduction['median']:.2f}",'RLOW':f"{reduction['ci'][0]:.2f}",'RHIGH':f"{reduction['ci'][1]:.2f}",'NORMALPHA':f"{normbest['alpha']:.1f}",'RESIDUAL':texsci(summary['max_dispatch_residual']),'SCANTIME':f"{1000*timing['scan']['median']:.2f}",'MATRIXTIME':f"{1000*timing['matrix']['median']:.2f}",'LPTIME':f"{1000*timing['LP']['median']:.2f}"}
for key in sorted(repl,key=len,reverse=True):text=text.replace(key,repl[key])
rows=[' & '.join([f"{100*r['eps']:.0f}",f"{r['E']:.3f}",f"{r['ENS']:.0f}",f"{r['shortage_hours']:.0f}",f"{r['longest_shortage_hours']:.0f}"])+r'\\' for r in a['reliability']]
text=text.replace('RELIABILITYROWS','\n'.join(rows))
(R/'annual_results.tex').write_text(text,encoding='utf-8')
print(json.dumps(summary,indent=2))
