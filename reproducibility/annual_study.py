"""Plan-aligned annual synthetic study; all reported values are computed."""
from pathlib import Path
from collections import deque
import json,time,tracemalloc
import numpy as np
from scipy.special import ndtr
from scipy.optimize import minimize_scalar
from analysis import capacities,solve_lp,bus_problem,zero_frontier

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'annual_results';OUT.mkdir(exist_ok=True)
SEEDS=list(range(20261017,20261047))

def ar(rng,n,rho):
    z=np.empty(n);z[0]=rng.normal();innov=rng.normal(size=n)
    for i in range(1,n):z[i]=rho*z[i-1]+np.sqrt(1-rho*rho)*innov[i]
    return z

def generate(seed,n=8760,dt=1.,mode='evening',rho=.95,cross=-.3):
    rng=np.random.default_rng(seed);t=np.arange(n)*dt;day=t/24;hour=t%24
    zw=ar(rng,n,rho**dt);zc=cross*zw+np.sqrt(1-cross**2)*ar(rng,n,.85**dt)
    daylight=12+2*np.sin(2*np.pi*(day-80)/365)
    solar=np.maximum(0,np.sin(np.pi*(hour-(12-daylight/2))/daylight))
    solar*=((hour>=12-daylight/2)&(hour<=12+daylight/2))
    cloud=np.clip(.78+.16*zc,.15,1.)
    pv=np.clip(solar*cloud*(.88+.12*np.sin(2*np.pi*(day-80)/365)),0,1)
    u=np.clip(ndtr(zw),1e-10,1-1e-10)
    speed=8*np.sqrt(-np.log(1-u))*(1+.15*np.cos(2*np.pi*day/365))
    wind=np.where(speed<3,0,np.where(speed<12,(speed**3-3**3)/(12**3-3**3),np.where(speed<=25,1,0)))
    if mode=='flat':base=np.ones(n)
    else:
        peak=19 if mode=='evening' else 7
        base=1+.35*np.exp(-((hour-peak)/2.5)**2)+.12*np.exp(-((hour-12)/4)**2)
    load=base*(1+.08*np.cos(2*np.pi*day/365))*np.where((np.floor(day)%7)>=5,.92,1.)
    load*=np.clip(1+.025*rng.normal(size=n),.8,1.2);load*=100/load.mean()
    return pv,wind,load,speed

def interval_frontier(g,l,pc,pd,ec=.95,ed=.95,dt=1.):
    a,b,d=capacities(g,l,pc,pd,ec,ed,dt);n=len(a)
    reasons=[]
    if a.sum()<b.sum()-1e-7:reasons.append('recharge')
    if np.any(d<b-1e-7):reasons.append('discharge')
    if reasons:return dict(E=float('inf'),reason='+'.join(reasons),start=None,length=None)
    pref=np.r_[0,np.cumsum(np.tile(b-a,2))];queue=deque([0]);best=0.;start=length=0
    for j in range(1,2*n+1):
        while queue and queue[0]<j-n:queue.popleft()
        value=pref[j]-pref[queue[0]]
        if value>best+1e-9:best=float(value);start=queue[0]%n;length=j-queue[0]
        while queue and pref[queue[-1]]>=pref[j]:queue.pop()
        queue.append(j)
    return dict(E=best,reason='finite',start=int(start),length=int(length))

def charge_floor(g,l,ec=.95,ed=.95):
    surplus=np.maximum(g-l,0);target=np.maximum(l-g,0).sum()/(ec*ed)
    if surplus.sum()<target-1e-7:return float('inf')
    lo=0.;hi=float(surplus.max())
    for _ in range(60):
        mid=(lo+hi)/2
        if np.minimum(surplus,mid).sum()>=target:hi=mid
        else:lo=mid
    return hi

def longest_cyclic(mask):
    if np.all(mask):return len(mask)
    best=run=0
    for v in np.tile(mask,2):
        run=run+1 if v else 0;best=max(best,run)
    return min(len(mask),best)

def dispatch_metrics(g,l,pc,pd,eps=0,ec=.95,ed=.95,dt=1.,tag=None,**kw):
    E,res=solve_lp(g,l,pc,pd,eps,ec,ed,dt,**kw)
    if not np.isfinite(E):return dict(E=E,status=res.status)
    n=len(g);e,c,d,u,s=np.split(res.x[:-1],5);u=np.maximum(u,0)
    balance=g+d+u-l-c-s;dynamic=np.roll(e,-1)-e-ec*c*dt+d*dt/ed
    residual=max(np.max(abs(balance)),np.max(abs(dynamic)),max(0.,np.max(e-E)),max(0.,-np.min(e)))
    out=dict(E=E,status=res.status,ENS=float(u.sum()*dt),fraction=float(u.sum()*dt/(sum(l)*dt)),shortage_hours=float(np.sum(u>1e-6)*dt),longest_shortage_hours=float(longest_cyclic(u>1e-6)*dt),curtailment=float(s.sum()*dt),max_residual=float(residual),e0=float(e[0]))
    if tag:np.savetxt(OUT/(tag+'.csv'),np.c_[np.arange(n)*dt,g,l,e,c,d,u,s],delimiter=',',header='hour,g_kw,load_kw,e_kwh,charge_kw,discharge_kw,unserved_kw,curtailed_kw',comments='')
    assert residual<1e-5,out
    return out


def run():
    mix=[];annual_checks=[];pc=200.;pd=160.;K=600.;alphas=np.linspace(0,1,11)
    for j,seed in enumerate(SEEDS):
        p,w,l,v=generate(seed)
        for alpha in alphas:
            g=K*(alpha*p+(1-alpha)*w);f=interval_frontier(g,l,pc,pd)
            gn=g*(1.4*l.sum()/g.sum());fn=interval_frontier(gn,l,pc,pd)
            mix.append(dict(seed=seed,alpha=float(alpha),yield_ratio=float(g.sum()/l.sum()),**f,normalized_E=fn['E']))
        g=K*(.5*p+.5*w);start=time.perf_counter();f=interval_frontier(g,l,pc,pd);scan=time.perf_counter()-start
        start=time.perf_counter();lp=dispatch_metrics(g,l,pc,pd,tag='annual_base_dispatch' if j==0 else None);elapsed=time.perf_counter()-start
        assert (np.isinf(f['E']) and np.isinf(lp['E'])) or abs(f['E']-lp['E'])<1e-5
        annual_checks.append(dict(seed=seed,scan_seconds=scan,lp_seconds=elapsed,**lp,analytical_E=f['E']))
        print('Annual seed',j+1,'/30',flush=True)
    p,w,l,v=generate(SEEDS[0]);g=K*(.5*p+.5*w)
    np.savetxt(OUT/'annual_base_profiles.csv',np.c_[np.arange(8760),p,w,v,l,g],delimiter=',',header='hour,pv_pu,wind_pu,wind_speed_ms,load_kw,generation_kw',comments='')
    # Refine the one-dimensional fixed-nameplate mixture, then verify the selected point.
    def objective(alpha):return interval_frontier(K*(alpha*p+(1-alpha)*w),l,pc,pd)['E']
    opt=minimize_scalar(objective,bounds=(0,1),method='bounded',options={'xatol':1e-6})
    optimum=dict(alpha=float(opt.x),E=float(opt.fun),lp=dispatch_metrics(K*(opt.x*p+(1-opt.x)*w),l,pc,pd))
    print('Mix optimum',optimum,flush=True)
    frontier=[]
    for power in np.linspace(30,300,55):frontier.append(dict(Pc=float(power),**interval_frontier(g,l,power,pd)))
    sensitivity=[]
    for ec in [.85,.9,.95,1.]:
        for ed in [.85,.9,.95,1.]:
            f=interval_frontier(g,l,pc,pd,ec,ed);lp=dispatch_metrics(g,l,pc,pd,ec=ec,ed=ed)
            assert abs(f['E']-lp['E'])<1e-5
            sensitivity.append(dict(ec=ec,ed=ed,**f,charging_floor=charge_floor(g,l,ec,ed),nominal_60=f['E']/.6,nominal_80=f['E']/.8,nominal_90=f['E']/.9))
    reliability=[dict(eps=eps,**dispatch_metrics(g,l,pc,pd,eps=eps,tag=f'dispatch_eps_{eps:.2f}')) for eps in [0,.01,.03,.05]]
    # Explicit finite scenario family: durations, demand patterns, wind persistence,
    # resource correlation and renewable-energy adequacy are varied separately.
    scenarios=[('base',g,l)]
    for hours in [6,12,24,48,72]:
        gg=g.copy()
        for start in [1200,4300,7000]:gg[start:start+hours]*=.1
        scenarios.append((f'joint_low_{hours}h',gg,l))
    for mode in ['morning','flat']:
        pp,ww,ll,_=generate(SEEDS[0],mode=mode);scenarios.append((mode+' demand',K*(pp+ww)/2,ll))
    for factor in [.7,1.3]:scenarios.append((f'wind_scale_{factor}',K*(.5*p+.5*np.clip(w*factor,0,1)),l))
    for rho in [.7,.99]:
        pp,ww,ll,_=generate(SEEDS[0],rho=rho);scenarios.append((f'persistence_{rho}',K*(pp+ww)/2,ll))
    for cross in [-.8,.8]:
        pp,ww,ll,_=generate(SEEDS[0],cross=cross);scenarios.append((f'latent_cross_{cross}',K*(pp+ww)/2,ll))
    scenario_results=[];robust_grid=[]
    for name,gg,ll in scenarios:
        f=interval_frontier(gg,ll,pc,pd);lp=dispatch_metrics(gg,ll,pc,pd)
        assert (np.isinf(f['E']) and np.isinf(lp['E'])) or abs(f['E']-lp['E'])<1e-5
        scenario_results.append(dict(name=name,**f,lp_E=lp['E']))
    for power in [50,75,100,125,150,175,200,250,300]:
        fs=[interval_frontier(gg,ll,power,pd) for _,gg,ll in scenarios];worst=int(np.argmax([f['E'] for f in fs]))
        robust_grid.append(dict(Pc=power,E=fs[worst]['E'],controlling=scenarios[worst][0],start=fs[worst]['start'],length=fs[worst]['length']))
    adequacy=[]
    for factor in [.95,1.05,1.2,1.4]:adequacy.append(dict(factor=factor,**interval_frontier(g*(factor*l.sum()/g.sum()),l,pc,pd)))
    # Explicit omissions relative to the matched model; no claim these are full
    # implementations of another publication.
    ablations=[]
    for name,gg,ll,pp,dd,ec,ed in [('baseline',g,l,pc,pd,.95,.95),('ignore_charge_limit',g,l,1e6,pd,.95,.95),('ignore_discharge_limit',g,l,pc,1e6,.95,.95),('unit_efficiency',g,l,pc,pd,1,1),('pv_only',K*p,l,pc,pd,.95,.95),('wind_only',K*w,l,pc,pd,.95,.95)]:ablations.append(dict(name=name,**interval_frontier(gg,ll,pp,dd,ec,ed)))
    # The plan's six-step example; save a witness dispatch and fixed-start check.
    tiny_l=np.full(6,3.);tiny_g=tiny_l+np.array([4,-3,1,-2,3,-1.])
    tiny=dispatch_metrics(tiny_g,tiny_l,2.5,3,ec=1,ed=1,tag='six_step_dispatch')
    assert abs(tiny['E']-4)<1e-8
    # Annual boundary with 15-minute forcing. Different draws imply a resolution
    # sensitivity, not a controlled discretisation-error estimate.
    pp,ww,ll,_=generate(SEEDS[0],n=35040,dt=.25);gg=K*(pp+ww)/2
    quarter=interval_frontier(gg,ll,pc,pd,dt=.25);qlp=dispatch_metrics(gg,ll,pc,pd,dt=.25)
    assert abs(quarter['E']-qlp['E'])<1e-5
    # Greedy cyclic dispatch feasibility for fixed-day autonomy heuristics is
    # checked using the exact zero-shortfall boundary.
    baseline=interval_frontier(g,l,pc,pd)
    autonomy=[dict(days=d,E=100*24*d/.95,feasible=100*24*d/.95>=baseline['E']) for d in [.5,1,2,3]]
    report=dict(config=dict(seeds=SEEDS,K=K,Pc=pc,Pd=pd,eta_c=.95,eta_d=.95,SOC_window=.8,horizon=8760,dt=1),mix=mix,annual_checks=annual_checks,optimum=optimum,frontier=frontier,sensitivity=sensitivity,reliability=reliability,scenarios=scenario_results,robust_grid=robust_grid,adequacy=adequacy,ablations=ablations,tiny=tiny,quarter_hour=dict(**quarter,lp=qlp),autonomy=autonomy,charging_floor=charge_floor(g,l),discharge_floor=float(np.max(np.maximum(l-g,0))))
    (OUT/'annual_results.json').write_text(json.dumps(report,indent=2))
    print('ANNUAL STUDY SAVED',flush=True)

if __name__=='__main__':run()
