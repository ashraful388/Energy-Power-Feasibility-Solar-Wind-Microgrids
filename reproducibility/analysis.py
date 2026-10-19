"""Reproducible synthetic study. Run: python analysis.py. No measured data."""
from pathlib import Path
import itertools, json, time, platform, sys, csv
import numpy as np
import scipy
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, vstack, coo_matrix
import networkx as nx

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
RESULTS.mkdir(exist_ok=True)

def capacities(g,l,pc,pd,ec=1.,ed=1.,dt=1.):
    g=np.asarray(g,float);l=np.asarray(l,float)
    a=ec*np.minimum(pc,np.maximum(g-l,0))*dt
    b=np.maximum(l-g,0)*dt/ed
    d=np.minimum(pd,np.maximum(l-g,0))*dt/ed
    return a,b,d

def cut(a,d,E):
    n=len(a);G=nx.DiGraph();G.add_nodes_from(range(n+2));s=n;z=n+1
    for t in range(n):
        G.add_edge(s,t,capacity=float(a[t]));G.add_edge(t,z,capacity=float(d[t]))
        G.add_edge(t,(t+1)%n,capacity=float(E))
    value,(left,right)=nx.minimum_cut(G,s,z,flow_func=nx.algorithms.flow.preflow_push)
    S=set(left)-{s};k=sum(t in S and (t+1)%n not in S for t in range(n))
    C=sum(a[t] for t in range(n) if t not in S)+sum(d[t] for t in S)
    return float(value),float(C),int(k),sorted(S)

def frontier(g,l,pc,pd,eps=0.,ec=1.,ed=1.,dt=1.):
    a,b,d=capacities(g,l,pc,pd,ec,ed,dt);R=max(0.,sum(b)-eps*sum(l)*dt/ed)
    if R<=1e-10:return 0.,[]
    if min(sum(a),sum(d))<R-1e-9:return float('inf'),[]
    E=0.;history=[]
    for _ in range(len(a)//2+2):
        F,C,k,S=cut(a,d,E);history.append(dict(E=E,F=F,C=C,k=k,S=S))
        if F>=R-1e-8:return E,history
        if k==0:return float('inf'),history
        E=max(E,(R-C)/k)
    raise RuntimeError('cut-update bound exceeded')

def enumeration(g,l,pc,pd,eps=0.,ec=1.,ed=1.,dt=1.):
    a,b,d=capacities(g,l,pc,pd,ec,ed,dt);n=len(a);R=max(0.,sum(b)-eps*sum(l)*dt/ed)
    if min(sum(a),sum(d))<R-1e-9:return float('inf')
    E=0.
    for mask in range(1,2**n-1):
        S=np.array([(mask>>i)&1 for i in range(n)],bool)
        k=np.sum(S & ~np.roll(S,-1));C=a[~S].sum()+d[S].sum()
        E=max(E,(R-C)/k)
    return float(E)

def bus_problem(g,l,pc,pd,eps=0.,ec=1.,ed=1.,dt=1.,Efixed=None,restricted=True):
    """Independent bus-power LP: e,c,d,u,s,E, with cyclic state indices."""
    g=np.asarray(g,float);l=np.asarray(l,float);n=len(g);N=5*n+1
    eq=lil_matrix((2*n,N));rhs=np.r_[l-g,np.zeros(n)]
    for t in range(n):
        eq[t,n+t]=-1;eq[t,2*n+t]=1;eq[t,3*n+t]=1;eq[t,4*n+t]=-1
        eq[n+t,(t+1)%n]=1;eq[n+t,t]-=1
        eq[n+t,n+t]=-ec*dt;eq[n+t,2*n+t]=dt/ed
    ub=lil_matrix((n+1,N));bound=np.zeros(n+1)
    for t in range(n):ub[t,t]=1;ub[t,-1]=-1
    ub[n,3*n:4*n]=dt;bound[n]=eps*sum(l)*dt
    cb=np.minimum(pc,np.maximum(g-l,0)) if restricted else np.full(n,pc)
    db=np.minimum(pd,np.maximum(l-g,0)) if restricted else np.full(n,pd)
    bounds=[(0,None)]*n+[(0,float(v)) for v in cb]+[(0,float(v)) for v in db]+[(0,float(v)) for v in l]+[(0,float(v)) for v in g]+[(0,None) if Efixed is None else (Efixed,Efixed)]
    objective=np.zeros(N);objective[-1]=1
    return objective,ub.tocsr(),bound,eq.tocsr(),rhs,bounds

def solve_lp(*args,**kwargs):
    c,A,b,Q,r,bounds=bus_problem(*args,**kwargs)
    res=linprog(c,A_ub=A,b_ub=b,A_eq=Q,b_eq=r,bounds=bounds,method='highs')
    if res.status==2:return float('inf'),res
    if not res.success:raise RuntimeError(res.message)
    return float(res.x[-1]),res

def zero_frontier(g,l,pc,pd,ec=1.,ed=1.,dt=1.):
    a,b,d=capacities(g,l,pc,pd,ec,ed,dt)
    if np.any(d<b-1e-9) or a.sum()<b.sum()-1e-9:return float('inf')
    # Two traversals of sequent-peak recurrence capture wrap-around.
    run=peak=0.
    for x in np.tile(b-a,2):run=max(0.,run+x);peak=max(peak,run)
    return peak

def same(x,y):return (np.isinf(x) and np.isinf(y)) or abs(x-y)<2e-7

def profiles(n=336,seed=20261017):
    rng=np.random.default_rng(seed);t=np.arange(n);h=t%24
    pv=np.maximum(0,np.sin(np.pi*(h-6)/12))
    cloud=np.clip(.82+.13*np.sin(2*np.pi*t/120)+.05*rng.normal(size=n),.3,1.)
    wind=np.clip(.47+.22*np.sin(2*np.pi*t/83+.8)+.11*np.cos(2*np.pi*t/37),.05,.9)
    load=1+.3*np.exp(-((h-19)/3)**2)+.12*np.exp(-((h-7)/2)**2)
    g=2.8*pv*cloud+1.7*wind
    return g,load,pv*cloud,wind

def robust_joint(scenarios,pc,pd,eps,ec,ed):
    probs=[bus_problem(g,l,pc,pd,eps,ec,ed) for g,l in scenarios]
    sizes=[len(p[0])-1 for p in probs];N=sum(sizes)+1;U=[];Q=[];bs=[];rs=[];bounds=[];offset=0
    for size,(_,A,b,B,r,bo) in zip(sizes,probs):
        for matrix,store in [(A,U),(B,Q)]:
            co=matrix.tocoo();cols=np.where(co.col==size,N-1,co.col+offset)
            store.append(coo_matrix((co.data,(co.row,cols)),shape=(matrix.shape[0],N)).tocsr())
        bs.extend(b);rs.extend(r);bounds.extend(bo[:-1]);offset+=size
    bounds.append((0,None));c=np.zeros(N);c[-1]=1
    res=linprog(c,A_ub=vstack(U),b_ub=bs,A_eq=vstack(Q),b_eq=rs,bounds=bounds,method='highs')
    if res.status==2:return float('inf')
    if not res.success:raise RuntimeError(res.message)
    return float(res.x[-1])

def run():
    rng=np.random.default_rng(20261017);checks=[];errors=[]
    for q in itertools.product([-1.,0.,1.],repeat=4):
        l=np.ones(4);g=l+q
        for pc,pd,eps in itertools.product([.5,1.],[.5,1.],[0.,.1,.25]):
            ex=enumeration(g,l,pc,pd,eps);cf,h=frontier(g,l,pc,pd,eps);lp,_=solve_lp(g,l,pc,pd,eps)
            assert same(ex,cf) and same(ex,lp),(q,pc,pd,eps,ex,cf,lp)
            checks.append(dict(kind='exhaustive',finite=bool(np.isfinite(cf)),error=abs(cf-lp) if np.isfinite(cf) else 0.,calls=len(h)))
    print('Exhaustive checks',len(checks),flush=True)
    for j in range(240):
        n=int(rng.integers(3,13));l=rng.uniform(.3,2.5,n);g=rng.uniform(0,4,n)
        pc,pd=rng.uniform(.2,3,2);ec,ed=rng.uniform(.8,1.,2);eps=float(rng.choice([0.,.02,.1,.3]))
        cf,h=frontier(g,l,pc,pd,eps,ec,ed);ex=enumeration(g,l,pc,pd,eps,ec,ed);lp,res=solve_lp(g,l,pc,pd,eps,ec,ed)
        broad,_=solve_lp(g,l,pc,pd,eps,ec,ed,restricted=False)
        assert same(cf,ex) and same(cf,lp) and same(cf,broad),(j,cf,ex,lp,broad)
        ks=[x['k'] for x in h[:-1]];assert all(a>b for a,b in zip(ks,ks[1:])),ks
        if eps==0:assert same(cf,zero_frontier(g,l,pc,pd,ec,ed))
        checks.append(dict(kind='random',finite=bool(np.isfinite(cf)),error=abs(cf-lp) if np.isfinite(cf) else 0.,calls=len(h)))
    print('Random checks completed',flush=True)
    examples=[]
    cases=[('wrap-around',[1,3,3,1],[2]*4,1.,1.,0.),('charge bottleneck',[5,0],[1,1],.5,2.,0.),('discharge bottleneck',[0,4],[2,1],4.,1.,0.),('two deficit intervals',[0,4,0,4],[2]*4,2.,2.,.125)]
    for name,g,l,pc,pd,eps in cases:
        E,h=frontier(g,l,pc,pd,eps);examples.append(dict(name=name,g=g,l=l,Pc=pc,Pd=pd,eps=eps,E=E,history=h));assert same(E,solve_lp(g,l,pc,pd,eps)[0])
    pv=np.array([0,0,2,2,0,0,2,2.]);wind=np.array([0,0,0,0,2,2,2,2.]);perm=np.array([0,6,2,4,1,7,3,5]);L=np.full(8,2.)
    timing=[]
    for name,p,w in [('clustered',pv,wind),('interleaved',pv[perm],wind[perm])]:
        E,h=frontier(p+w,L,2,2);timing.append(dict(name=name,pv=p.tolist(),wind=w.tolist(),load=L.tolist(),E=E,correlation=float(np.corrcoef(p,w)[0,1])))
    g,l,p,w=profiles();np.savetxt(RESULTS/'synthetic_profiles.csv',np.c_[np.arange(len(g)),p,w,g,l],delimiter=',',header='hour,pv_per_unit,wind_per_unit,generation_kw,load_kw',comments='')
    grid=[]
    for pc,pd,eps in itertools.product(np.linspace(.4,2.2,10),np.linspace(.5,1.6,8),[0.,.01,.05]):
        start=time.perf_counter();E,h=frontier(g,l,pc,pd,eps,.95,.95);duration=time.perf_counter()-start
        stlp=time.perf_counter();ep,_=solve_lp(g,l,pc,pd,eps,.95,.95);lp_duration=time.perf_counter()-stlp;assert same(E,ep),(pc,pd,eps,E,ep)
        grid.append(dict(Pc=pc,Pd=pd,eps=eps,E=E,lp_E=ep,calls=len(h),seconds=duration,lp_seconds=lp_duration))
    print('Frontier grid completed',flush=True)
    scenarios=[(g,l),(np.roll(g,6),l),(g*.9,l*1.05),(g*np.where((np.arange(len(g))>=120)&(np.arange(len(g))<156),.4,1.),l)]
    er=[frontier(gs,ls,2,1.6,.02,.95,.95)[0] for gs,ls in scenarios];joint=robust_joint(scenarios,2,1.6,.02,.95,.95);assert same(max(er),joint)
    timing_runs=[]
    for n in [168,720,2160,8760]:
        gg,ll,_,_=profiles(n);f=[];z=[];lp=[]
        for rep in range(3):
            st=time.perf_counter();ez=zero_frontier(gg,ll,2,1.6,.95,.95);z.append(time.perf_counter()-st)
            st=time.perf_counter();el,_=solve_lp(gg,ll,2,1.6,0,.95,.95);lp.append(time.perf_counter()-st);assert same(ez,el)
        timing_runs.append(dict(T=n,E=ez,scan_seconds=float(np.median(z)),lp_seconds=float(np.median(lp))))
        print('Runtime',n,timing_runs[-1],flush=True)
    summary=dict(seed=20261017,versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,networkx=nx.__version__,platform=platform.platform(),processor=platform.processor()),checks=checks,examples=examples,timing=timing,grid=grid,scenario_E=er,robust_E=joint,runtimes=timing_runs)
    (RESULTS/'results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print('RESULTS',json.dumps(dict(checks=len(checks),finite=sum(c['finite'] for c in checks),max_error=max(c['error'] for c in checks),max_calls=max(c['calls'] for c in checks),grid_finite=sum(bool(np.isfinite(c['E'])) for c in grid),robust_E=joint,scenario_E=er)),flush=True)

if __name__=='__main__':run()
