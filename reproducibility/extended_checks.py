"""Published cumulative baseline, frontier certificates and sensitivity audits."""
import json,time,tracemalloc
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from annual_study import generate,interval_frontier,charge_floor,SEEDS,OUT
from analysis import capacities,bus_problem,solve_lp

def ren_matrix(g,l,pc,pd,ec=.95,ed=.95,dt=1.):
    """Restricted zero-leakage Eq. (5)-(6) difference-matrix reimplementation.
    Uses the same power-clipped increment as the analytical/LP benchmark.
    Does not reproduce Ren's leakage or capacity-dependent power iteration.
    """
    a,b,d=capacities(g,l,pc,pd,ec,ed,dt)
    if a.sum()<b.sum()-1e-7 or np.any(d<b-1e-7):return float('inf'),0
    inc=a-b;cum=np.r_[0,np.cumsum(inc)];nz=np.flatnonzero(abs(inc)>1e-12)
    turn=nz[1:][np.sign(inc[nz[1:]])!=np.sign(inc[nz[:-1]])]
    levels=cum[np.unique(np.r_[0,turn,len(inc)])];m=len(levels)
    D=levels[None,:]-levels[:,None]
    D+=np.tril(np.full((m,m),cum[-1]),-1)
    return max(0.,float(-D.min())),m

def benchmark():
    rows=[]
    for seed in SEEDS:
        p,w,l,v=generate(seed);g=300*(p+w)
        st=time.perf_counter();E,m=ren_matrix(g,l,200,160);dur=time.perf_counter()-st
        target=interval_frontier(g,l,200,160)['E'];assert abs(E-target)<1e-5
        rows.append(dict(seed=seed,E=E,critical_points=m,seconds=dur))
    memory=[]
    from analysis import zero_frontier
    for n in [168,720,2160,8760]:
        p,w,l,v=generate(SEEDS[0],n=n);g=300*(p+w)
        for name,fn in [('cumulative_scan',lambda:zero_frontier(g,l,200,160,.95,.95)),('critical_matrix',lambda:ren_matrix(g,l,200,160)),('LP',lambda:solve_lp(g,l,200,160,ec=.95,ed=.95))]:
            tracemalloc.start();st=time.perf_counter();fn();duration=time.perf_counter()-st;_,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
            memory.append(dict(T=n,method=name,seconds_with_tracing=duration,python_peak_bytes=peak))
    return rows,memory

def modified_lp(g,l,pc,pd,ec=.95,ed=.95,initial=None,cyclic=True,retention=1.):
    c,A,b,Q,r,bo=bus_problem(g,l,pc,pd,ec=ec,ed=ed)
    n=len(g);Q=Q.tolil()
    for t in range(n):Q[n+t,t]=-retention
    if not cyclic:
        # Retain a genuine terminal inventory via one added variable, rather
        # than deleting the last dynamics equation.
        from scipy.sparse import hstack,csr_matrix
        Q=hstack([Q.tocsr(),csr_matrix((2*n,1))]).tolil();Q[2*n-1,0]=0;Q[2*n-1,-1]=1
        A=hstack([A,csr_matrix((n+1,1))]).tolil()
        from scipy.sparse import vstack
        row=np.zeros(n*5+2);row[-1]=1;row[-2]=-1
        A=vstack([A.tocsr(),row]).tocsr();b=np.r_[b,0];c=np.r_[c,0];bo=bo+[(0,None)]
    if initial is not None:bo[0]=(initial,initial)
    res=linprog(c,A_ub=A,b_ub=b,A_eq=Q.tocsr(),b_eq=r,bounds=bo,method='highs')
    if res.status==2:return dict(E=float('inf'),status=res.status)
    if not res.success:raise RuntimeError(res.message)
    return dict(E=float(res.x[5*n]),status=res.status)

def exact_small_breakpoints():
    # Enumerate all circular intervals and their affine expressions between
    # clipping thresholds; then retain actual upper-envelope changes.
    rng=np.random.default_rng(741)
    for trial in range(1000):
        q=rng.integers(-5,12,12).astype(float);l=np.full(12,6.);g=l+q
        if q.sum()<=0:continue
        floor=charge_floor(g,l,1,1);hi=max(q)
        if not np.isfinite(floor) or hi-floor<.2:continue
        clips=sorted(set([floor,hi]+[float(x) for x in q if floor<x<hi]));pieces=[]
        intervals=[[(s+i)%12 for i in range(k)] for s in range(12) for k in range(1,12)]
        for lo,up in zip(clips,clips[1:]):
            mid=(lo+up)/2;lines=[(0.,0.,None)]
            for ids in intervals:
                v=q[ids];m=-float(np.sum(v>mid));b=float(np.maximum(-v,0).sum()-v[(v>0)&(v<mid)].sum());lines.append((m,b,ids))
            pts=[lo,up]
            for i,(m,b,ids) in enumerate(lines):
                for mm,bb,_ in lines[i+1:]:
                    if m!=mm:
                        x=(bb-b)/(m-mm)
                        if lo+1e-9<x<up-1e-9:pts.append(x)
            pts=sorted(set(pts))
            for x,y in zip(pts,pts[1:]):
                m,b,ids=max(lines,key=lambda row:row[0]*(x+y)/2+row[1]);new=dict(lo=x,hi=y,slope=m,intercept=b,interval=ids)
                if pieces and (pieces[-1]['slope'],pieces[-1]['intercept'])==(m,b):pieces[-1]['hi']=y
                else:pieces.append(new)
        if len(pieces)>=2:
            for piece in pieces:
                x=(piece['lo']+piece['hi'])/2;E=interval_frontier(g,l,x,6,1,1)['E'];lp,_=solve_lp(g,l,x,6)
                assert abs(E-(piece['slope']*x+piece['intercept']))<1e-7 and abs(E-lp)<1e-7
            return dict(q=q.tolist(),demand=l.tolist(),charging_floor=floor,pieces=pieces)
    raise RuntimeError('No diagnostic breakpoint case found')

def run():
    matrix,memory=benchmark();p,w,l,v=generate(SEEDS[0]);g=300*(p+w)
    extras={}
    extras['noncyclic']=modified_lp(g,l,200,160,cyclic=False)
    extras['empty_initial']=modified_lp(g,l,200,160,initial=0)
    extras['self_discharge']=[dict(retention=r,**modified_lp(g,l,200,160,retention=r)) for r in [1.,.9999,.999]]
    extras['power_stress']=[dict(Pc=pc,Pd=pd,actual=interval_frontier(g,l,pc,pd),omit_charge=interval_frontier(g,l,1e6,pd),omit_discharge=interval_frontier(g,l,pc,1e6)) for pc,pd in [(40,160),(200,100),(75,160)]]
    # Fixed-E feasibility classifications on both sides of annual boundaries.
    classifications=[]
    for seed in SEEDS[:5]:
        pp,ww,ll,_=generate(seed);gg=300*(pp+ww);E=interval_frontier(gg,ll,200,160)['E']
        for factor in [.99,1.,1.01]:
            val,res=solve_lp(gg,ll,200,160,ec=.95,ed=.95,Efixed=E*factor)
            assert bool(np.isfinite(val))==(factor>=1)
            classifications.append(dict(seed=seed,factor=factor,status=res.status))
    # Joint convexity checks at random feasible design pairs, using zero-shortfall
    # result and endpoint feasibility (not an assumption of globally finite E).
    rng=np.random.default_rng(991);convex=[]
    for _ in range(100):
        alpha1,alpha2=rng.uniform(0,1,2);pc1,pc2=rng.uniform(80,300,2);lam=rng.uniform(0,1)
        E1=interval_frontier(600*(alpha1*p+(1-alpha1)*w),l,pc1,160)['E'];E2=interval_frontier(600*(alpha2*p+(1-alpha2)*w),l,pc2,160)['E']
        if not(np.isfinite(E1) and np.isfinite(E2)):continue
        am=lam*alpha1+(1-lam)*alpha2;pm=lam*pc1+(1-lam)*pc2;Em=interval_frontier(600*(am*p+(1-am)*w),l,pm,160)['E'];gap=Em-lam*E1-(1-lam)*E2
        assert gap<1e-6
        convex.append(dict(gap=float(gap)))
    (OUT/'extended_checks.json').write_text(json.dumps(dict(matrix=matrix,memory=memory,extensions=extras,classifications=classifications,convexity=convex,breakpoints=exact_small_breakpoints()),indent=2))
    print('EXTENDED CHECKS SAVED',flush=True)
if __name__=='__main__':run()
