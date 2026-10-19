"""Export original seeded inputs and verify selected archived results."""
from pathlib import Path
import sys, json, io, zipfile, gzip, hashlib, platform
import numpy as np
import scipy, networkx, matplotlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'reproducibility'))
from annual_study import generate, SEEDS, interval_frontier
from analysis import solve_lp

def csv_text(a,header):
    s=io.StringIO(); np.savetxt(s,a,delimiter=',',header=header,comments=''); return s.getvalue()

def main():
    out=ROOT/'datasets';out.mkdir(exist_ok=True)
    archived=json.loads((ROOT/'reproducibility/annual_results/annual_results.json').read_text())
    errors=[]
    def check(actual,expected):
        if np.isinf(expected): assert actual==expected; return
        error=abs(actual-expected);errors.append(error);assert error<1e-5,(actual,expected)
    header='hour,pv_pu,wind_pu,wind_speed_ms,load_kw,generation_kw'
    with zipfile.ZipFile(out/'annual_seed_profiles.zip','w',zipfile.ZIP_DEFLATED) as z:
        for seed in SEEDS:
            p,w,l,v=generate(seed);g=300*(p+w)
            data=np.c_[np.arange(len(l)),p,w,v,l,g]
            z.writestr(f'annual_seed_{seed}.csv',csv_text(data,header))
            for row in (r for r in archived['mix'] if r['seed']==seed):
                gg=600*(row['alpha']*p+(1-row['alpha'])*w)
                check(interval_frontier(gg,l,200,160)['E'],row['E'])
                check(interval_frontier(gg*(1.4*l.sum()/gg.sum()),l,200,160)['E'],row['normalized_E'])
    p,w,l,v=generate(SEEDS[0]);g=300*(p+w)
    original=np.loadtxt(ROOT/'reproducibility/annual_results/annual_base_profiles.csv',delimiter=',',skiprows=1)
    regenerated=np.c_[np.arange(len(l)),p,w,v,l,g]
    base_error=float(np.max(np.abs(original-regenerated)));assert base_error<1e-10
    scenarios=[('base',g,l)]
    for hours in [6,12,24,48,72]:
        gg=g.copy()
        for start in [1200,4300,7000]:gg[start:start+hours]*=.1
        scenarios.append((f'joint_low_{hours}h',gg,l))
    for mode in ['morning','flat']:
        pp,ww,ll,_=generate(SEEDS[0],mode=mode);scenarios.append((mode+' demand',300*(pp+ww),ll))
    for factor in [.7,1.3]:scenarios.append((f'wind_scale_{factor}',300*(p+np.clip(w*factor,0,1)),l))
    for rho in [.7,.99]:
        pp,ww,ll,_=generate(SEEDS[0],rho=rho);scenarios.append((f'persistence_{rho}',300*(pp+ww),ll))
    for cross in [-.8,.8]:
        pp,ww,ll,_=generate(SEEDS[0],cross=cross);scenarios.append((f'latent_cross_{cross}',300*(pp+ww),ll))
    with zipfile.ZipFile(out/'scenario_profiles.zip','w',zipfile.ZIP_DEFLATED) as z:
        for name,gg,ll in scenarios:
            z.writestr(name.replace(' ','_')+'.csv',csv_text(np.c_[np.arange(len(ll)),gg,ll],'hour,generation_kw,load_kw'))
            row=next(r for r in archived['scenarios'] if r['name']==name)
            check(interval_frontier(gg,ll,200,160)['E'],row['E'])
    pp,ww,ll,vv=generate(SEEDS[0],n=35040,dt=.25);gg=300*(pp+ww)
    with gzip.open(out/'quarter_hour_profiles.csv.gz','wt') as f:f.write(csv_text(np.c_[np.arange(len(ll))*.25,pp,ww,vv,ll,gg],header))
    check(interval_frontier(gg,ll,200,160,dt=.25)['E'],archived['quarter_hour']['E'])
    E,res=solve_lp(g,l,200,160,ec=.95,ed=.95)
    assert res.success
    check(E,archived['annual_checks'][0]['E'])
    report=dict(python=platform.python_version(),dependencies={m.__name__:m.__version__ for m in [np,scipy,networkx,matplotlib]},annual_seeds=len(SEEDS),mixture_frontiers_checked=660,scenarios_checked=len(scenarios),quarter_hour_steps=35040,base_profile_max_absolute_difference=base_error,max_finite_frontier_difference=max(errors),base_lp_energy_kwh=E,base_lp_success=bool(res.success),scope='Input regeneration and analytical frontier checks plus one base-case LP; full archived experiment suite and timings were not rerun.')
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    sums={str(f.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(ROOT.rglob('*')) if f.is_file() and '.git' not in f.parts and '__pycache__' not in f.parts and f.name!='SHA256SUMS.json'}
    (ROOT/'SHA256SUMS.json').write_text(json.dumps(sums,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
