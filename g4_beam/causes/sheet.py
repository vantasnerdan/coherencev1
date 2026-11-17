"""Pd/PdD sheet occupancy profiles at fixed original m=1 bulk.
Absolute ML axes include0; all three causes use the shared joint likelihood.
"""
from pathlib import Path
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv
import json
import os
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.merit.api import Target,reference,bestfits,CAUSES,ENERGIES,FOUR,evaluate_fit,write_csv
from g4_beam.merit.responses import Response
from g4_beam.forward import model as M
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.causes.likelihood import joint_fit
HERE=Path(__file__).resolve().parent
THETA=np.linspace(0,1,6);NS=1.53e19
DATA=M.load_data(str(HERE.parent/'forward/errbars.csv'))
WORKERS=int(os.environ.get('G5_SHEET_WORKERS','4'))
if not 1<=WORKERS<=4:raise ValueError('G5_SHEET_WORKERS must be1–4')


@lru_cache(maxsize=16)
def occupied_target(material,theta):
    base=reference(material)
    if material not in ('Pd','PdD'):return base
    def profile(e):
        layers,nads,state=base.profiles(e)
        return layers,theta*NS,dict(state,sheet_occupancy_ML=theta,bulk_profile='original m1 boundary; diagnostic not self-consistent new transport')
    return Target(f'{material}_theta{theta}',profile,material,transport_stopping)


@lru_cache(maxsize=56)
def bulk_response(material,e):
    layers,_,state=reference(material).profiles(e)
    return Response(e,(layers,0.,state))


class SheetResponse(Response):
    def __init__(self,material,e,theta):
        super().__init__(e,occupied_target(material,theta).profiles(e))
        self.bulk=bulk_response(material,e)
    def build_neutral_table(self):
        if self.bulk.neutral_tables is None:self.bulk.build_neutral_table()
        self.neutral_tables={} # completed linear bulk+sheet composition, no duplicate stopping table
    def neutral(self,e,u,ch):
        if self.neutral_tables is None:self.build_neutral_table()
        return self.bulk.neutral(e,u,ch)+self.nads*float(M.reaction_sigma(max(e/2,.05),u,ch))*1e-31


@lru_cache(maxsize=16)
def occupied_responses(material,theta):
    if material in ('Pd','PdD'):return {float(e):SheetResponse(material,float(e),theta) for e in ENERGIES}
    return {float(e):Response(float(e),reference(material).profiles(float(e))) for e in ENERGIES}


@lru_cache(None)
def zero_sheet_fits():
    return {r['cause']:r['fit'] for r in json.loads((HERE/'density1d_fits.json').read_text()) if r['kind']=='sheet0'}


def point(job):
    cause,tp,td=job;start=time.monotonic();baseline=bestfits()[cause]
    targets={m:occupied_target(m,tp if m=='Pd' else td if m=='PdD' else 1.) for m in M.MATS}
    responses={m:occupied_responses(m,tp if m=='Pd' else td if m=='PdD' else 1.) for m in M.MATS}
    # Interpolate only optimizer seeds between completed full/zero-sheet fits;
    # this is not interpolation of a likelihood, observation or physical yield.
    full=baseline['params'];zero=zero_sheet_fits()[cause]['params'];mean=(tp+td)/2
    initial={n:mean*full[n]+(1-mean)*zero[n] for n in full}
    for material,theta in [('Pd',tp),('PdD',td)]:
        for prefix in ('Ue_','Ues_'):
            key=prefix+material
            if key in initial:initial[key]=theta*full[key]+(1-theta)*zero[key]
    result=joint_fit(cause,targets,responses,baseline,DATA,initial_params=initial)
    result.update(cause=cause,theta_Pd=float(tp),theta_PdD=float(td),Pd_sheet_D_m2=tp*NS,PdD_sheet_D_m2=td*NS,
                  bulk_multiplier_Pd=1.,bulk_multiplier_PdD=1.,wall_s=time.monotonic()-start,
                  Pd_sheet_fuller_than_PdD=bool(tp>td))
    return result


def flat(row):
    out={k:v for k,v in row.items() if k not in ('fit','model_ratio_by_energy','SI1_E_cm_keV')}
    out.update({f'ratio_E{e:g}':v for e,v in zip(row['SI1_E_cm_keV'],row['model_ratio_by_energy'])})
    return out


def summarize(rows):
    summaries=[]
    for cause in CAUSES:
        rr=[r for r in rows if r['cause']==cause];best=min(rr,key=lambda r:r['joint_chi2']);shape_min=min(r['shape_chi2'] for r in rr)
        for row in rr:
            row['shape_delta_chi2']=row['shape_chi2']-shape_min
            row['requested_region']=bool(row['shape_delta_chi2']<=4 and row['ratio_all_energies_in_measured_band'])
            row['geomean_region']=bool(row['shape_delta_chi2']<=4 and row['ratio_geomean_in_measured_band'])
        accepted=[r for r in rr if r['requested_region']]
        summaries.append(dict(cause=cause,joint_minimum=flat(best),shape_minimum=shape_min,region_count=len(accepted),
                         region_theta_pairs=[[r['theta_Pd'],r['theta_PdD']] for r in accepted],
                         region_Pd_fuller_pairs=[[r['theta_Pd'],r['theta_PdD']] for r in accepted if r['Pd_sheet_fuller_than_PdD']],
                         geomean_region_count=sum(r['geomean_region'] for r in rr),all_solver_success=all(r['success'] for r in rr)))
    write_csv(HERE/'sheet_grid.csv',[flat(r) for r in rows])
    with (HERE/'sheet_region.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(flat(rows[0])));writer.writeheader();writer.writerows(flat(r) for r in rows if r['requested_region'])
    (HERE/'sheet_summary.json').write_text(json.dumps(summaries,indent=2))
    (HERE/'sheet_fits.json').write_text(json.dumps(rows,indent=2))


def main():
    start=time.monotonic();path=HERE/'sheet_fits.json';rows=json.loads(path.read_text()) if path.exists() else []
    if '--pilot-sheet' in sys.argv:
        row=next((r for r in rows if (r['cause'],r['theta_Pd'],r['theta_PdD'])==('surface',0.,0.)),None)
        if row is None:
            row=point(('surface',0.,0.));rows.append(row);path.write_text(json.dumps(rows,indent=2))
        print(json.dumps(flat(row),indent=2),flush=True);return
    done={(r['cause'],r['theta_Pd'],r['theta_PdD']) for r in rows}
    jobs=[(c,float(a),float(b)) for c in CAUSES for a in THETA for b in THETA if (c,float(a),float(b)) not in done]
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(point,j):j for j in jobs}
        for future in as_completed(futures):
            row=future.result();rows.append(row);path.write_text(json.dumps(rows,indent=2))
            print('SHEET',row['cause'],row['theta_Pd'],row['theta_PdD'],'shape',row['shape_chi2'],'ratio',row['ratio_chi2'],'seconds',row['wall_s'],flush=True)
    summarize(rows)
    print('SHEET COMPLETE',len(rows),'wall_s',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
