"""Generate only absent pass5 energies with existing transport/BCA cache."""
from pathlib import Path
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from g4_beam.transport import transport as T
from g4_beam.merit.api import ENERGIES
HERE=Path(__file__).resolve().parent
OUT=HERE/'profiles.json'


def one(job):
    mat,e=job; start=time.monotonic()
    m=T.Model(dict(material='Pd' if mat.startswith('Pd') else 'Ti',loaded=mat in ('PdD','TiD'),E_cm=e))
    # Keep the original run history and its12h/2h boundary; no energy/time alias.
    tr=m.c['t_run'];te=np.unique(np.r_[600.,1800.,3600.,np.linspace(0,tr,13)[1:]])
    m.run_te(te[te<=tr]);layers,nads,state=m.profile()
    row=dict(case=dict(condition=mat,E_cm_keV=e,beam_time_s=float(m.t_beam[-1]),T_K=float(m.T)),
             layers=np.asarray(layers).tolist(),nads_m2=float(nads),
             state={k:np.asarray(v).tolist() for k,v in state.items()},
             source='pass5 generated missing energy; retained pass4 surface inputs',wall_s=time.monotonic()-start)
    return row


def main():
    start=time.monotonic()
    supplied=json.loads((HERE.parent/'forward/inputs/profiles_baseline.json').read_text())
    records=supplied['records']
    if OUT.exists():records=json.loads(OUT.read_text())['records']
    found={(r['case']['condition'],r['case']['E_cm_keV']) for r in records}
    jobs=[(m,float(e)) for m in ('Pd','PdD','Ti','TiD') for e in ENERGIES if (m,float(e)) not in found]
    if '--pilot' in sys.argv:jobs=[j for j in jobs if j==('PdD',1.)]
    with ProcessPoolExecutor(max_workers=1 if '--pilot' in sys.argv else 4) as pool:
        for row in pool.map(one,jobs):
            records.append(row)
            OUT.write_text(json.dumps(dict(units='SI; energies CMkeV; states bound to case run end',records=records),separators=(',',':')))
            print(row['case'],'wall_s',row['wall_s'],flush=True)
    print('generated',len(jobs),'wall_s',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
