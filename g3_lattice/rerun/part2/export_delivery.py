"""Preserve every scalar/configuration; omit only unplotted long curves.
Original full outputs stay on the producing host. The two figure-input NPZs
retain every plotted numeric curve in float64, and display-only slices have
explicitly recorded float32/complex64 export errors. No dynamics is run.
"""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def compact(value,path,omitted):
    if isinstance(value,dict):
        return {key:compact(item,path+'/'+key,omitted) for key,item in value.items()}
    if isinstance(value,list):
        curve=any(part in path for part in ('/spectrum/','/half_pump_correlations/',
            '/matched_legacy_filter/','/joint_harmonic_filter/','/equal_time_16cycle_envelopes/'))
        if curve and len(value)>64 and all(isinstance(v,(int,float)) for v in value):
            digest=hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()
            omitted.append({'path':path,'points':len(value),'values_SHA256':digest})
            return {'unplotted_long_curve_omitted':True,'points':len(value),'values_SHA256':digest}
        return [compact(item,path+'/'+str(i),omitted) for i,item in enumerate(value)]
    return value


def export(source,destination):
    original=(ROOT/source).read_bytes();omitted=[]
    data=compact(json.loads(original),'',omitted)
    result={'full_original_file':source,'full_original_SHA256':hashlib.sha256(original).hexdigest(),
        'storage_scope':'All scalar/short summaries and configurations preserved verbatim. Only long numeric correlation/spectral curves omitted; originals retained on producing host. Every plotted curve is included in figure_input.npz or ni_figure_input.npz.',
        'omitted_unplotted_curves':omitted,'result':data}
    (ROOT/destination).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return len(omitted)


def main():
    counts={
        'Pd':export('updated/pd_results.json','updated/pd_summary.json'),
        'Ni':export('updated/ni_epoch/phase.json','updated/ni_epoch/phase_summary.json')}
    pd=json.loads((ROOT/'updated/pd_results.json').read_text())
    rows=pd['results'] if isinstance(pd,dict) else pd
    ni=json.loads((ROOT/'updated/ni_epoch/phase.json').read_text())
    temperature={'definition':'Configured Langevin bath temperature, not an inferred driven-state/quantum temperature',
        'PdD':[{'label':r['label'],'temperature_K':r['configuration']['temperature_K']} for r in rows],
        'NiH':{'label':ni['summary']['label'],'temperature_K':ni['summary']['configuration']['temperature_K']}}
    (ROOT/'updated/temperatures.json').write_text(json.dumps(temperature,indent=2)+'\n')
    print(json.dumps({'curve_exports':counts,'Pd_cases':len(rows),'Ni_bath_K':temperature['NiH']['temperature_K']}))


if __name__=='__main__':main()
