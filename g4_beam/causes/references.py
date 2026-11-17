"""Explicit diagnostic reference scenarios, without changing canonical merit.

Both sides use the SAME refitted surface parameters. The zero-sheet projection
removes only an actually exposed Pd adsorbate, independently in each parallel
column. It does not alter bulk, film D, stopping, area, beam or material identity.
An oxide's separately assumed outer adsorbate is not relabelled as a Pd sheet.
This projection is NOT the separately authorized joined transport-boundary run.
"""
from functools import lru_cache
from pathlib import Path
import json
from g4_beam.merit.api import Target,reference,evaluate_fit

HERE=Path(__file__).resolve().parent
ZERO_SHEET_SURFACE='refitted_zero_sheet_surface'


def at(profiles,e):
    return profiles(float(e)) if callable(profiles) else profiles[float(e)]


def exposed_pd(profile):
    layers,_,_=profile
    return bool(len(layers) and float(layers[0][4])==46.)


def without_exposed_pd_sheet(target):
    def project(profiles):
        def profile(e):
            layers,nads,state=at(profiles,e)
            remove=exposed_pd((layers,nads,state))
            return layers,0. if remove else nads,dict(state,reference_scenario=ZERO_SHEET_SURFACE,
                    exposed_Pd_sheet_removed=remove,bulk_and_film_D_unchanged=True)
        return profile
    columns=tuple((weight,project(profiles)) for weight,profiles in target.columns)
    return Target(target.name,project(target.profiles),target.material,target.stopping_fn,columns,
                  dict(target.measurement),target.Ue_shift_keV)


@lru_cache(None)
def zero_sheet_fit():
    rows=json.loads((HERE/'density1d_fits.json').read_text())
    return next(row['fit'] for row in rows if row['cause']=='surface' and row['kind']=='sheet0')


@lru_cache(None)
def zero_sheet_reference():
    return without_exposed_pd_sheet(reference('PdD'))


def zero_sheet_yields(target,energies):
    fit=zero_sheet_fit()
    return evaluate_fit(without_exposed_pd_sheet(target),energies,'surface',fit)['physical_yield_per_C']


def zero_sheet_merits(target,energies):
    fit=zero_sheet_fit()
    denominator=evaluate_fit(zero_sheet_reference(),energies,'surface',fit)['physical_yield_per_C']
    return zero_sheet_yields(target,energies)/denominator


def pd_sheet_identity(target,e):
    columns=target.columns or ((1.,target.profiles),)
    exposed=retained=inventory=0.;outer=[]
    for weight,profiles in columns:
        layers,nads,state=at(profiles,e)
        is_pd=exposed_pd((layers,nads,state))
        if is_pd:
            exposed+=weight;inventory+=weight*nads
            if nads>0:retained+=weight
        if weight>0:outer.append(dict(area_fraction=weight,outer_Z=float(layers[0][4]),
                                     exposed_Pd=is_pd,adsorbed_D_m2=float(nads)))
    return dict(beam_side_exposed_Pd_area_fraction=float(exposed),
                retains_beam_side_Pd_sheet_area_fraction=float(retained),
                beam_side_Pd_sheet_D_m2=float(inventory),
                outer_column_identity=json.dumps(outer,separators=(',',':')))
