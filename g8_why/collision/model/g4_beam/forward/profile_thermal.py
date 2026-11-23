"""Thermal-track density moments of the supplied frozen transport profiles.

Track assumption: volume-averaged N_D^2 along straight-ahead CSDA stopping
path, excluding the separately represented surface adsorbate. The unresolved
track radius, lifetime and branching belong to fitted COMMON amplitude A.
At unsupplied energies, log-linear interpolation of the positive density moment
is an explicit steady-profile interpolation, not target time evolution.
"""
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from g4_beam.transport import bca
from . import model as M
from . import thermal

HOT_KEV={'Pd':.400,'PdD':.400,'Ti':.370,'TiD':.370,'Zr':.340}


def transport_stopping(E_lab_keV,layer):
    """Supplier BCA/LS/ZBL cross sections eV nm² -> eV/m using SI density."""
    _,_,nh,nd,z,m,_=layer
    e=E_lab_keV*1000
    return float((nh*(bca.se_ls(e,z)+bca.sn_zbl(e,z,m))
                  +nd*(bca.se_ls(e,1)+bca.sn_zbl(e,1,M.M_D)))*1e-18)


def density_squared_average(layers,depth_m):
    layers=np.asarray(layers,float)
    if depth_m<=0:raise ValueError('Positive track depth required')
    weight=np.maximum(0,np.minimum(layers[:,1],depth_m)-layers[:,0])
    return float(np.sum(layers[:,3]**2*weight)/depth_m)


def track_density_moment(layers,E_lab_keV,E_cut=.1):
    energy=float(E_lab_keV);depth=0.
    if energy<=E_cut:return 0.
    for layer in layers:
        zt,zb=layer[:2]
        def rhs(z,y):return [-transport_stopping(max(y[0],E_cut),layer)/1000]
        def event(z,y):return y[0]-E_cut
        event.terminal=True;event.direction=-1
        sol=solve_ivp(rhs,(zt,zb),[energy],rtol=2e-6,atol=1e-7,events=event,max_step=10e-9)
        if not sol.success:raise RuntimeError(sol.message)
        energy=float(sol.y[0,-1]);depth=float(sol.t[-1])
        if sol.t_events[0].size:break
    return density_squared_average(layers,depth)


@lru_cache(maxsize=1)
def density_table():
    data=json.load(open(Path(__file__).resolve().parent/'inputs/profiles_baseline.json'))
    gx,gw=np.polynomial.hermite_e.hermegauss(5);gw/=gw.sum()
    result={m:[] for m in M.MATS}
    for r in data['records']:
        mat=r['case']['condition'];ec=r['case']['E_cm_keV'];layers=r['layers']
        # Beam-current fractions, not atom fractions; density moment is per e.
        moment=sum(w*(.85*track_density_moment(layers,2*ec*(1+.05/2.3548*g))
                      +.30*track_density_moment(layers,ec*(1+.05/2.3548*g))) for g,w in zip(gx,gw))
        result[mat].append((ec,moment))
    return {m:np.array(sorted(a)) for m,a in result.items()}


def charge_probability(material,E_cm_keV,A,T_track=1080.,E_R_eV=1.):
    """Per true incident charge e, not per detected count or cup charge."""
    table=density_table()[material];e=np.asarray(E_cm_keV,float)
    if np.any(e<table[0,0]-1e-12) or np.any(e>table[-1,0]+1e-12):
        raise ValueError('Thermal density requires supplied profile energy range .25-6.5keV CM')
    moment=np.exp(np.interp(np.log(e),np.log(table[:,0]),np.log(table[:,1])))
    host='Pd' if material.startswith('Pd') else 'Ti'
    return thermal.probability(np.sqrt(moment),Ue_hot_keV=HOT_KEV[material],E_A_eV=thermal.BARRIERS[host],exposure_m3_s=A,T_K=T_track,E_R_eV=E_R_eV)
