"""Frozen PdD local well. Cartesian Å -> eV and eV/Å, no neighbour term.

The radial Hermite interpolants anchor both E and -F. Old quartic is used only
below .03 Å, retaining the measured near K and a regular polynomial origin.
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'inputs'
OLD = json.loads((INPUT / 'old-readout.json').read_text())
WIDE = json.loads((INPUT / 'bulk-wells.json').read_text())
K = OLD['coupling']['on_site_k2_eV_A2']
A = OLD['directions']['100']['k4_eV_A4']
Q111 = OLD['directions']['111']['k4_eV_A4']
B = (3*Q111-A)/2
RAYS = ('100', '110', '111')
VECTORS = np.array([[1.,0,0], [1.,1,0], [1.,1,1]])
VECTORS /= np.linalg.norm(VECTORS,axis=1)[:,None]
K4 = np.array([A,(A+3*Q111)/4,Q111])
VARIANTS = ('old','quadratic','quartic')
R0, RMAX = .03, .48


def anchors():
    result = []
    for j,axis in enumerate(RAYS):
        knots = []
        for r in (.12,.24,.48):
            rows = [v for v in WIDE['rows'] if v['material']=='PdD' and
                    v['direction']==axis and abs(abs(v['displacement_A'])-r)<1e-8]
            native = [v for v in rows if v['sample_kind']=='native_scf']
            if native:
                E = np.mean([v['delta_energy_free_eV'] for v in native])
                slope = np.mean([-np.sign(v['displacement_A'])*v['native_force_parallel_eV_A'] for v in native])
                source = [v['case'] for v in native]
            else:
                old = next(v for v in OLD['directions'][axis]['points'] if v['displacement_A']==r)
                E, slope = old['delta_energy_free_eV'], -old['force_parallel_eV_A']
                source = [old['point_id']]
            knots.append(dict(r_A=r,energy_eV=float(E),slope_eV_A=float(slope),native_cases=source,
                              independent_native_count=len(native) if native else 1))
        result.append(knots)
    return result


ANCHORS = anchors()
SPLINES = []
for j,points in enumerate(ANCHORS):
    r = np.array([R0]+[p['r_A'] for p in points])
    e = [.5*K*R0**2+K4[j]*R0**4/4]+[p['energy_eV'] for p in points]
    s = [K*R0+K4[j]*R0**3]+[p['slope_eV_A'] for p in points]
    SPLINES.append(CubicHermiteSpline(r,e,s,extrapolate=False))


def radial(r,variant):
    r = np.asarray(r,float)
    E = np.empty(r.shape+(3,)); S = np.empty_like(E)
    for j,spline in enumerate(SPLINES):
        if variant=='old':
            E[...,j]=.5*K*r*r+K4[j]*r**4/4
            S[...,j]=K*r+K4[j]*r**3
            continue
        if variant not in ('quadratic','quartic'):
            raise ValueError(f'Unknown variant {variant!r}')
        rc=np.clip(r,R0,RMAX)
        e,s=spline(rc),spline(rc,1)
        small=r<R0
        e=np.where(small,.5*K*r*r+K4[j]*r**4/4,e)
        s=np.where(small,K*r+K4[j]*r**3,s)
        t=np.maximum(r-RMAX,0.)
        er,sr=ANCHORS[j][-1]['energy_eV'],ANCHORS[j][-1]['slope_eV_A']
        if variant=='quadratic':
            ee=er+sr*t+.5*K*t*t; ss=sr+K*t
        else:
            # Old quartic curvature and all higher Taylor terms, shifted to
            # match the new outer value and slope; no old linear offset.
            c=K+3*K4[j]*RMAX**2
            ee=er+sr*t+.5*c*t*t+K4[j]*RMAX*t**3+K4[j]*t**4/4
            ss=sr+c*t+3*K4[j]*RMAX*t*t+K4[j]*t**3
        E[...,j]=np.where(r>RMAX,ee,e)
        S[...,j]=np.where(r>RMAX,ss,s)
    return E,S


def energy_force(x,variant):
    """x (...,3) in Å; explicit variant; E (...) eV, F (...,3) eV/Å."""
    x=np.asarray(x,float)
    r2=np.sum(x*x,axis=-1);r=np.sqrt(r2)
    safe=np.where(r2>0,r2,1.)
    z=x*x/safe[...,None]
    q=z[...,0]*z[...,1]+z[...,1]*z[...,2]+z[...,2]*z[...,0]
    p=np.prod(z,axis=-1)
    w=np.stack([1-4*q+9*p,4*q-36*p,27*p],axis=-1)
    E,S=radial(r,variant)
    energy=np.sum(w*E,axis=-1)
    dq_dz=np.sum(z,axis=-1)[...,None]-z
    dp_dz=np.stack([z[...,1]*z[...,2],z[...,0]*z[...,2],z[...,0]*z[...,1]],axis=-1)
    dq=2*x/safe[...,None]*(dq_dz-2*q[...,None])
    dp=2*x/safe[...,None]*(dp_dz-3*p[...,None])
    dEdq=-4*E[...,0]+4*E[...,1]
    dEdp=9*E[...,0]-36*E[...,1]+27*E[...,2]
    gradient=x/np.where(r>0,r,1.)[...,None]*np.sum(w*S,axis=-1)[...,None]
    gradient+=dEdq[...,None]*dq+dEdp[...,None]*dp
    # Use the exact polynomial near the origin, avoiding cancellation of
    # formally absent angular terms and making smoothness manifest.
    near=(r<=R0) | (variant=='old')
    pe=.5*K*r2+A*np.sum(x**4,axis=-1)/4+B/2*(x[...,0]**2*x[...,1]**2+x[...,1]**2*x[...,2]**2+x[...,2]**2*x[...,0]**2)
    pg=K*x+A*x**3+B*x*(r2[...,None]-x*x)
    return np.where(near,pe,energy),-np.where(near[...,None],pg,gradient)
