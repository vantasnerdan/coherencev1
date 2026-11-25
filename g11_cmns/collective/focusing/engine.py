"""Frozen PdD focusing Hamiltonian; eV, Å, fs. No live-source edits."""
from pathlib import Path
import csv
import importlib.util
import itertools
import json
import numpy as np
from numba import njit
from scipy.constants import atomic_mass, electron_volt
from scipy.interpolate import CubicHermiteSpline

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

M = module('focusing_original_model', ROOT/'g10_quantum_site/model.py')
WELL = module('focusing_updated_well', ROOT/'g10_quantum_site/bulk-refit-080/well.py')
BASE = M.full()
MASS = M.MASS['D']*atomic_mass*1e10/electron_volt
K, A, B = WELL.K, WELL.A, WELL.B
R0, RMAX = WELL.R0, .8
CELL = 8.168
HALF = CELL/4
AXIS = np.array([1., 1., 0.])/np.sqrt(2)
TI, TJ = [int(np.flatnonzero(BASE['original_indices'] == i)[0]) for i in (25, 63)]
JN = float(BASE['J'][TI, TJ])
PAIR_PATH = HERE.parent/'.pair-input/g10_quantum_site/pdd-nearest-pair/pair-cut.csv'
with PAIR_PATH.open() as stream:
    PAIR = list(csv.DictReader(stream))
QR = np.array([float(p['displacement_per_D_A']) for p in PAIR])
QE = np.array([float(p['connected_free_energy_eV']) for p in PAIR])
QS = np.array([float(p['connected_dE_du_eV_A']) for p in PAIR])
PAIR_SPLINE = CubicHermiteSpline(QR, QE, QS, extrapolate=False)
RESID = QE[0]-JN*QR[0]**2
SRESID = QR[0]*QS[0]-2*JN*QR[0]**2
C6 = (SRESID-4*RESID)/2
C4 = RESID-C6
RX = WELL.SPLINES[0].x.copy()
RC = np.stack([s.c for s in WELL.SPLINES])
PC = PAIR_SPLINE.c.copy()

@njit(cache=True)
def pair_correction(s):
    """Return ΔW(s) and its derivative; zero quadratic Hessian at origin."""
    if s <= 0:
        return 0., 0.
    if s < QR[0]**2:
        z = s/QR[0]**2
        return C4*z*z+C6*z*z*z, (2*C4*z+3*C6*z*z)/QR[0]**2
    q = np.sqrt(s)
    j = min(np.searchsorted(QR, q)-1, len(QR)-2)
    j = max(j, 0)
    t = q-QR[j]
    e = ((PC[0,j]*t+PC[1,j])*t+PC[2,j])*t+PC[3,j]
    slope = (3*PC[0,j]*t+2*PC[1,j])*t+PC[2,j]
    return e-JN*s, slope/(2*q)-JN

@njit(cache=True)
def local(x, force, harmonic):
    energy = 0.
    ee, ss = np.empty(3), np.empty(3)
    for i in range(len(x)):
        r2 = x[i,0]**2+x[i,1]**2+x[i,2]**2
        r = np.sqrt(r2)
        if harmonic or r <= R0:
            energy += .5*K*r2
            if not harmonic:
                energy += A*(x[i,0]**4+x[i,1]**4+x[i,2]**4)/4
                energy += B*(x[i,0]**2*x[i,1]**2+x[i,1]**2*x[i,2]**2+x[i,2]**2*x[i,0]**2)/2
            for a in range(3):
                force[i,a] = -K*x[i,a]
                if not harmonic:
                    force[i,a] -= A*x[i,a]**3+B*x[i,a]*(r2-x[i,a]**2)
            continue
        z0, z1, z2 = x[i,0]**2/r2, x[i,1]**2/r2, x[i,2]**2/r2
        aq = z0*z1+z1*z2+z2*z0
        ap = z0*z1*z2
        weights = (1-4*aq+9*ap, 4*aq-36*ap, 27*ap)
        j = max(0, min(np.searchsorted(RX, min(r,RMAX))-1, len(RX)-2))
        t = min(r,RMAX)-RX[j]
        for a in range(3):
            ee[a] = ((RC[a,0,j]*t+RC[a,1,j])*t+RC[a,2,j])*t+RC[a,3,j]
            ss[a] = (3*RC[a,0,j]*t+2*RC[a,1,j])*t+RC[a,2,j]
            # A continuation only makes the first crossing step numerically defined.
            # No post-exit trajectory is scored or interpreted.
            if r > RMAX:
                dr = r-RMAX
                ee[a] += ss[a]*dr+.5*K*dr*dr
                ss[a] += K*dr
        radial = 0.
        for a in range(3):
            energy += weights[a]*ee[a]
            radial += weights[a]*ss[a]
        deq = 4*(ee[1]-ee[0])
        dep = 9*ee[0]-36*ee[1]+27*ee[2]
        zs = (z0,z1,z2)
        for a in range(3):
            dq = 2*x[i,a]/r2*(1-zs[a]-2*aq)
            dp = 2*x[i,a]/r2*(zs[(a+1)%3]*zs[(a+2)%3]-3*ap)
            force[i,a] = -x[i,a]/r*radial-deq*dq-dep*dp
    return energy

@njit(cache=True)
def energy_force(x, force, edges, weights, near, directions, harmonic=False):
    energy = local(x, force, harmonic)
    for b in range(len(edges)):
        i,j = edges[b,0], edges[b,1]
        w = weights[b]
        for a in range(3):
            energy -= w*x[i,a]*x[j,a]
            force[i,a] += w*x[j,a]
            force[j,a] += w*x[i,a]
    if not harmonic:
        for b in range(len(near)):
            i,j = near[b,0],near[b,1]
            u,v = 0.,0.
            for a in range(3):
                u += directions[b,a]*x[i,a]
                v -= directions[b,a]*x[j,a]
            if u > 0 and v > 0:
                e,d = pair_correction(u*v)
                energy += e
                for a in range(3):
                    force[i,a] -= d*v*directions[b,a]
                    force[j,a] += d*u*directions[b,a]
    return energy


def stencil():
    """Disaggregate each original folded bond among equal shortest images."""
    p = BASE['positions']
    out = []
    for j in range(1,len(p)):
        delta = p[j]-p[0]
        values = []
        for a in range(3):
            d = delta[a]-np.rint(delta[a]/CELL)*CELL
            values.append((-CELL/2,CELL/2) if np.isclose(abs(d),CELL/2) else (d,))
        images = list(itertools.product(*values))
        for d in images:
            out.append((np.rint(np.array(d)/HALF).astype(int), float(BASE['J'][0,j])/len(images)))
    return out


class Lattice:
    def __init__(self, m):
        if m not in (2,4,8):
            raise ValueError('conventional cell count must be 2,4,8')
        self.m, self.box = m, m*CELL/2
        tiles = list(itertools.product(range(m//2), repeat=3))
        self.positions = np.concatenate([BASE['positions']+CELL*np.array(t) for t in tiles])
        self.grid = np.rint(self.positions/HALF).astype(int)
        self.n = len(self.positions)
        side = 2*m
        lookup = {tuple(p):i for i,p in enumerate(self.grid)}
        edge_map, near_map = {}, {}
        self.stencil = stencil()
        for i,g in enumerate(self.grid):
            for d,w in self.stencil:
                j = lookup[tuple((g+d)%side)]
                if i < j:
                    edge_map[i,j] = edge_map.get((i,j),0.)+w
                    if np.dot(d,d) == 2:
                        near_map[i,j] = d*HALF
        self.edges = np.array(list(edge_map),dtype=np.int64)
        self.weights = np.array(list(edge_map.values()))
        self.near = np.array(list(near_map),dtype=np.int64)
        self.bonds = np.array(list(near_map.values()))
        self.directions = self.bonds/np.linalg.norm(self.bonds,axis=1)[:,None]
        self.target = (TI,TJ)
        self.target_bond = BASE['positions'][TJ]-BASE['positions'][TI]
        self.target_axis = self.target_bond/np.linalg.norm(self.target_bond)
        kk = np.meshgrid(*([np.fft.fftfreq(side)*2*np.pi]*3),indexing='ij')
        symbol = np.full((side,side,side),K)
        for d,w in self.stencil:
            symbol -= w*np.cos(sum(kk[a]*d[a] for a in range(3)))
        self.symbol = symbol
        self.omega = np.sqrt(symbol/MASS)
        seed = np.zeros((side,side,side,3))
        seed[tuple(self.grid[TI])] = self.target_axis
        seed[tuple(self.grid[TJ])] = -self.target_axis
        self.seed_fft = np.fft.fftn(seed, axes=(0,1,2))

    def predicted_kernel(self, t):
        field = np.fft.ifftn(self.seed_fft*np.cos(self.omega*t)[...,None],axes=(0,1,2)).real
        return field[tuple(self.grid.T)]

    def energy(self,x,v=None,harmonic=False):
        f = np.empty_like(x)
        e = energy_force(x,f,self.edges,self.weights,self.near,self.directions,harmonic)
        return float(e+(0 if v is None else .5*MASS*np.sum(v*v)))

@njit(cache=True)
def evolve(x0,v0,edges,weights,near,directions,bonds,target_bond,ti,tj,dt,steps,damping,harmonic):
    """All-step support checks; stop before unsupported points enter scoring."""
    x,v = x0.copy(),v0.copy()
    f, fnew = np.empty_like(x),np.empty_like(x)
    potential = energy_force(x,f,edges,weights,near,directions,harmonic)
    # Plus actual local V and old harmonic V at the most-displaced site.
    trace = np.empty((steps+1,12))
    wall_force = np.empty((1,3))
    exit_time, exit_kind = -1.,0
    count = 0
    for step in range(steps+1):
        maxsite,minbond,whole = 0.,1e30,0.
        imax = 0
        for i in range(len(x)):
            radius_i = np.sqrt(np.dot(x[i],x[i]))
            if radius_i > maxsite:
                maxsite,imax = radius_i,i
        for b in range(len(near)):
            i,j = near[b,0],near[b,1]
            rr,vv = 0.,0.
            for a in range(3):
                dr = bonds[b,a]+x[j,a]-x[i,a]
                rr += dr*dr
                vv += (v[i,a]-v[j,a])**2
            minbond = min(minbond,np.sqrt(rr))
            whole = max(whole,.25*MASS*vv)
        dv2,rr,rad = 0.,0.,0.
        for a in range(3):
            dr = target_bond[a]+x[tj,a]-x[ti,a]
            rr += dr*dr
            dv2 += (v[ti,a]-v[tj,a])**2
            rad += (v[ti,a]-v[tj,a])*dr
        radius = np.sqrt(rr)
        closing = rad/radius
        kinetic = .5*MASS*np.sum(v*v)
        local_wall = local(x[imax:imax+1],wall_force,harmonic)
        trace[step] = (step*dt,.25*MASS*dv2,.25*MASS*closing**2,radius,maxsite,minbond,kinetic+potential,potential,whole,closing,local_wall,.5*K*maxsite*maxsite)
        count = step+1
        if maxsite > RMAX or minbond < 1.28782409436586:
            exit_time,exit_kind = step*dt,1 if maxsite>RMAX else 2
            break
        if step == steps:
            break
        decay = np.exp(-dt/damping) if damping > 0 else 1.
        for i in range(len(x)):
            for a in range(3):
                v[i,a] += .5*dt/MASS*f[i,a]
                x[i,a] += .5*dt*v[i,a]
                v[i,a] *= decay
                x[i,a] += .5*dt*v[i,a]
        potential = energy_force(x,fnew,edges,weights,near,directions,harmonic)
        for i in range(len(x)):
            for a in range(3):
                v[i,a] += .5*dt/MASS*fnew[i,a]
        f,fnew = fnew,f
    return x,v,trace[:count],exit_time,exit_kind


def trajectory(lattice,q,v,dt=.5,horizon=2400.,damping=0.,harmonic=False):
    return evolve(q,v,lattice.edges,lattice.weights,lattice.near,lattice.directions,lattice.bonds,
                  lattice.target_bond,*lattice.target,dt,int(round(horizon/dt)),damping,harmonic)


def launch_energy_diagnostics(lattice,q,v):
    kinetic = .5*MASS*np.sum(v*v,axis=1)
    local_energy = WELL.energy_force(q,'wide_quadratic')[0]
    positive = kinetic+local_energy
    total = float(positive.sum())
    return dict(
        participation_positive_site_energy=float(total*total/np.sum(positive*positive)) if total else 0.,
        target_positive_site_energy_fraction=float(positive[list(lattice.target)].sum()/total) if total else 0.,
        positive_site_energy_definition='K_i+updated local V_i; signed bond and connected energies excluded',
        launch_max_site_kinetic_eV=float(kinetic.max()),
        launch_max_site_local_potential_eV=float(local_energy.max()))


def summarize(lattice,q,v,result,energy,pattern,name,scope='target'):
    x,vf,trace,exit_time,exit_kind = result
    supported = trace if exit_time < 0 else trace[:-1]
    if not len(supported):
        supported = trace[:1]
    column = 1 if scope=='target' else 8
    peak = int(np.argmax(supported[:,column]))
    value = supported[peak,column]
    lo,hi = peak,peak
    while lo>0 and supported[lo-1,column]>=value/2: lo-=1
    while hi+1<len(supported) and supported[hi+1,column]>=value/2: hi+=1
    site_k = .5*MASS*np.sum(v*v,axis=1)
    participation = site_k.sum()**2/np.sum(site_k**2) if np.any(site_k) else 0.
    total0 = trace[0,6]
    preload = .5*MASS*np.sum(v[list(lattice.target)]**2)/total0 if total0>0 else 0.
    return dict(n_D=lattice.n,energy_per_D_eV=energy,pattern=pattern,candidate=name,scope=scope,
        peak_relative_kinetic_eV=float(value),relative_radial_kinetic_eV_at_peak=float(supported[peak,2]) if scope=='target' else None,
        first_site_or_bond_exit_fs=None if exit_time<0 else float(exit_time),exit_kind={0:None,1:'site_radius',2:'bond_distance'}[exit_kind],
        candidate_supported=bool(exit_time<0 or supported[peak,0]<exit_time),trajectory_supported=bool(exit_time<0),
        peak_before_support_exit=bool(exit_time<0 or supported[peak,0]<exit_time),
        peak_time_fs=float(supported[peak,0]),half_peak_hold_fs=float(supported[hi,0]-supported[lo,0]),
        hold_left_censored=bool(lo==0),hold_right_censored=bool(hi==len(supported)-1),
        participation_at_launch=float(participation),participation_definition='kinetic inverse participation',
        target_kinetic_preload_fraction=float(preload),potential_fraction_at_launch=float(trace[0,7]/total0),
        pair_gain=float(value/energy),fraction_total_energy_at_peak=float(value/(lattice.n*energy)),
        separation_A_at_peak=float(supported[peak,3]) if scope=='target' else None,
        signed_inward_speed_A_per_fs_at_peak=float(supported[peak,9]) if scope=='target' else None,
        max_site_radius_A_before_exit=float(supported[:,4].max()),min_bond_A_before_exit=float(supported[:,5].min()),
        horizon_fs=float(trace[-1,0]),relative_energy_drift=float(np.max(np.abs(supported[:,6]-total0))/total0),
        launch_total_energy_eV=float(total0),peak_index=peak)
