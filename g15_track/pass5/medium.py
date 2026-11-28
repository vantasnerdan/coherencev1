"""Pass 5 heat and native optical landscape on the accepted finite foil.

Heat is a conditional Poisson mean, not a sampled microscopic temperature.
Carrier clocks use seconds since the electrical parent crossed the front face.
The deterministic single-parent electronic excitation is a separate observable.
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import csv
import hashlib
import importlib
import json
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal

BASELINE = 'f6010129f11ba8010815b78a67b3b6bed1db7194'
ENERGIES = (.25, .5, 1., 1.5, 2., 2.5)
BEAMS = ('even_steady', 'hot_steady', 'hot_bursts')
ALLOCATIONS = ('onsite', 'coupling')
MOUNTS = {'weak_edge': 1e3, 'nominal': 1e4, 'strong_edge': 1e5}
SPECIES = np.array([.85, .15])  # electrical parent fractions
MODULES = ('sources', 'fields', 'run', 'carriers', 'native', 'acoustic', 'launch', 'response', 'retarded')
MANIFESTS = {
    'manifest_core.json': '47f2c543f484419248afa5989d9cb4edbd4ccce9384c1c44d5f9b8e30a495883',
    'manifest_evidence.json': 'fad2f67b94f706dc14a77bfb0d36e561b556d3a8dd9a77c4e3c383d1c1fbcdf6',
}
SOURCE_FILES = tuple(f'source_E{ec:g}_nominal_retained_sp{sp}.npz' for ec in ENERGIES for sp in range(2))
NATIVE_FILES = tuple(f'native_E{ec:g}_nominal_retained_{a}_sp{sp}_m16.npz'
                     for ec in ENERGIES for a in ALLOCATIONS for sp in range(2))
PREDECESSOR_FILES = (*MANIFESTS, *(n+'.py' for n in MODULES), 'calendar.npz', 'native_summary.json',
                     'paths_nominal_retained_onsite.npz', 'frozen/pass1_inputs.json',
                     'frozen/tracks.json', 'frozen/g16_response.csv', 'frozen/guided_branches.csv',
                     *SOURCE_FILES, *NATIVE_FILES)
OMEGA_COX = 2*np.pi*1e5
DT = .001
LAG_STEP = 5e-5
T_AMBIENT = 300.
H_FACE = 5.
H_RAD = 4*.1*5.670374419e-8*T_AMBIENT**3
V0 = 17.029316176e-30
QE = 1.602176634e-19
KB = 8.617333262145e-5


def _sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _bind_inputs(path):
    entries = {}
    for name, sha in MANIFESTS.items():
        if _sha(path/name) != sha:
            raise ValueError(f'{name} does not match the accepted Pass 4 manifest')
        entries.update({r['path'].removeprefix('g15_track/pass4/'): r['sha256']
                        for r in json.loads((path/name).read_text())['files']})
    bindings = dict(MANIFESTS)
    for name in PREDECESSOR_FILES[len(MANIFESTS):]:
        expected = entries[name]
        if _sha(path/name) != expected:
            raise ValueError(f'{name} does not match the accepted Pass 4 bytes')
        bindings[name] = expected
    return bindings


@contextmanager
def _pass4(path):
    """Isolate the legacy bare imports and forbid writing pyc into inherited inputs."""
    aliases = {name: sys.modules.pop(name, None) for name in MODULES}
    bytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(path))
    try:
        yield {name: importlib.import_module(name) for name in MODULES}
    finally:
        sys.path.pop(0)
        sys.dont_write_bytecode = bytecode
        for name, old in aliases.items():
            sys.modules.pop(name, None)
            if old is not None:
                sys.modules[name] = old


class Clock:
    def __init__(self, data):
        self.edges = data['edges_s']
        self.state = data['burst_on']
        self.phase = float(data['phase'])
        self.norm = float(data['normalization'])
        self.seconds = float(self.edges[-1])
        self.burst = np.r_[0., np.cumsum(np.diff(self.edges)*self.state)]

    def index(self, t):
        return np.clip(np.searchsorted(self.edges, t, side='right')-1, 0, len(self.state)-1)

    def value(self, t):
        return np.where((t >= 0) & (t <= self.seconds),
                        (1+.7*np.cos(OMEGA_COX*t+self.phase)+9*self.state[self.index(t)])/self.norm, 0.)

    def primitive(self, t, beam):
        t = np.clip(t, 0., self.seconds)
        if beam != 'hot_bursts':
            return t
        i = self.index(t)
        return (t + .7/OMEGA_COX*(np.sin(OMEGA_COX*t+self.phase)-np.sin(self.phase))
                + 9*(self.burst[i]+(t-self.edges[i])*self.state[i]))/self.norm


def _fold(x, length):
    return length - abs(np.mod(x, 2*length)-length)


def _transport_xy(carrier, maps, delta):
    """Exact inherited finite-cell lateral translation; depth is independent."""
    n = maps.shape[-1]
    ix, iy = np.nonzero(maps[0])
    x, wx = carrier.axis_transfer(ix, delta[0]); y, wy = carrier.axis_transfer(iy, delta[1])
    weight = maps[:, ix, iy]
    result = np.zeros_like(maps)
    for a in (0, 1):
        for b in (0, 1):
            dest = x[a]*n+y[b]
            for i in range(len(maps)):
                result[i] += np.bincount(dest, weights=weight[i]*wx[a]*wy[b], minlength=n*n).reshape(n, n)
    return result


def _birth_rows(track, sp, origin, src):
    """Same depth GL2 as Pass 4, retaining birth depth/age correlation and both molecular starts."""
    z = np.r_[0., [r['depth_A']*1e-10 for r in track['timeline']]]
    t = np.r_[0., [r['time_fs']*1e-15 for r in track['timeline']]]
    if origin == 0:
        depth, age, w = z[-1:], t[-1:], np.ones(1)
    else:
        key = ('', 'host_recoil_eV', 'D_recoil_eV', 'electron_eV')[origin]
        energy = np.diff(np.r_[0., [r[key] for r in track['timeline']]])
        fraction = .5*(1+np.array([-1., 1.])/np.sqrt(3.))
        depth = (z[:-1, None]+np.diff(z)[:, None]*fraction).ravel()
        age = (t[:-1, None]+np.diff(t)[:, None]*fraction).ravel()
        w = np.repeat(energy/(2*energy.sum()), 2)
    if sp:
        delta = .74e-10/np.sqrt(2*track['initial_eV']*src.QE/src.MD)
        return np.tile(depth, 2), np.r_[age, age+delta], np.tile(w/2, 2)
    return depth, age, w


def _zbin(depth, zedges):
    return np.minimum(np.searchsorted(zedges[1:], _fold(depth, zedges[-1]), side='right'), len(zedges)-2)


def _lag_add(groups, age, distribution, power):
    """Two-bin temporal projection: preserve zeroth and first age moments."""
    integer, part = divmod(age/LAG_STEP, 1.)
    left = int(integer)
    for index, share in ((left, 1-part), (left+1, part)):
        if share == 0:
            continue
        if index not in groups:
            groups[index] = np.zeros_like(distribution)
        groups[index] += (share*power)*distribution


def _refined_acoustic(adapter, field, launch, zedges):
    """Apply the accepted loss characteristic and bath Green to the refined enthalpy grid.

    Inherit source, standing-guided/ballistic rays, loading and outgoing port;
    only the terminal deposition is changed from the 3-slab grid.
    """
    class Endpoints(adapter.Operator):
        def __init__(self, arm, ec, beam, tracks):
            self.groups = {}
            self.bath_e = np.zeros((field.N, field.N, len(zedges)-1))
            self.bath_e_age = 0.
            super().__init__(arm, ec, beam, tracks)

        def add_loss(self, age, points, amount, energies, native_fraction, port, pol, di):
            self.volume.append((age, field.deposit(points, amount), energies.copy(), native_fraction, port, pol, di))
            power = float(energies@SPECIES) * adapter.s.current(self.ec)/QE
            weights = np.repeat(amount/6, 6)
            directions = launch.DIRECTIONS
            for channel, part, tau, path in ((6, 1-.2*native_fraction, self.arm['pd'], self.arm['lambda_pd']),
                                              (8, .2*native_fraction, self.arm['te'], self.arm['lambda_e'])):
                if part == 0:
                    continue
                q = ((np.arange(6)+.5)/6)**(1/3)*np.pi/2
                speed = np.full(6, 1e6) if channel == 8 else np.array([2440.,3300.,4460.]*2)*np.cos(q)
                for duration in ((tau, tau/2) if channel == 8 else (tau,)):
                    flight = np.minimum(duration, path/speed)
                    distance = np.minimum(speed*duration, speed*flight+np.sqrt(2*speed*path*(duration-flight)))
                    loc = (points[:, None, :]+directions[None, :, :]*distance[None, :, None]).reshape(-1, 3)
                    xy = _fold(loc[:, :2]+adapter.s.L/2, adapter.s.L)
                    ij = np.minimum((xy/field.DX).astype(int), field.N-1)
                    k = _zbin(loc[:, 2], zedges)
                    shape = self.bath_e.shape
                    index = (ij[:,0]*field.N+ij[:,1])*shape[-1]+k
                    distribution = np.bincount(index, weights=weights, minlength=np.prod(shape)).reshape(shape)
                    if duration == tau:
                        _lag_add(self.groups, age+tau, distribution, power*part)
                        # Base class retains a record of the final release clock and port mass.
                        self.heat.append((age+tau, np.array([weights.sum()]), energies*part, native_fraction))
                    else:
                        self.bath_e += power*part*tau*distribution
                        self.bath_e_age += power*part*tau*distribution.sum()*(age+tau/2)
    return Endpoints


def _carrier_fields(mod, path, arm, tracks, xy, zedges):
    src, model, carrier, native = (mod[k] for k in ('sources', 'run', 'carriers', 'native'))
    with np.load(path/'paths_nominal_retained_onsite.npz') as data:
        ends, segments = data['heat_endpoints'], data['segments']
    residence = np.linalg.inv(-model.network(arm))
    transport = {}
    for o in range(4):
        for ch in (3, 8):
            records = segments[(segments[:,0] == o) & (segments[:,2] == ch)]
            if residence[ch, o] == 0 or len(records) == 0:
                continue
            duration = records[:,4]-records[:,3]
            weights = duration*records[:,11]
            for row, weight in zip(records, weights/weights.sum()*residence[ch,o]):
                offset = (row[5:8]+row[8:11])/2
                transport.setdefault((o,ch), []).append((offset[2], (row[3]+row[4])/2, weight,
                                                          _transport_xy(carrier, xy, offset)))
    # The spatial mark is independent of beam energy; time is expanded over
    # birth nodes (not replaced by a single averaged age for Cox boundaries).
    exy = np.stack([event[3] for group in transport.values() for event in group], axis=1)
    hx = [_transport_xy(carrier, xy, row[3:6]) for row in ends]
    all_heat, all_e, rows = [], [], []
    wdepth = np.diff(zedges)
    for ec in ENERGIES:
        rate = src.current(ec)/QE
        thermal = [{}, {}]
        ep = []
        source = 0.
        e_channels = np.zeros((2, 2))
        maximum_endpoint = 0.
        maxima = {'maximum_short_age_s': 0., 'maximum_e_age_s': 0.}
        for sp in range(2):
            tr = tracks[f'{ec:g}_{sp}']
            dots = model.source_totals(tr, sp, arm)[:4]
            birth = [_birth_rows(tr, sp, o, src) for o in range(4)]
            with np.load(path/f'source_E{ec:g}_nominal_retained_sp{sp}.npz') as archive:
                sound_J = float(archive['energy_J'].sum())
            parent_heat = ec*2000*QE-sound_J
            source += SPECIES[sp]*rate*parent_heat
            for row, distribution_xy in zip(ends, hx):
                o = int(row[0]); depth, age, weight = birth[o]
                ibin = _zbin(depth+row[5], zedges)
                release = age+row[2]
                maxima['maximum_short_age_s'] = max(maxima['maximum_short_age_s'], float(release.max()))
                maximum_endpoint = max(maximum_endpoint, float(_fold(depth+row[5], src.H).max()))
                # Clock projection with true correlated depth and birth age.
                fractional = release/LAG_STEP
                left = np.floor(fractional).astype(int)
                power = parent_heat*rate*SPECIES[sp]*row[6]*dots[o]/dots.sum()
                for lag in np.unique(np.r_[left, left+1]):
                    wt = weight*np.maximum(0., 1-abs(lag-fractional))
                    if not np.any(wt):
                        continue
                    dz = np.bincount(ibin, weights=wt, minlength=len(wdepth))
                    for b in range(2):
                        if lag not in thermal[b]:
                            thermal[b][lag] = np.zeros((*xy.shape[1:], len(wdepth)))
                        thermal[b][lag] += power*distribution_xy[b,...,None]*dz
            ratios = native.constitutive(arm, ec)['species_residence_ratio'][sp]
            number = 0
            for (o,ch), group in transport.items():
                depth, age, wt = birth[o]
                for offset, delay, residence_mark, _ in group:
                    release = age+delay
                    ibin = _zbin(depth+offset, zedges)
                    total = dots[o]*rate*SPECIES[sp]*residence_mark*ratios[ch]
                    weight = total*V0/(src.L/xy.shape[-1])**2*wt/wdepth[ibin]
                    ep.extend(zip(np.full(len(wt), number), ibin, release, weight))
                    e_channels[sp,int(ch==8)] += total
                    maxima['maximum_e_age_s'] = max(maxima['maximum_e_age_s'], float(release.max()))
                    number += 1
        e_array = np.asarray(ep)
        all_e.append(e_array)
        heat_defect = [float(sum(m.sum() for m in group.values())-source) for group in thermal]
        rows.append(dict(E_cm_keV=ec, parent_rate_s=rate, current_A=src.current(ec),
                         short_heat_W=source, short_source_map_defect_W=heat_defect,
                         deepest_short_endpoint_m=maximum_endpoint,
                         electron_residence_eV_population_by_parent_hot_then_thermal=e_channels.tolist(), **maxima))
        all_heat.append(thermal)
    return all_heat, exy, all_e, rows, segments


def _operators(xedges, zedges, k, cv, edge_h):
    dx = np.diff(xedges)[0]; dz = np.diff(zedges)
    lateral = k/dx*np.ones(len(xedges)-2)
    edge = 1/(1/edge_h+dx/(2*k))
    faces = 1/(1/(H_FACE+H_RAD)+dz[[0,-1]]/(2*k))
    axial = k/(.5*(dz[1:]+dz[:-1]))

    def modes(link, ends, cap):
        diagonal = np.r_[ends[0],link]+np.r_[link,ends[1]]
        _, vectors = eigh_tridiagonal(diagonal/cap, -link/np.sqrt(cap[:-1]*cap[1:]), lapack_driver='stemr')
        basis = vectors/np.sqrt(cap[:,None])
        # Re-evaluate lambda as a positive energy form, not the cancellation
        # of O(1e12) diagonal/off-diagonal contributions to the slow mode.
        rate = link@np.diff(basis, axis=0)**2 + ends@basis[[0,-1],:]**2
        return rate, vectors

    ax, qz = modes(axial, faces, cv*dz)
    lat, qx = modes(lateral, np.array([edge,edge]), np.full(len(xedges)-1, cv*dx))
    # Physical FV face flux uses conductance/area, not eigenvalue arithmetic.
    return lat[:,None,None]+lat[None,:,None]+ax[None,None,:], qx, qz, edge, faces, axial


def _solve_heat(groups, beam, clock, cell_mass, op, times, sampled):
    rates, qx, qz, *_ = op
    cap = np.sqrt(cell_mass)
    shape = rates.shape
    lags = np.array(sorted(groups))
    # FV operator in per-cell physical units: sqrt(M)^{-1} G sqrt(M)^{-1}.
    def project(p):
        return np.einsum('ia,jb,kc,ijk->abc', qx,qx,qz,p/cap,optimize=True)
    def recover(p):
        return np.einsum('ia,jb,kc,abc->ijk', qx,qx,qz,p,optimize=True)/cap
    modal = np.stack([project(groups[lag]) for lag in lags])
    stationary = recover(modal.sum(axis=0)/rates)
    cumulative = np.zeros((len(sampled),)+shape)
    state = np.zeros(shape)
    time_integral = np.zeros(shape)
    contraction = 1/(1+DT*rates)
    arrivals = np.diff(clock.primitive(times[:,None]-lags*LAG_STEP, beam), axis=0)
    input_j = 0.
    powers = np.array([groups[lag].sum() for lag in lags])
    j = 1
    for i, count in enumerate(arrivals, 1):
        state *= contraction
        state += np.einsum('l,lijk->ijk', count, modal, optimize=False)*contraction
        time_integral += DT*state
        input_j += float(count@powers)
        if j < len(sampled) and i == sampled[j]:
            cumulative[j] = recover(state)
            j += 1
    return stationary, cumulative, recover(time_integral), float(input_j), float(sum(powers))


def _split_temperature(history, depths):
    # Do not round a 10-K lateral background and a microkelvin depth contrast
    # together; gradients are evaluated on the separated numbers.
    column = (history@depths)/depths.sum()
    return column.astype('float32'), (history-column[...,None]).astype('float32')


def _cross_flux(temp, dx, dz, op):
    _, _, _, edge_g, face_g, _ = op
    edge = edge_g*dx*((temp[0,:,:]+temp[-1,:,:]+temp[:,0,:]+temp[:,-1,:]).sum(axis=0)@dz)
    faces = dx*dx*np.dot(face_g, temp[:,:,[0,-1]].sum(axis=(0,1)))
    return np.array([edge, faces*H_FACE/(H_FACE+H_RAD), faces*H_RAD/(H_FACE+H_RAD)])


def _front_budget(temp, dx, dz, op, n):
    """Outward flux through the illuminated layer's rim, front and lower face."""
    _, _, _, edge_g, face_g, axial = op
    slab = temp[...,:n]
    side = edge_g*dx*((slab[0]+slab[-1]+slab[:,0]+slab[:,-1]).sum(axis=0)@dz[:n])
    front = dx*dx*face_g[0]*slab[...,0].sum()
    interior = dx*dx*axial[n-1]*(temp[...,n-1]-temp[...,n]).sum()
    return [float(side),float(front),float(interior)]


def _native_payload(path, src, native, segments):
    p = src.inputs()['lattice']
    lattice_bonds, lattice_j = src.alias_bonds(p)
    old_onsite = p['Kii_eV_A2']+(2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4'])*p['zero_point_axis_variance_A2']
    material = native.Material(16)
    q = native.primitive_basis(p)
    q_cut = np.min(np.linalg.norm(q, axis=1)[np.linalg.norm(q, axis=1)>1.])
    launch_q = np.unique(segments[segments[:,2]==4, 14:17], axis=0)
    data = dict(bonds=lattice_bonds, bonds_j=lattice_j, onsite=old_onsite, q_cut=q_cut,
                native_reference=material.reference, electronic_abscissa=material.u, electronic_deltaK=material.dk,
                first_shell=np.asarray(p['shells'][0]['periodic_alias_offsets_A'])*1e-10,
                optical_mass=p['mass_eVc2']*QE/299792458.**2, linear_strain_slope=p['collective_strain_slope_eV_A2'],
                covariance_coefficient=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4'], alpha=src.ALPHA,
                optical_tau=src.NOMINAL['opt'], extra_dephasing_tau=src.NOMINAL['dephase'], launch_q=launch_q)
    for ec in ENERGIES:
        for allocation in ALLOCATIONS:
            for sp in (0,1):
                with np.load(path/f'native_E{ec:g}_nominal_retained_{allocation}_sp{sp}_m16.npz') as archive:
                    prefix = f'parent_{ec:g}_{allocation}_{sp+1}_'
                    for k in ('time_s','maximum_U_eV_formula','hot_e_eV','thermal_e_excess_eV'):
                        data[prefix+k] = archive[k]
                    if ec == 2.5:
                        for k in ('map_time_s','map_U','positions_r_z_nm','focus_weight'):
                            data[prefix+k] = archive[k]
    return data


def build(out, path4):
    """Build the resolved medium and compact diagnostics without writing Pass 4."""
    out = Path(out); path4 = Path(path4).resolve()
    bindings = _bind_inputs(path4)
    out.mkdir(parents=True, exist_ok=True)
    with _pass4(path4) as mod, np.load(path4/'calendar.npz') as calendar:
        clock = Clock(calendar)
        src = mod['sources']; field = mod['fields']; native = mod['native']; native.load()
        tracks = json.loads((path4/'frozen/tracks.json').read_text())
        zedges = np.r_[np.linspace(0.,.25e-6,51), np.geomspace(.35e-6,src.H,18)]
        dz = np.diff(zedges); xedges = np.linspace(-src.L/2,src.L/2,field.N+1)
        dx = float(np.diff(xedges)[0]); xy = np.stack([src.beam_cells(xedges, b) for b in BEAMS[:2]])
        front_count = int(np.searchsorted(zedges, .25e-6))
        shape = (field.N,field.N,len(dz)); mass = np.broadcast_to(src.CV*dx**2*dz,shape)
        heat, exy, e_events, source_rows, segments = _carrier_fields(mod, path4, src.NOMINAL, tracks, xy, zedges)
        EndpointOperator = _refined_acoustic(mod['acoustic'], field, mod['launch'], zedges)
        times = np.linspace(0.,clock.seconds,int(round(clock.seconds/DT))+1)
        sample_indices = np.unique(np.r_[0,1,2,4,8,12,20,np.arange(25,1990,25),1990,1995,1999,2000])
        ops = {m:_operators(xedges,zedges,src.KTH,src.CV,h) for m,h in MOUNTS.items()}
        thermal_t2 = np.empty((len(ENERGIES),len(BEAMS),len(MOUNTS),*shape))
        steady = np.empty_like(thermal_t2)
        archive = {k:calendar[k] for k in ('edges_s','burst_on','phase','normalization')}
        archive.update(xedges=xedges, zedges=zedges, cross_section=shape, mounts=np.array(list(MOUNTS)),
                       edge_h=np.array(list(MOUNTS.values())), k_total=src.KTH, time_nodes=times[sample_indices],
                       electron_xy=exy, **_native_payload(path4,src,native,segments))
        conservation = []
        for i, ec in enumerate(ENERGIES):
            archive[f'e_events_{i}'] = e_events[i]
            rate = src.current(ec)/QE
            for shape_id in range(2):
                acoustic = EndpointOperator(src.NOMINAL,ec,BEAMS[shape_id],tracks)
                merged = {j:p.copy() for j,p in heat[i][shape_id].items()}
                for j, p in acoustic.groups.items():
                    if j in merged:
                        merged[j] += p
                    else:
                        merged[j] = p
                sound_exit_W = acoustic._exit.total*rate
                sound_delays = acoustic.maximum_delay
                acoustic_heat_W = float(sum(p.sum() for p in acoustic.groups.values()))
                archived_bath = acoustic.bath_e*V0/(QE*dx*dx*dz)
                aged = float(acoustic.bath_e_age/acoustic.bath_e.sum()) if acoustic.bath_e.sum() else 0.
                archive[f'bath_{i}_{shape_id}'] = archived_bath
                archive[f'bath_age_{i}_{shape_id}'] = aged
                del acoustic
                front_source = float(sum(p[...,:front_count].sum() for p in merged.values()))
                for beam_id in ((0,) if shape_id == 0 else (1,2)):
                    beam = BEAMS[beam_id]
                    for mount_id,(mount,op) in enumerate(ops.items()):
                        stat, history, integral, deposited_J, source_W = _solve_heat(merged,beam,clock,mass,op,times,sample_indices)
                        final = history[-1]
                        thermal_t2[i,beam_id,mount_id] = final
                        steady[i,beam_id,mount_id] = stat
                        prefix=f'T_{i}_{beam_id}_{mount_id}_'
                        archive[prefix+'col'], depth = _split_temperature(history, dz)
                        archive[prefix+'depth'] = depth
                        stored = float((final*mass).sum())
                        exits_J = _cross_flux(integral,dx,dz,op); exits_W = _cross_flux(stat,dx,dz,op)
                        inc_W = src.current(ec)*2000*ec
                        front_exit = _front_budget(stat,dx,dz,op,front_count)
                        front_arrival = float(sum(float(clock.primitive(clock.seconds-j*LAG_STEP,beam))*p[...,:front_count].sum()
                                                  for j,p in merged.items()))
                        front_stored = float((final[...,:front_count]*mass[...,:front_count]).sum())
                        front_integrated = _front_budget(integral,dx,dz,op,front_count)
                        conservation.append(dict(E_cm_keV=ec,beam=beam,mount=mount,incident_W=inc_W,
                            short_heat_W=source_rows[i]['short_heat_W'], acoustic_thermalized_W=acoustic_heat_W,
                            outgoing_acoustic_W=sound_exit_W, max_carrier_or_acoustic_delay_s=max(sound_delays,source_rows[i]['maximum_short_age_s']),
                            inherited_acoustic_bath_Ue_max_eV_formula=float(archived_bath.max()), total_heat_W=source_W,
                            input_partition_defect_W=inc_W-source_W-sound_exit_W, thermal_sink_W=exits_W.tolist(),
                            steady_balance_W=float(source_W-exits_W.sum()), deposited_to_heat_2s_J=deposited_J,
                            front_250nm_source_W=front_source, front_250nm_exits_W=front_exit,
                            front_250nm_steady_balance_W=front_source-sum(front_exit),
                            front_250nm_transient_balance_J=front_arrival-front_stored-sum(front_integrated),
                            thermal_sink_2s_J=exits_J.tolist(),stored_heat_2s_J=stored,
                            transient_balance_J=float(deposited_J-exits_J.sum()-stored),
                            maximum_T_2s_K=T_AMBIENT+float(final.max()), maximum_T_steady_K=T_AMBIENT+float(stat.max())))
        archive.update(steady=steady, thermal_t2=thermal_t2)
        bundle = out/'medium.npz'
        np.savez_compressed(bundle, **archive)
        del archive
    metadata = dict(baseline_commit=BASELINE, consumed_sha256=bindings,
        geometry=dict(foil_side_m=src.L, thickness_m=src.H, beam_area_m2=1e-4, n_lateral=field.N,
                      front_resolved_cell_depth_m=float(dz[0]), front_resolved_depth_m=.25e-6),
        transport=dict(k_total_W_m_K=src.KTH,k_e_Pd_reference_W_m_K=2.44e-8*300/1.08e-7,
                       k_e_note='pure Pd, 300 K, Wiedemann–Franz with inherited 10.80 µΩ cm; not measured loaded PdD and not added to total',
                       rho_kg_m3=src.RHO,cv_J_m3_K=src.CV, bath_K=T_AMBIENT,
                       edge_mount_W_m2_K=MOUNTS, face_W_m2_K=H_FACE,radiation_linearized_W_m2_K=H_RAD,
                       electronic_formula_volume_m3=V0),
        conventions=dict(heat='whole-foil conservative Cartesian FV, shared faces/edges, path thermalization endpoints; signed outgoing acoustic power is not thermalized',
                         finite_clock='actual normalized archived Cox schedule integrated in 1 ms conditioned-mean backward Euler heat steps; endpoint lags projected to 50 µs positive moment-preserving bins; off-grid heat samples interpolate retained thermal snapshots, no claim of sub-ms temperature',
                         mean_clock='unit-rate steady intensity mean, each beam and named mount',
                         electronic='retarded short-path normalized hot + thermal excess, P4 birth GL2/molecule clock with specular folds; tiny distributed acoustic bath evaluated at its residence-weighted mean delay; not electron temperature',
                         conditional='single-parent archived native age histories for six energies, both allocations; 2.5-keV archived space/age nodes only, not added to mean',
                         mechanical='dynamic linear strain and its gradient supplied separately by the signed-wave consumer; no independent static clamp strain law assumed',
                         below_resolution='no nonlinear pressure from the ensemble mean assumed; no sub-5-nm steady thermal resolution or sub-ms conditioned temperature fluctuations inferred'),
        sources=source_rows)
    _json(out/'medium_parameters.json',metadata)
    _json(out/'medium_conservation.json',dict(thermal_exit_order=['mount','face conductance','radiation'],
          front_250nm_exit_order=['rim mount','illuminated face conduction plus radiation','into deeper foil'],
          cases=conservation))
    _readout(Medium(bundle),out)
    return bundle


def _json(path, record):
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


def _interpolate(values, axes, points):
    """Trilinear value and derivative of the same interpolant, in SI."""
    indices, weights, widths = [], [], []
    for axis, coordinate in zip(axes, points.T):
        i = np.clip(np.searchsorted(axis,coordinate,side='right')-1,0,len(axis)-2)
        width = axis[i+1]-axis[i]
        f = (coordinate-axis[i])/width
        indices.append(i); weights.append((1-f,f)); widths.append(width)
    i,j,k = indices
    x,y,z = weights
    corners = {(a,b,c):values[i+a,j+b,k+c] for a in (0,1) for b in (0,1) for c in (0,1)}
    value = sum(v*x[a]*y[b]*z[c] for (a,b,c),v in corners.items())
    gx = sum((corners[1,b,c]-corners[0,b,c])*y[b]*z[c] for b in (0,1) for c in (0,1))/widths[0]
    gy = sum((corners[a,1,c]-corners[a,0,c])*x[a]*z[c] for a in (0,1) for c in (0,1))/widths[1]
    gz = sum((corners[a,b,1]-corners[a,b,0])*x[a]*y[b] for a in (0,1) for b in (0,1))/widths[2]
    return value, np.column_stack((gx,gy,gz))


class Medium:
    """Read-only SI field consumer; `.sample` never adds a tagged parent's excitation."""
    def __init__(self, bundle, *, native_projection=None):
        self.native_projection = native_projection
        path = Path(bundle)
        if path.is_dir():
            path = path/'medium.npz'
        self.data = np.load(path)
        self.k_total = float(self.data['k_total'])
        self.time = self.data['time_nodes']; self.clock = Clock(self.data)
        self.xedges = self.data['xedges']; self.zedges = self.data['zedges']
        self.x = (self.xedges[:-1]+self.xedges[1:])/2
        self.z = (self.zedges[:-1]+self.zedges[1:])/2
        self.dz = np.diff(self.zedges); self.dx = float(np.diff(self.xedges)[0])
        self.mounts = self.data['mounts'].tolist(); self.h_edge = self.data['edge_h']
        self.midaxes = (np.r_[self.xedges[0],self.x,self.xedges[-1]],)*2
        self.thermal_z = np.empty(2*len(self.z)+1)
        self.thermal_z[::2] = self.zedges; self.thermal_z[1::2] = self.z
        self.e_z = np.r_[self.zedges[0],self.z,self.zedges[-1]]
        self.steady = self.data['steady']; self.two_seconds = self.data['thermal_t2']
        self.e_xy = self.data['electron_xy']
        self.e_cache = {}
        self.t_history_key = None
        self.t_history = None
        self.temperature_cache = {}
        self.branches = {k:self.data[k] for k in ('bonds','bonds_j','onsite','q_cut','native_reference',
            'electronic_abscissa','electronic_deltaK','first_shell','optical_mass',
            'linear_strain_slope','covariance_coefficient','alpha','optical_tau','extra_dephasing_tau','launch_q')}

    def close(self):
        self.data.close()

    def _thermal_nodes(self, cell, m):
        # Boundary node is at the resistance-weighted temperature on the
        # physical surface. A shared depth face has a single temperature and
        # identical one-sided normal flux from the two half-cells.
        z = np.empty(cell.shape[:-1]+(len(self.thermal_z),))
        z[...,1::2] = cell
        b = self.dz[1:]/(self.dz[1:]+self.dz[:-1])
        z[...,2:-1:2] = b*cell[...,:-1]+(1-b)*cell[...,1:]
        z[...,0] = cell[...,0]/(1+(H_FACE+H_RAD)*self.dz[0]/(2*self.k_total))
        z[...,-1] = cell[...,-1]/(1+(H_FACE+H_RAD)*self.dz[-1]/(2*self.k_total))
        side = 1/(1+self.h_edge[m]*self.dx/(2*self.k_total))
        z = np.pad(z,((1,1),(1,1),(0,0)),mode='edge')
        z[[0,-1],:,:] *= side; z[:,[0,-1],:] *= side
        return z

    def _history_array(self, name):
        if name in self.data.files:
            return self.data[name]
        import lzma
        payload=self.data[name+'_lzma'].tobytes()
        shape=tuple(self.data[name+'_shape'])
        shuffled=np.frombuffer(lzma.decompress(payload),dtype=np.uint8).reshape(4,-1).T.copy()
        delta=shuffled.reshape(-1).view(np.uint32).reshape(shape)
        return np.cumsum(delta,axis=-1,dtype=np.uint32).view(np.float32)

    def _temperature(self, i, b, m, t, point):
        value = np.empty(len(t)); grad = np.empty((len(t),3))
        axes = (*self.midaxes,self.thermal_z)
        def stationary_or_last(mask, source, key):
            if not np.any(mask):
                return
            cached = self.temperature_cache.get(key)
            if cached is None:
                cached = self._thermal_nodes(source[i,b,m],m)
                self.temperature_cache[key] = cached
            value[mask], grad[mask] = _interpolate(cached,axes,point[mask])

        stationary_or_last(np.isposinf(t), self.steady, (i,b,m,True))
        stationary_or_last(t == self.clock.seconds, self.two_seconds, (i,b,m,False))
        other = np.flatnonzero(np.isfinite(t) & (t != self.clock.seconds))
        if len(other):
            if self.t_history_key != (i,b,m):
                prefix = f'T_{i}_{b}_{m}_'
                self.t_history = (self.data[prefix+'col'],self._history_array(prefix+'depth'))
                self.t_history_key = (i,b,m)
            col, depth = self.t_history
            lower = np.clip(np.searchsorted(self.time,t[other],side='right')-1,0,len(self.time)-2)
            for j in np.unique(lower):
                selected = other[lower==j]
                mix = (t[selected]-self.time[j])/(self.time[j+1]-self.time[j])
                aa,gg = _interpolate(self._thermal_nodes(col[j,...,None].astype(float)+depth[j],m), axes, point[selected])
                bb,hh = _interpolate(self._thermal_nodes(col[j+1,...,None].astype(float)+depth[j+1],m), axes, point[selected])
                value[selected] = (1-mix)*aa + mix*bb
                grad[selected] = (1-mix[:,None])*gg + mix[:,None]*hh
        return T_AMBIENT+value, grad

    @staticmethod
    def _mirror_nodes(cells):
        return np.pad(cells,((1,1),(1,1),(1,1)),mode='edge')

    def _electronic_grids(self, i, b):
        shape = 0 if b == 0 else 1
        key = (i, shape)
        if key not in self.e_cache:
            events = self.data[f'e_events_{i}']
            ix,z = events[:,:2].astype(int).T
            age, w = events[:,2:].T
            xy = self.e_xy[shape]
            depth_shape = (xy.shape[0], len(self.z))
            address = ix*len(self.z)+z
            scale = np.bincount(address,weights=w,minlength=np.prod(depth_shape)).reshape(depth_shape)
            ripple = np.zeros(depth_shape,complex)
            np.add.at(ripple,(ix,z),w*np.exp(-1j*OMEGA_COX*age))
            mean = np.einsum('kxy,kz->xyz',xy,scale,optimize=True)
            harmonic = np.einsum('kxy,kz->xyz',xy,ripple,optimize=True)
            bath = self._mirror_nodes(self.data[f'bath_{i}_{shape}'])
            bath_age = float(self.data[f'bath_age_{i}_{shape}'])
            self.e_cache[key] = (self._mirror_nodes(mean),self._mirror_nodes(harmonic),bath,bath_age,
                                 age, w, address, depth_shape, float(age.max()), xy)
        return self.e_cache[key]

    def _electronic(self, i, b, t, point):
        mean, harmonic, bath, bath_age, age, w, address, shape, maximum, xy = self._electronic_grids(i,b)
        axes = (*self.midaxes,self.e_z)
        value, grad = _interpolate(mean,axes,point)
        bath_v, bath_g = _interpolate(bath,axes,point)
        bath_scale = np.ones(len(t))
        finite = np.flatnonzero(np.isfinite(t))
        if len(finite):
            at = t[finite]
            exact = at<maximum
            if b == 2:
                state = self.clock.index(at)
                exact |= (at-self.clock.edges[state]) < maximum
                c = (1+9*self.clock.state[state])/self.clock.norm
                r = .7*np.exp(1j*(OMEGA_COX*at+self.clock.phase))/self.clock.norm
                dv,dg = _interpolate(harmonic,axes,point[finite])
                value[finite] = value[finite]*c+(r*dv).real
                grad[finite] = grad[finite]*c[:,None]+(r[:,None]*dg).real
                bath_scale[finite] = self.clock.value(at-bath_age)
            else:
                bath_scale[finite] = at>=bath_age
            # Only a ps-wide strip at startup/switches needs the full selected
            # delayed population. Away from a switch, the harmonic moment is exact.
            for now in np.unique(at[exact]):
                selected = finite[exact & (at==now)]
                retarded = now-age
                c = self.clock.value(retarded) if b==2 else (retarded>=0)
                depth = np.bincount(address,weights=w*c,minlength=np.prod(shape)).reshape(shape)
                grid = np.einsum('kxy,kz->xyz',xy,depth,optimize=True)
                value[selected], grad[selected] = _interpolate(self._mirror_nodes(grid),axes,point[selected])
        return value+bath_v*bath_scale, grad+bath_g*bath_scale[:,None]

    def sample(self, ec, beam, xyz_m, mount='nominal', time_s=2):
        """Return (T, Ue, ∇T, ∇Ue) for xyz_m[...,3] and broadcastable time_s.

        x/y origin is foil centre; z=0 is the illuminated face. `np.inf`
        uses the unit mean-rate steady heat and electron field. Finite time is
        the two-second inherited source aperture, not a repeated Cox history.
        """
        i,b,m = ENERGIES.index(ec), BEAMS.index(beam), self.mounts.index(mount)
        point = np.asarray(xyz_m, dtype=float)
        if point.ndim == 0 or point.shape[-1] != 3 or not np.isfinite(point).all():
            raise ValueError('xyz_m must have shape (...,3) and finite SI coordinates')
        time = np.asarray(time_s,dtype=float)
        shape = np.broadcast_shapes(point.shape[:-1],time.shape)
        at = np.broadcast_to(time,shape).ravel()
        if not np.all((0<=at)&((at<=self.clock.seconds)|np.isposinf(at))):
            raise ValueError('Finite medium time must lie in accepted parent aperture [0,2] s')
        point = np.broadcast_to(point,(*shape,3)).reshape(-1,3)
        bounds = np.array([self.xedges[[0,-1]],self.xedges[[0,-1]],self.zedges[[0,-1]]]).T
        if np.any((point<bounds[0])|(point>bounds[1])):
            raise ValueError('xyz_m outside original foil')
        temperature, grad_T = self._temperature(i,b,m,at,point)
        electronic, grad_U = self._electronic(i,b,at,point)
        if not shape:
            return float(temperature[0]),float(electronic[0]),grad_T[0],grad_U[0]
        return (temperature.reshape(shape),electronic.reshape(shape),
                grad_T.reshape((*shape,3)),grad_U.reshape((*shape,3)))

    def _hamiltonian(self, q, temperature, electronic, grad_T, grad_U, strain, grad_strain, allocation, frozen):
        p = self.branches
        phase = q@p['bonds'].T
        bare = p['onsite']-np.cos(phase)@p['bonds_j']
        dq = np.einsum('...b,b,bi->...i',np.sin(phase),p['bonds_j'],p['bonds'])
        if frozen:
            temperature = np.full_like(temperature,T_AMBIENT)
        else:
            gamma = p['onsite']-p['bonds_j'].sum()
            eta = gamma/p['native_reference']-1
            shell = np.maximum(0.,1-(q*q).sum(axis=-1)/p['q_cut']**2)
            taper = 1+eta*shell**3
            dq = dq/taper[...,None]+(6*eta*bare*shell**2/taper**2/p['q_cut']**2)[...,None]*q
            bare = bare/taper
            if self.native_projection is not None:
                de, slope, electronic_dq = self.native_projection.project(q,electronic,allocation)
                dq = dq+electronic_dq
            else:
                u, d = p['electronic_abscissa'], p['electronic_deltaK']
                if np.any(electronic < -1e-13) or np.any(electronic>u[-1]*(1+1e-12)):
                    raise ValueError('Electronic field outside supported G16 abscissa')
                de = np.interp(electronic,u,d)
                index = np.clip(np.searchsorted(u,electronic,side='right')-1,0,len(u)-2)
                slope = np.diff(d)[index]/np.diff(u)[index]
                if allocation == 'coupling':
                    argument = q@p['first_shell'].T
                    weight = np.cos(argument).mean(axis=-1)
                    dq = dq-de[...,None]*np.einsum('...b,bi->...i',np.sin(argument),p['first_shell'])/len(p['first_shell'])
                    de = de*weight
                    slope = slope*weight
            g = p['linear_strain_slope']
            bare = bare+de+g*(strain+p['alpha']*(temperature-T_AMBIENT))
        if np.any(bare<=0):
            raise ValueError('Optical branch outside retained stable support')
        C = p['covariance_coefficient']
        cov = 2*KB*temperature/(bare+np.sqrt(bare*bare+4*C*KB*temperature))
        dD = bare+2*C*cov
        evolved = bare+C*cov
        omega = np.sqrt(evolved*(QE/1e-20)/p['optical_mass'])
        curvature_to_omega = omega/(2*evolved)
        vg = (curvature_to_omega*evolved/dD)[...,None]*dq
        if frozen:
            return omega, vg, np.zeros_like(vg)
        direct = (g*p['alpha']*evolved/dD+C*KB/dD)[...,None]*grad_T
        field = (evolved/dD)[...,None]*(slope[...,None]*grad_U+g*grad_strain)
        return omega, vg, -curvature_to_omega[...,None]*(direct+field)

    def optical(self, ec, beam, xyz_m, q_m_inverse, mount='nominal', time_s=2, allocation='onsite',
                strain=0., grad_strain=0., frozen=False):
        """Angular Ω, ∂qΩ and −∂xΩ in rad/s, m/s and m⁻¹/s.

        `strain` is linear strain (trace/3), `grad_strain` is per metre.
        An optional wave actuator may set both; frozen removes all gradients.
        """
        allowed = ALLOCATIONS if self.native_projection is None else self.native_projection.allocations
        if allocation not in allowed:
            raise ValueError(f'Unknown optical allocation {allocation}')
        fields = self.sample(ec,beam,xyz_m,mount,time_s)
        return self._hamiltonian(np.asarray(q_m_inverse,float),*(np.asarray(a) for a in fields),
                                 np.asarray(strain,float),np.asarray(grad_strain,float),allocation,frozen)

    def conditional_electron(self, ec, species, age_s, allocation='onsite'):
        """Archived single-parent peak and integrated eV, NOT local mean U or a temperature."""
        if (ec not in ENERGIES) or (species not in (1,2)) or (allocation not in ALLOCATIONS):
            raise ValueError('Expected accepted energy, species 1 (D+) or 2 (D₂+), and allocation')
        k = f'parent_{ec:g}_{allocation}_{species}_'
        at = np.asarray(age_s)
        times = self.data[k+'time_s']
        if np.any(~np.isfinite(at) | (at<0) | (at>times[-1])):
            raise ValueError('Single-parent age outside archived support; no tail extrapolation')
        t = np.r_[0.,times]
        return {name: np.interp(at,t,np.r_[0.,self.data[k+name]]) for name in
                ('maximum_U_eV_formula','hot_e_eV','thermal_e_excess_eV')}

    def conditional_map(self, ec, species, allocation='onsite'):
        """2.5 keV native, tagged-parent exact mesh nodes and ages, not interpolated."""
        if ec!=2.5 or species not in (1,2) or allocation not in ALLOCATIONS:
            raise ValueError('Conditional spatial snapshots are inherited only at Ecm=2.5 keV')
        prefix = f'parent_{ec:g}_{allocation}_{species}_'
        return {k:self.data[prefix+k] for k in ('map_time_s','map_U','positions_r_z_nm','focus_weight')}


def _csv(path, rows):
    with path.open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def _trace(medium, p0, q0, beam, time_s, allocation, frozen):
    """Characteristic over one and five energy lifetimes; all axes are in SI."""
    tau = float(medium.branches['optical_tau'])
    r, q = p0.copy(), q0.copy()
    rs, qs, freqs = [r.copy()], [q.copy()], []
    size = np.array([2*medium.xedges[-1],2*medium.xedges[-1],medium.zedges[-1]])
    lower = np.array([medium.xedges[0],medium.xedges[0],0.])
    reflections = 0
    def evaluate(x,k):
        return medium.optical(2.5,beam,x,k,time_s=time_s,allocation=allocation,frozen=frozen)
    def reflect(x):
        offset = x-lower
        sense = np.where(np.floor(offset/size)%2, -1, 1)
        return lower+_fold(offset,size), sense
    step = tau/16
    for _ in range(80):
        w,v,f = evaluate(r,q)
        pm, sense = reflect(r+.5*step*v)
        _,v_mid,f_mid = evaluate(pm,sense*(q+.5*step*f))
        r, end_sense = reflect(r+step*sense*v_mid)
        q = end_sense*(q+step*sense*f_mid)
        reflections += int(np.count_nonzero(end_sense<0))
        rs.append(r.copy()); qs.append(q.copy()); freqs.append(float(w))
    freqs.append(float(evaluate(r,q)[0]))
    return np.array(rs),np.array(qs),np.array(freqs),reflections


def _readout(medium, out):
    profiles, rays = [], []
    depths = np.r_[0.,5e-9,10e-9,25e-9,50e-9,100e-9,150e-9,250e-9,1e-6,125e-6,250e-6]
    for ec in ENERGIES:
        for beam in BEAMS:
            for label,at in (('two_seconds',2.),('mean_steady',np.inf)):
                for spot, xy in (('center',(0.,0.)),('primary_hotspot',(.0025,-.0025)),('outside_beam',(.0075,0.))):
                    p = np.column_stack((np.broadcast_to(xy,(len(depths),2)),depths))
                    T,U,gT,gU = medium.sample(ec,beam,p,time_s=at)
                    for z,t,u,gt,gu in zip(depths,T,U,gT,gU):
                        profiles.append(dict(E_cm_keV=ec,beam=beam,mount='nominal',snapshot=label,spot=spot,depth_m=z,
                            T_K=t,Ue_eV_formula=u,dTdx_K_m=gt[0],dTdy_K_m=gt[1],dTdz_K_m=gt[2],
                            normal_heat_flux_W_m2=-medium.k_total*gt[2],
                            dUdx_eV_formula_m=gu[0],dUdy_eV_formula_m=gu[1],dUdz_eV_formula_m=gu[2]))
    _csv(out/'medium_profiles.csv',profiles)
    qbank = medium.branches['launch_q']
    tau = float(medium.branches['optical_tau'])
    qlen = np.linalg.norm(qbank,axis=1)
    _,v0,_ = medium.optical(2.5,'hot_bursts',np.array([.0025,-.0025,50e-9]),qbank)
    r0s = {'depth':np.array([.0025,-.0025,50e-9]),'lateral':np.array([.0035,-.0025,50e-9])}
    for direction, axis in (('depth',2),('lateral',0)):
        candidate = np.flatnonzero((abs(v0[:,axis])>1e-6)&(qlen>0.))
        selected = np.unique([candidate[np.argmin(qlen[candidate])],
                              candidate[np.argmax(abs(v0[candidate,axis]))],
                              candidate[np.argmax(abs(qbank[candidate,axis]))]])
        r0 = r0s[direction]
        scan = np.broadcast_to(r0,(67,3)).copy()
        scan[:,axis] = depths[0]+np.linspace(0.,250e-9,67) if axis==2 else np.linspace(-.005,.005,67)
        for i in selected:
            q = qbank[i].copy()
            if v0[i,axis] < 0.:
                q *= -1
            for beam in BEAMS:
                for state, at in (('two_seconds',2.),('mean_steady',np.inf)):
                    for a in ALLOCATIONS:
                        T,U,grad_T,grad_U = (np.asarray(x) for x in medium.sample(2.5,beam,r0,time_s=at))
                        h,c,g = medium._hamiltonian(q,T,U,grad_T,grad_U,0.,0.,a,False)
                        _,_,g_T = medium._hamiltonian(q,T,U,grad_T,np.zeros(3),0.,0.,a,False)
                        _,_,g_U = medium._hamiltonian(q,T,U,np.zeros(3),grad_U,0.,0.,a,False)
                        cold = medium.optical(2.5,beam,r0,q,time_s=at,allocation=a,frozen=True)[0]
                        w,_,_ = medium.optical(2.5,beam,scan,q,time_s=at,allocation=a)
                        q0n = q.copy(); q0n[axis] = 0.
                        turn_w, _, _ = medium.optical(2.5,beam,scan,q0n,time_s=at,allocation=a)
                        start_turn = float(medium.optical(2.5,beam,r0,q0n,time_s=at,allocation=a)[0])
                        # Necessary direction-branch candidate, not a claim that
                        # a lifetime-limited characteristic reaches the contour.
                        scan_ahead = scan[:,axis]>=r0[axis]
                        threshold = (h-start_turn)*(h-turn_w[scan_ahead])
                        candidate_turn = bool(np.any(threshold<0.))
                        path, qpath, ws, reflections = _trace(medium,r0,q,beam,at,a,False)
                        f_path, _, _, frozen_reflections = _trace(medium,r0,q,beam,at,a,True)
                        vg0 = medium.optical(2.5,beam,r0,q,time_s=at,allocation=a,frozen=True)[1]
                        _, vg_path, _ = medium.optical(2.5,beam,path,qpath,time_s=at,allocation=a)
                        wave = 2*np.pi/qlen[i]
                        rays.append(dict(E_cm_keV=2.5,beam=beam,snapshot=state,allocation=a,direction=direction,
                            initial_qx_m_inverse=q[0],initial_qy_m_inverse=q[1],initial_qz_m_inverse=q[2],
                            frozen_frequency_Hz=float(cold/(2*np.pi)),evolved_frequency_Hz=float(h/(2*np.pi)),
                            scan_frequency_span_Hz=float(np.ptp(w)/(2*np.pi)),
                            local_signed_kdot_normal_m_inv_s=float(g[axis]),
                            kdot_thermal_normal_m_inv_s=float(g_T[axis]),
                            kdot_electronic_normal_m_inv_s=float(g_U[axis]),
                            local_evolved_vg_normal_m_s=float(c[axis]),local_frozen_vg_normal_m_s=float(vg0[axis]),
                            one_lifetime_travel_m=float(np.linalg.norm(path[16]-r0)),
                            evolved_minus_frozen_one_lifetime_normal_m=float(path[16,axis]-f_path[16,axis]),
                            five_lifetime_relative_dq=float(np.linalg.norm(qpath[-1]-q)/qlen[i]),
                            five_lifetime_frequency_fraction_residual=float(np.ptp(ws)/h),
                            energy_lifetime_s=tau,amplitude_decay_lifetime_s=2*tau,
                            pure_dephasing_s=float(medium.branches['extra_dephasing_tau']),
                            physical_boundary_reflections=reflections,frozen_boundary_reflections=frozen_reflections,
                            packet_width_q_m_inverse=5e7,wavelength_m=float(wave),
                            local_WKB_lambda_times_log_gradient=float(wave*np.linalg.norm(g)/h),
                            local_1D_zero_normal_q_barrier_candidate=candidate_turn,
                            normal_vg_reversed_in_five_lifetimes=bool(not reflections and np.any(vg_path[:,axis]*c[axis]<0))))
    _csv(out/'medium_refraction.csv',rays)
    shot = [dict(species=sp, allocation=a, peak_eV_formula=float(medium.data[f'parent_2.5_{a}_{sp}_maximum_U_eV_formula'].max()))
            for a in ALLOCATIONS for sp in (1,2)]
    nominal_hot = [p for p in profiles if (p['E_cm_keV'],p['beam'],p['spot'])==(2.5,'hot_bursts','primary_hotspot')]
    d = {label:{p['depth_m']:p for p in nominal_hot if p['snapshot']==label} for label in ('two_seconds','mean_steady')}
    summary = {label:dict(surface_T_K=nodes[0.]['T_K'], delta_T_0_to_250nm_K=nodes[250e-9]['T_K']-nodes[0.]['T_K'],
                          front_flux_into_foil_W_m2=nodes[0.]['normal_heat_flux_W_m2'],
                          passing_flux_at_250nm_W_m2=nodes[250e-9]['normal_heat_flux_W_m2'],
                          peak_profile_mean_Ue_eV_formula=max(p['Ue_eV_formula'] for p in nodes.values()),
                          peak_profile_abs_electron_depth_gradient_eV_formula_m=max(abs(p['dUdz_eV_formula_m']) for p in nodes.values()))
               for label,nodes in d.items()}
    _json(out/'medium_readout.json',dict(nominal_hot_Cox_2p5keV=summary,tagged_parent_peaks=shot,
        maximum_observed_refraction_relative_dq_5tau=max(r['five_lifetime_relative_dq'] for r in rays),
        maximum_ray_integration_frequency_fraction_residual=max(r['five_lifetime_frequency_fraction_residual'] for r in rays),
        local_1D_zero_normal_q_barrier_candidates=sum(r['local_1D_zero_normal_q_barrier_candidate'] for r in rays),
        observed_normal_group_reversals=sum(r['normal_vg_reversed_in_five_lifetimes'] for r in rays),
        ray_scope='genuine retained path q marks, nominal two-second and mean static fields; fixed linear strain=0, not a statement about moving acoustic strain or conditional cloud; 1D slice criterion is heuristic outside locally axial propagation and does not establish trapping',
        frozen_control='unbridged frozen G3 q curvature with 300-K equilibrium covariance, no spatial gradients or carrier work; differences in speed include native bridge, whereas force and scan span isolate spatial variation',
        limitations='5-nm endpoint histogram and 17×17 xy finite volumes; no subcell steady electronic structure or inferred temperature; dynamic strain is a separate actuator. The native q projection is not measured and G16 low-U curvature sign is mesh-sensitive (16³ retained here). No crossing probability inferred.'))
    medium.close()
