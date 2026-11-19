"""Stream compact, matched observables from the actual corrected lattice engine.

No upstream file changes; one process runs six cases serially. Arrays retain
joint same-time nearest-neighbour information and sampled target velocities.
"""
from pathlib import Path
import hashlib
import json
import os
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT/'g3_lattice/rerun/part2'
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ENGINE))
from lattice import run_case, CELL, NN, HC, HBAR, KB


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Observer:
    def __init__(self, p, case):
        self.p, self.case = p, case
        eref = np.sqrt(HC*HC*p['Kii_eV_A2']/p['mass_eVc2'])
        self.length = HC/np.sqrt(p['mass_eVc2']*eref)
        self.dt_ps = case['dt']*HBAR/eref*1e12
        A = p['K4_100_eV_A4']; B = (3*p['K4_111_eV_A4']-A)/2
        self.kgamma = p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']-sum(s['J_sum_eV_A2'] for s in p['shells'])
        self.omega = np.sqrt(HC*HC*self.kgamma/p['mass_eVc2'])/HBAR*1e-12
        self.mass_eVs2_A2 = p['mass_eVc2']/(2.99792458e18)**2
        self.gamma_s = 2/p['reference_amplitude_decay_ps']*1e12
        self.count = 0
        self.old = []
        self.velocity_samples = []
        self.velocity_blocks = []
        self.times = []
        self.v2_blocks = np.zeros(8)
        self.q2_blocks = np.zeros(8)
        self.pump_blocks = np.zeros(8)
        self.block_counts = np.zeros(8, int)
        self.r_edges = np.linspace(0., 5., 501)
        self.r_hist = np.zeros(500, np.int64)
        self.inward_hist = np.zeros(500)
        self.last_pair = None
        self.r_sum = None
        self.r2_sum = None
        self.pair_count = 0
        self.max_displacement = 0.
        self.max_pair_shift = 0.

    def initialize(self, occupied):
        shape = occupied.shape
        rng = np.random.default_rng(613)
        sites = np.flatnonzero(occupied)
        self.sites = rng.choice(sites, min(128, len(sites)), replace=False)
        starts = rng.choice(sites, min(512, len(sites)), replace=False)
        coordinates = np.array(np.unravel_index(starts, shape)).T
        ends = []; first = []; directions = []
        for shift in NN:
            if next(int(v) for v in shift if v) < 0:
                continue
            shifted = (coordinates + shift) % shape
            index = np.ravel_multi_index(shifted.T, shape)
            keep = occupied.ravel()[index]
            first.extend(starts[keep]); ends.extend(index[keep])
            directions.extend([self.p['a_A']*(CELL@shift)]*int(keep.sum()))
        self.first = np.array(first)
        self.ends = np.array(ends)
        self.reference = np.array(directions)
        self.rmin = np.full(len(first), np.inf)
        self.speed_min = np.zeros(len(first))
        self.time_min = np.zeros(len(first))
        self.r_sum = np.zeros(len(first)); self.r2_sum = np.zeros(len(first))
        self.occupied = occupied

    def __call__(self, q, occupied, time_ps):
        if not hasattr(self, 'sites'):
            self.initialize(occupied)
        flat = q.reshape(3, -1)
        vectors = self.reference + (flat[:,self.ends]-flat[:,self.first]).T
        distance = np.linalg.norm(vectors, axis=1)
        self.r_sum += distance; self.r2_sum += distance*distance
        self.pair_count += 1
        self.r_hist += np.histogram(distance, self.r_edges)[0]
        self.max_displacement = max(self.max_displacement, float(np.linalg.norm(flat[:,occupied.ravel()], axis=0).max()))
        self.max_pair_shift = max(self.max_pair_shift, float(np.linalg.norm(vectors-self.reference, axis=1).max()))
        smaller = distance < self.rmin
        self.rmin[smaller] = distance[smaller]
        self.time_min[smaller] = time_ps
        if self.last_pair is not None:
            velocity = (vectors-self.last_pair)/(time_ps-self.last_time)*100  # A/ps -> m/s
            self.speed_min[smaller] = np.linalg.norm(velocity[smaller], axis=1)
            inward = np.maximum(0., -(distance-self.last_distance)/((time_ps-self.last_time)*1e-12))
            self.inward_hist += np.histogram(distance, self.r_edges, weights=inward)[0]
        self.last_pair = vectors
        self.last_distance = distance
        self.last_time = time_ps
        self.old.append((q, time_ps))
        if len(self.old) < 3:
            return
        left, middle, right = self.old
        velocity = (right[0]-left[0])/((right[1]-left[1])*1e-12)  # A/s
        middle_q = middle[0][:,occupied]
        t = middle[1]
        first_step = self.case['steps']//2
        time_index = round(t/self.dt_ps)-1-first_step
        block = min(7, max(0, 8*time_index//(self.case['steps']-first_step)))
        v2 = float(np.mean(np.sum(velocity[:,occupied]**2, axis=0)))
        self.v2_blocks[block] += v2
        self.q2_blocks[block] += float(np.mean(np.sum(middle_q**2, axis=0)))
        # Hamiltonian explicit-time power; h is fractional Gamma-curvature modulation.
        pump = -self.case['pump_h']*self.kgamma*self.omega*np.sin(2*self.omega*t)*float(np.mean(middle_q[0]**2))*1e12
        self.pump_blocks[block] += pump
        self.block_counts[block] += 1
        if self.count % 16 == 0:
            self.velocity_samples.append((velocity.reshape(3,-1)[:,self.sites].T*1e-10).astype(np.float32))
            self.velocity_blocks.extend([block]*len(self.sites))
            self.times.append(t)
        self.count += 1
        self.old.pop(0)

    def finish(self, path):
        counts = self.block_counts
        mean = lambda values: np.divide(values, counts, out=np.zeros(8), where=counts>0)
        velocity = np.concatenate(self.velocity_samples)
        np.savez_compressed(path, velocities_m_s=velocity,
            velocity_block=np.array(self.velocity_blocks, np.int8),
            pair_start=self.first, pair_end=self.ends, pair_reference_A=self.reference,
            pair_rmin_A=self.rmin, pair_speed_at_rmin_m_s=self.speed_min,
            pair_tmin_ps=self.time_min, pair_mean_r_A=self.r_sum/self.pair_count,
            pair_mean_r2_A2=self.r2_sum/self.pair_count,
            r_edges_A=self.r_edges, r_counts=self.r_hist,
            inward_flux_A_per_s=self.inward_hist/self.pair_count,
            samples_per_block=counts)
        drag = self.gamma_s*self.mass_eVs2_A2*mean(self.v2_blocks)
        thermal = 3*self.gamma_s*KB*self.case['temperature_K']
        return dict(observation_samples=self.pair_count, dt_ps=self.dt_ps,
            record_ps=self.pair_count*self.dt_ps, retained_velocity_vectors=len(velocity),
            sampled_joint_pairs=len(self.rmin), min_pair_A=float(self.rmin.min()),
            pair_min_quantiles_A=np.quantile(self.rmin,[0,.01,.1,.5,.9,.99,1]).tolist(),
            maximum_observed_site_displacement_A=self.max_displacement,
            maximum_observed_pair_displacement_A=self.max_pair_shift,
            omega_rad_ps=self.omega, gamma_per_s=self.gamma_s,
            ratio4_omega_over_gamma=self.omega*1e12/self.gamma_s,
            q2_A2_blocks=mean(self.q2_blocks).tolist(),
            mean_v2_m2_s2_blocks=(mean(self.v2_blocks)*1e-20).tolist(),
            drag_eV_per_s_per_occupied_site_blocks=drag.tolist(),
            classical_thermal_injection_eV_per_s_per_occupied_site=thermal,
            net_drag_minus_thermal_eV_per_s_per_occupied_site_blocks=(drag-thermal).tolist(),
            explicit_pump_work_eV_per_s_per_occupied_site_blocks=mean(self.pump_blocks).tolist(),
            velocity_estimator='Central difference at every native dt; no change to force/noise. Harmonic v amplitude attenuation sinc(Omega dt); stochastic unresolved kinetic power is not a measured thermostat balance.',
            harmonic_velocity_squared_attenuation=float(np.sinc(self.omega*self.dt_ps/np.pi)**2),
            zero_point_motion='Quantum variance enters Hartree curvature only; classical trajectories are not nuclear quantum samples.')


def main():
    pd_path = ENGINE/'updated/parameters.json'
    ni_path = ENGINE/'updated/ni_epoch/parameters.json'
    pd = json.loads(pd_path.read_text()); ni = json.loads(ni_path.read_text())
    pd_old = json.loads((ENGINE/'updated/pd_summary.json').read_text())['result']['results']
    ni_old = json.loads((ENGINE/'updated/ni_epoch/phase_summary.json').read_text())['result']['summary']
    jobs = []
    for h in (0., .08):
        c = dict(ni_old['configuration'], label=f'NiH_cold_h{h:g}', pump_h=h, stride=1)
        jobs.append((ni, ni_path, c))
    for label in ('random_x1_h0','random_x1_h0.08','thermal_x0.875_h0','thermal_x0.875_h0.08'):
        original = next(r for r in pd_old if r['label']==label)
        c = dict(original['configuration'], label='PdD_'+label, stride=1)
        jobs.append((pd, pd_path, c))
    output = []
    old_cwd = Path.cwd()
    try:
        os.chdir(ENGINE)
        for parameters, source, case in jobs:
            path = HERE/(case['label']+'.npz')
            result_path = HERE/(case['label']+'.json')
            observer = Observer(parameters['materials'][case['material']], case)
            summary, fields = run_case(parameters, case, observe=observer)
            stats = observer.finish(path)
            # Keep scalar, physically useful native summary; do not mirror spectra.
            keep = ('label','configuration','sites','occupied','actual_loading','elapsed_s','state',
                    'static_rms_A','driven_demodulated_rms_A','pump_frequency_rms_A',
                    'incident_E_squared_carrier_THz','peak_fraction','occupation_SHA256',
                    'material_status','loading_response')
            item = dict(summary={k:summary[k] for k in keep}, statistics=stats,
                        parameters_file=str(source.relative_to(ROOT)), parameters_SHA256=digest(source),
                        arrays_file=path.name, arrays_SHA256=digest(path),
                        engine_SHA256=digest(ENGINE/'lattice.py'))
            result_path.write_text(json.dumps(item, indent=2, allow_nan=False)+'\n')
            output.append(result_path.name)
            print(json.dumps(dict(case=case['label'],elapsed_s=summary['elapsed_s'],
                                  rmin_A=stats['min_pair_A'],v2=np.mean(stats['mean_v2_m2_s2_blocks']))), flush=True)
            del observer, fields
    finally:
        os.chdir(old_cwd)
    (HERE/'trajectories.json').write_text(json.dumps(dict(cases=output,
        scope='Same seed, occupancy, physical duration and parameters within each pump contrast; stride changed to1 only for observation.'), indent=2)+'\n')


if __name__ == '__main__':
    main()
