"""One algebraic DC+2Omega / DC+2Omega+Omega probe on run_case's actual sampler.
No trajectory, changed material, parameter sweep, relaxation or SCF is run.
"""
import json
from hashlib import sha256
from pathlib import Path
import numpy as np
from lattice import OBSERVATION_DEFAULTS,harmonic_projection


def main():
    root=Path(__file__).resolve().parent
    quantum=json.loads((root/'pd_x1_quantum'/'result.json').read_text())
    row=next(value for value in quantum['rows'] if value['material']=='PdD')
    pumped=row['pumped'];units=row['units_bath']
    # run_case's clock uses classical Gamma Hartree curvature, not the distinct RPA pole.
    HF_energy=pumped['mode_quantum_meV']/pumped['RPA_clock_over_classical_Hartree_clock']
    omega=HF_energy/units['bare_harmonic_energy_unit_meV']
    profile=dict(OBSERVATION_DEFAULTS,pump_ratio=1.)
    steps=profile['steps'];dt=profile['dt'];stride=profile['stride']
    sampled_steps=np.arange(steps//2,steps)
    sampled_steps=sampled_steps[sampled_steps%stride==0]
    times=(sampled_steps+1)*dt
    phase=np.exp(-1j*omega*profile['pump_ratio']*times);phase2=phase*phase
    basis=np.stack((phase.real,-phase.imag,phase2.real,-phase2.imag,np.ones(len(phase))))
    gram=basis@basis.T
    config=json.loads((root/'new_inputs'/'loading_b_793'/'config.json').read_text())
    direction=np.asarray(next(value for value in config['probes'] if value['probe_id']=='near')['direction_vector'])
    DC=.03*direction
    pump=.006*(1-.35j)*direction
    half=.008*(1+.4j)*direction
    forced=DC[:,None]+np.real(pump[:,None]*phase2.conj())
    mixed=forced+np.real(half[:,None]*phase.conj())
    def project(record):
        return harmonic_projection(record.sum(axis=1),(record*phase).sum(axis=1),(record*phase2).sum(axis=1),gram)
    dc0,half0,pump0=project(forced)
    dc1,half1,pump1=project(mixed)
    old_basis=basis[[0,1,4]]
    old_coefficients=np.linalg.solve(old_basis@old_basis.T,old_basis@forced.T)
    old_false_half=old_coefficients[0]-1j*old_coefficients[1]
    result={'scope':'Algebraic measurement probe ONLY; amplitudes below are injected controls, not material predictions',
        'clock':'run_case classical Gamma Hartree reference, pump_ratio1; NOT the distinct quantum RPA pole',
        'sampling_profile':profile,'sample_count':int(len(times)),
        'first_last_sample_steps':[int(sampled_steps[0]),int(sampled_steps[-1])],
        'time_interval_ps':[float(times[0]*units['time_unit_ps']),float(times[-1]*units['time_unit_ps'])],
        'angular_half_pump_frequency_code':float(omega),'Gram_rank':int(np.linalg.matrix_rank(gram)),
        'Gram_condition_number':float(np.linalg.cond(gram)),
        'pure_DC_plus_pump':{'legacy_DC_only_false_half_amplitude_A':float(np.linalg.norm(old_false_half)),
            'joint_false_half_amplitude_A':float(np.linalg.norm(half0)),
            'DC_max_error_A':float(np.max(abs(dc0-DC))),'pump_complex_max_error_A':float(np.max(abs(pump0-pump)))},
        'plus_true_half_pump':{'injected_half_vector_norm_A':float(np.linalg.norm(half)),
            'half_complex_max_error_A':float(np.max(abs(half1-half))),
            'DC_max_error_A':float(np.max(abs(dc1-DC))),'pump_complex_max_error_A':float(np.max(abs(pump1-pump)))},
        'injected_control':{'DC_vector_A':DC.tolist(),'pump_real_A':pump.real.tolist(),'pump_imag_A':pump.imag.tolist(),
            'half_real_A':half.real.tolist(),'half_imag_A':half.imag.tolist()},
        'mathematical_scope':'Full-rank joint least squares annihilates the DC/2Omega nuisance columns exactly in exact arithmetic; this is not an alias-free certificate for arbitrary unmodelled frequencies',
        'lattice_source_SHA256':sha256((root/'lattice.py').read_bytes()).hexdigest()}
    errors=[result['pure_DC_plus_pump']['joint_false_half_amplitude_A'],result['pure_DC_plus_pump']['DC_max_error_A'],
        result['pure_DC_plus_pump']['pump_complex_max_error_A'],result['plus_true_half_pump']['half_complex_max_error_A'],
        result['plus_true_half_pump']['DC_max_error_A'],result['plus_true_half_pump']['pump_complex_max_error_A']]
    if max(errors)>1e-13:raise AssertionError('Actual joint sampler does not separate the prescribed tones')
    destination=root/'new_inputs'/'harmonic_separation';destination.mkdir(parents=True,exist_ok=True)
    (destination/'signal_proof.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
