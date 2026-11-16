"""Finite-q RPA from saved PW wavefunctions; no Thomas-Fermi substitution."""
import time

import numpy as np
from ase.units import Bohr
from gpaw import GPAW
from gpaw.mpi import serial_comm
from gpaw.response.df import DielectricFunction

from screen_io import manifest, require, write_json, write_npz


def response_fingerprint(config, output):
    gpw = output/'ground.gpw'
    calc = GPAW(str(gpw), legacy_gpaw=False,
                communicator=serial_comm, txt=None)
    cell = np.asarray(config['geometry']['cell_A'])
    reciprocal = 2*np.pi*np.linalg.inv(cell).T
    rows = []
    for index, settings in enumerate(config['responses']):
        start = time.monotonic()
        require(settings['nbands'] <= config['ground_state']['saved_nbands'],
                'Response requests unsaved empty states')
        q = np.asarray(settings['q_fractional'], float)
        require(np.linalg.norm(q) > 0, 'Optical q=0 is not this fingerprint')
        require(np.allclose(q*np.asarray(config['ground_state']['kmesh']),
                            np.rint(q*np.asarray(config['ground_state']['kmesh']))),
                'q does not map the saved uniform k grid onto itself')
        tag = f'response-{index:02d}'
        df = DielectricFunction(
            calc, frequencies=[0.0], ecut=settings['ecut_eV'],
            hilbert=False, eta=settings['eta_eV'], nbands=settings['nbands'],
            intraband=True, rate='eta', nblocks=1, world=serial_comm,
            integrationmode='point integration',
            txt=str(output/(tag+'.txt')))
        result = df.get_inverse_dielectric_function(q_c=q.tolist(), xc='RPA')
        vectors = result.qpd.get_reciprocal_vectors(add_q=True)/Bohr
        require(np.allclose(vectors[0], q@reciprocal, rtol=1e-10, atol=1e-10),
                'The presumed dielectric head is not the supplied physical q')
        chi0 = result.vchi0_symm_wGG[0]
        chi = result.vchi_symm_wGG[0]
        identity = np.eye(len(vectors), dtype=complex)
        epsilon = identity-chi0
        inverse = identity+chi
        require(epsilon.shape == inverse.shape == (len(vectors), len(vectors)),
                'Unexpected response matrix shape')
        require(np.isfinite(epsilon).all() and np.isfinite(inverse).all(),
                'Nonfinite response matrix')
        residual = float(np.linalg.norm(epsilon@inverse-identity, ord=np.inf))
        require(residual < 1e-8, 'Dielectric/Dyson inverse identity failed')
        inverse_nlf = 1/np.diag(epsilon)
        inverse_lf = np.diag(inverse)
        arrays = write_npz(
            output/(tag+'.npz'), Q_Ainv=vectors,
            epsilon_coulomb_symmetrized=epsilon,
            inverse_epsilon_coulomb_symmetrized=inverse,
            inverse_diagonal_no_local_fields=inverse_nlf,
            inverse_diagonal_local_fields=inverse_lf)
        row = dict(settings=settings, Q_head_Ainv=float(np.linalg.norm(vectors[0])),
                   matrix_dimension=len(vectors),
                   inverse_identity_infinity_residual=residual,
                   epsilon_head_no_local_fields=[float(epsilon[0, 0].real),
                                                  float(epsilon[0, 0].imag)],
                   epsilon_macro_local_fields=[float((1/inverse[0, 0]).real),
                                                float((1/inverse[0, 0]).imag)],
                   arrays=arrays, elapsed_wall_s=time.monotonic()-start)
        write_json(output/(tag+'.json'), row)
        rows.append(row)
        print(f'RPA_FINGERPRINT_POINT_COMPLETED index={index}', flush=True)
    record = dict(schema='coherence-screen-response-fingerprint-v1',
                  case=config['case'], model=config['model'],
                  gpw=manifest(gpw), points=rows,
                  eligible=False,
                  reason='API/cost fingerprint; k, eta, angular and q-tail errors unbounded')
    write_json(output/'response.json', record)
    return record
