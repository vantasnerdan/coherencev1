"""Effective static exponential widths: fixed Pd10/D1, no density normalization."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from project import evaluate

ATOMIC_WIDTHS_A = np.array([.45, .2645885])
SEARCH_A = (.005, 20.)


def atomic_basis(geometry, support):
    grid = support.shape[0]
    lengths = geometry['lengths_A'][:2]
    q = np.fft.fftfreq(grid) * grid
    kx, ky = np.meshgrid(2*np.pi*q/lengths[0], 2*np.pi*q/lengths[1], indexing='ij')
    basis = []
    for label, valence in [('Pd', 10), ('D', 1)]:
        coefficient = np.zeros(support.shape, complex)
        rows = geometry['rows'][label]
        for xy, line in zip(rows['xy_A'], rows['line_number_A_inv']):
            coefficient += valence * line / np.prod(lengths) * np.exp(-1j*(kx*xy[0]+ky*xy[1]))
        coefficient[~support] = 0
        basis.append(coefficient)
    return np.asarray(basis), kx*kx+ky*ky


def model_coefficients(widths_A, geometry, support):
    widths = np.asarray(widths_A)
    if widths.shape != (2,) or not np.isfinite(widths).all() or np.any(widths <= 0):
        raise ValueError('Require two finite positive exponential widths in Angstrom')
    basis, q2 = atomic_basis(geometry, support)
    return np.sum(basis / (1 + widths[:,None,None]**2*q2)**2, axis=0)


def fit_domain(target, basis, q2, mask, starts):
    from scipy.optimize import least_squares
    grid = target.shape[0]
    scale = float(np.sqrt(np.mean(target[mask]**2)))
    cache = {}

    def calculate(logb):
        if 'logb' not in cache or not np.array_equal(logb, cache['logb']):
            b2q = np.exp(2*logb)[:,None,None]*q2
            coefficients = np.sum(basis/(1+b2q)**2, axis=0)
            image = (np.fft.ifft2(coefficients)*grid**2).real
            derivatives = np.fft.ifft2(-4*basis*b2q/(1+b2q)**3,axes=(-2,-1)).real*grid**2
            cache.update(logb=logb.copy(), residual=(image[mask]-target[mask])/scale,
                         jacobian=np.column_stack([derivative[mask]/scale for derivative in derivatives]))
        return cache

    runs, best = [], None
    for start in starts:
        solution = least_squares(lambda p: calculate(p)['residual'], np.log(start),
                                 jac=lambda p: calculate(p)['jacobian'],
                                 bounds=np.log(SEARCH_A), ftol=1e-11,xtol=1e-11,gtol=1e-11,max_nfev=300)
        row = dict(start_A=list(start),widths_A=np.exp(solution.x).tolist(),success=bool(solution.success),
                   nfev=solution.nfev, relative_mean_square=float(np.mean(solution.fun**2)),
                   optimality=float(solution.optimality))
        runs.append(row)
        if best is None or row['relative_mean_square'] < best[0]:
            best = (row['relative_mean_square'], solution)
    solution = best[1]
    singular = np.linalg.svd(solution.jac, compute_uv=False)
    widths = np.exp(solution.x)
    return widths, dict(multistarts=runs, jacobian_log_width_singular_values=singular.tolist(),
        jacobian_condition=float(singular[0]/singular[-1]),
        computational_search_A=list(SEARCH_A),
        near_computational_boundary=bool(np.any(widths < SEARCH_A[0]*1.01) or np.any(widths > SEARCH_A[1]/1.01)),
        note='Jacobian/multistart sensitivity is deterministic model identifiability, not a statistical or physical confidence interval.')


def residuals(coefficients, target, mask, geometry):
    image = (np.fft.ifft2(coefficients)*target.size).real
    difference = image-target
    fractions = np.asarray(geometry['centres_xy_A'])/geometry['lengths_A'][:2]
    return dict(exterior_RMS_e_A3=float(np.sqrt(np.mean(difference[mask]**2))),
        exterior_relative_RMS=float(np.linalg.norm(difference[mask])/np.linalg.norm(target[mask])),
        full_plane_RMS_e_A3=float(np.sqrt(np.mean(difference**2))),
        full_plane_relative_RMS=float(np.linalg.norm(difference)/np.linalg.norm(target)),
        model_centres_n_e_A3=evaluate(coefficients,fractions).tolist(),
        model_mean_n_e_A3=float(image.mean()))


def fit_case(path,geometry):
    with np.load(path,allow_pickle=False) as data:
        target=data['n_perp_e_A3'];support=data['supported_modes'];mask=data['outside_Pd_augmentation']
        target_centres=data['centres_n_e_A3']
    basis,q2=atomic_basis(geometry,support)
    starts=[(pd,d) for pd in (.15,.45,.9) for d in (.12,.2645885,.6)]
    widths,diagnostics=fit_domain(target,basis,q2,mask,starts)
    coefficients=model_coefficients(widths,geometry,support)
    metrics=residuals(coefficients,target,mask,geometry)
    metrics['centre_relative_residual']=(np.abs(np.asarray(metrics['model_centres_n_e_A3'])-target_centres)/target_centres).tolist()
    full_widths,full_diagnostics=fit_domain(target,basis,q2,np.ones_like(mask),[tuple(ATOMIC_WIDTHS_A),tuple(widths),(.2,.2)])
    full_metrics=residuals(model_coefficients(full_widths,geometry,support),target,mask,geometry)
    return dict(widths_A=widths.tolist(),widths_nm=(widths/10).tolist(),
        ratios_to_atomic_length_table=(widths/ATOMIC_WIDTHS_A).tolist(),
        residuals=metrics,diagnostics=diagnostics,
        full_plane_fit_diagnostic=dict(widths_A=full_widths.tolist(),residuals=full_metrics,diagnostics=full_diagnostics))


def kernel_qualification():
    from scipy.integrate import quad
    errors=[]
    for b in (.2,.45,.8):
        for q in (0.,1.,3.,6.):
            radial=quad(lambda r:r*r*np.exp(-r/b)*np.sinc(q*r/np.pi)/(2*b**3),0,np.inf,
                        epsabs=1e-11,epsrel=1e-11)[0]
            errors.append(abs(radial-(1+b*b*q*q)**-2))
    if max(errors)>1e-9:
        raise ValueError('Real-space exponential radial Fourier oracle failed')
    return dict(independent_radial_quadrature_max_error=max(errors),tolerance=1e-9)


def run(readout,report):
    cpu,wall=time.process_time(),time.perf_counter()
    geometry=json.loads((readout/'geometry.json').read_bytes())
    qualification=kernel_qualification()
    cases={name:fit_case(readout/(name+'-transverse.npz'),geometry) for name in ('h020','h016')}
    coarse,fine=[np.asarray(cases[name]['widths_A']) for name in ('h020','h016')]
    result=dict(schema='pdd-exponential-width-fit-v1',species_order=['Pd','D'],cases=cases,
        actual_cell_a_A=geometry['a_A'],upstream_model_cell_a_A=4.084,
        h_relative_width_change=(np.abs(coarse-fine)/fine).tolist(),
        atomic_length_table_A=ATOMIC_WIDTHS_A.tolist(),kernel_qualification=qualification,
        contract='Static exp(-r/b)/(8*pi*b^3), fixed10/1valence weights; no Debye-Waller, amplitude, charge transfer or count correction. Primary fit outside Pd matching cylinders; centre predictions are extrapolations when inside those cylinders.',
        process_cpu_s=time.process_time()-cpu,wall_s=time.perf_counter()-wall)
    with report.open('x') as handle:
        handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--readout',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    run(args.readout,args.report)
