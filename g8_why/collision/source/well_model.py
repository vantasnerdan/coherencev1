"""Force-calibrated cubic anharmonic BO well with physical outer nuclear walls.

The commensurate native forces determine quadratic/quartic/sextic local terms.
All five BO energy changes are unfitted holdouts. Outside the sampled region,
a declared C-infinity continuation joins the same44-electron physical-nucleus
field. No measured fusion yield, density fraction or beam profile is fitted.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from well_field import FrozenHostField, file_binding


def verified(item):
    actual = file_binding(item['path'])
    if (actual['bytes'],actual['sha256']) != (item['bytes'],item['sha256']):
        raise ValueError(f'Changed original input:{item["path"]}')
    return Path(item['path'])


def local_basis(points):
    points = np.asarray(points)
    x2 = points*points
    pairs = x2[...,0]*x2[...,1]+x2[...,0]*x2[...,2]+x2[...,1]*x2[...,2]
    return np.stack([x2.sum(axis=-1)/2,(x2*x2).sum(axis=-1)/4,
                     (x2*x2*x2).sum(axis=-1)/6,pairs/2],axis=-1)


def local_force_basis(points):
    points = np.asarray(points)
    x2 = points*points
    return np.stack([-points,-points**3,-points**5,
                     -points*(x2.sum(axis=-1,keepdims=True)-x2)],axis=-1)


def local_fraction(points, inner, outer):
    radius = np.linalg.norm(points,axis=-1)
    fraction = np.ones_like(radius)
    fraction[radius >= outer] = 0.
    transition = (radius > inner)&(radius < outer)
    t = (radius[transition]-inner)/(outer-inner)
    left,right = np.exp(-1/t),np.exp(-1/(1-t))
    fraction[transition] = right/(left+right)
    return fraction




class CalibratedWell:
    def __init__(self, path):
        self.binding = file_binding(path)
        self.model = json.loads(Path(path).read_bytes())
        self.field = FrozenHostField(verified(self.model['density']))
        self.origin = np.asarray(self.model['origin_A'])
        self.coefficients = np.asarray(self.model['local_coefficients'])
        self.inner,self.outer = self.model['blend_inner_A'],self.model['blend_outer_A']
        self.zero = self.field.points([self.origin]).item()

    def values(self, points):
        points = np.asarray(points)
        result = local_basis(points)@self.coefficients
        fraction = local_fraction(points,self.inner,self.outer)
        outside = fraction < 1
        if np.any(outside):
            wall = self.field.points(points[outside]+self.origin)-self.zero
            result[outside] = fraction[outside]*result[outside]+(1-fraction[outside])*wall
        return result

    def cartesian_grid(self, x, y, z):
        axes = [np.atleast_1d(v) for v in (x,y,z)]
        points = np.stack(np.broadcast_arrays(axes[0][:,None,None],axes[1][None,:,None],
                                             axes[2][None,None,:]),axis=-1)
        fraction = local_fraction(points,self.inner,self.outer)
        local = local_basis(points)@self.coefficients
        wall = self.field.cartesian_grid(*[v+o for v,o in zip(axes,self.origin)])-self.zero
        return fraction*local+(1-fraction)*wall




def calibrate(spec, output):
    started = time.process_time()
    baseline_data = json.loads(verified(spec['baseline']).read_bytes())
    baseline = baseline_data['points'][spec['baseline_row']]
    origin = np.asarray(baseline['positions_A'][baseline_data['target_index']])
    field = FrozenHostField(verified(spec['density']))
    source_rows = []
    for binding in spec['calibration']:
        data = json.loads(verified(binding).read_bytes())
        if data['specification']['hydrogen_setup'] != 'default_PAW':
            raise ValueError('The commissioned BO calibration uses standard H PAW')
        source_rows.extend(data['points'])
    points = np.asarray([row['offset_A'] for row in source_rows])
    actual_forces = np.asarray([row['target_force_eV_A'] for row in source_rows])
    force_basis = local_force_basis(points).reshape(-1,4)
    coefficients,_,rank,singular = np.linalg.lstsq(force_basis,actual_forces.ravel(),rcond=None)
    if rank != 4:
        raise ValueError('Native forces do not identify the declared cubic local well')
    model = dict(schema='g8-force-calibrated-physical-wall-well-v2',origin_A=origin.tolist(),
        density=spec['density'],baseline=spec['baseline'],calibration=spec['calibration'],
        local_coefficients=coefficients.tolist(),
        coefficient_units=['eV/A^2','eV/A^4','eV/A^6','eV/A^4'],
        local_basis=['sum(x_i^2)/2','sum(x_i^4)/4','sum(x_i^6)/6','sum(pair x_i^2*x_j^2)/2'],
        blend_inner_A=spec['blend_inner_A'],blend_outer_A=spec['blend_outer_A'],
        continuation='C-infinity exp(-1/t) partition; pure BO polynomial inside inner radius, physical frozen field outside outer radius',
        specification=spec,singular_values=singular.tolist(),field_diagnostics=field.diagnostics(),
        scope='Clamped one-body BO model; unfitted continuation is a sensitivity, not an exact correlated PdD momentum distribution')
    output.write_text(json.dumps(model,indent=2,sort_keys=True)+'\n')
    well = CalibratedWell(output)
    dF = np.array([row['energy']['free_energy_eV']-baseline['energy']['free_energy_eV']
                   for row in source_rows])
    predicted_forces = np.einsum('ijk,k->ij',local_force_basis(points),coefficients)
    check = dict(schema='g8-force-fit-energy-holdout-consumer-v2',model_binding=file_binding(output),
        points_A=points.tolist(),native_delta_F_eV=dF.tolist(),
        predicted_unfitted_delta_V_eV=well.values(points).tolist(),
        native_forces_eV_A=actual_forces.tolist(),fitted_forces_eV_A=predicted_forces.tolist(),
        energy_holdout_residual_max_eV=float(np.max(abs(well.values(points)-dF))),
        force_fit_residual_max_eV_A=float(np.max(abs(predicted_forces-actual_forces))),
        process_cpu_s=time.process_time()-started,
        limitations=['Commensurate sampling preserves nuclear grid phase, not all finite-grid error',
                     'Outer continuation requires explicit sensitivity; it is not an imposed wavefunction/FFT box cutoff',
                     'Removal energy and exact off-shell fusion amplitudes are not computed'])
    output.with_suffix('.check.json').write_text(json.dumps(check,indent=2,sort_keys=True)+'\n')
    print('G8_CALIBRATED_'+'WELL_COMPLETE',flush=True)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Do not overwrite a published well model')
    calibrate(json.loads(args.input.read_bytes()),args.output)
