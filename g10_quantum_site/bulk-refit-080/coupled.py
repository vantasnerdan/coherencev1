"""Same archived paths/callback, one changed site; paired .48/.80 comparisons."""
import csv
import hashlib
import json
import sys
import numpy as np
import well

sys.path.insert(0, str(well.HERE.parent))
c = well.load_module('pdd048_reweight', well.PRIOR / 'coupled_sensitivity.py')
c.EXTERIORS = (.24, .48, .80)
c.DISPLACEMENTS = np.array([0., .25, .50])
VARIANTS = ('old', 'quadratic', 'quartic', *well.VARIANTS)

def main():
    out = well.HERE / 'results'
    out.mkdir(exist_ok=True)
    old = c.full('Pd', 'joint')
    site = old['site']
    records, common, curves = [], [], []
    for temperature, beads in ((300,32), (77,128)):
        path = well.HERE.parent / f'output/coupled/PdD_{temperature}K_paths.npz'
        with np.load(path) as data:
            x, field = data['beads_A'], data['neighbor_field_eV_A']
        if x.shape != (801, beads, 3) or field.shape != x.shape:
            raise ValueError('Archived path dimensions differ from accepted source')
        packed = np.stack((x,field), axis=2)
        beta = 1/(c.KB_EV_PER_K*temperature)
        columns, histogram = c.path_observables(x)
        names = list(columns)
        obs = np.column_stack([*columns.values(), histogram])
        old_bead = well.energy_force(x,'old')[0]
        old_local = dict(old,K=old['K'][site:site+1], f=old['f'][site:site+1], J=np.zeros((1,1)))
        control = np.max(abs(c.potential(x[...,None,:],old_local)[0]-old_bead))
        if control > 1e-12:
            raise ValueError('Old callback no longer matches the fixed graph Hamiltonian')
        estimates, chars, variants = {}, {}, {}
        for variant in VARIANTS:
            local = lambda q, v=variant: well.energy_force(q,v)[0]
            logw = -beta*(local(x)-old_bead).mean(axis=1)
            weights = np.exp(logw-logw.max())
            val, blocks = c.ratio_stats(obs,weights)
            estimates[variant] = val,blocks
            participation = dict(frame_ess=c.ess(weights),block_ess=c.ess(weights.reshape(9,89).sum(axis=1)),
                                 log_weight_span=float(np.ptp(logw)))
            if participation['frame_ess'] < 50 or participation['block_ess'] < 4:
                raise ValueError('Insufficient sampled weight overlap; do not continue to characteristic')
            metrics = {}
            for j,name in enumerate(names):
                positive = int(np.count_nonzero((weights*obs[:,j]).reshape(9,89).sum(axis=1)))
                metrics[name] = dict(value=float(val[j]), block_se=float(c.se(blocks)[j]),
                    nonzero_numerator_blocks=positive,
                    change_vs_old=float(val[j]-estimates['old'][0][j]),
                    paired_change_block_se=float(c.se(blocks-estimates['old'][1])[j]))
            char = c.characteristic(packed,temperature,local,logw)
            chars[variant] = char
            direct, share = c.direct_witness(x,field,old_bead,local,logw,beta)
            direct_error = abs(direct-char.value[1])
            if abs(char.value[0]-1)>1e-14 or direct_error>1e-12:
                raise ValueError('Finite-P shifted-target witness disagrees')
            for i,a in enumerate(c.DISPLACEMENTS):
                curves.append(dict(temperature_K=temperature,variant=variant,separation_A=float(a),
                    characteristic=float(char.value[i]),block_se=float(char.se[i]),
                    change_vs_old=float(char.value[i]-chars['old'].value[i]),
                    paired_change_block_se=float(c.se(char.blocks-chars['old'].blocks)[i]),
                    numerator_block_ess=float(char.numerator_block_ess[i])))
            variants[variant] = dict(position=metrics,overlap=participation,
                shifted_target_control_abs_error=float(direct_error),
                witness_numerator_share_any_shifted_bead_gt_048=float(share))
            tail80 = metrics['bead_fraction_gt_0.80_A']
            common.append(dict(material='PdD',isotope='D',temperature_K=temperature,variant=variant,
                msd_A2=metrics['MSD_A2']['value'],kinetic_eV=None,
                probability_r_gt_048=metrics['bead_fraction_gt_0.48_A']['value'],
                probability_r_gt_080=tail80['value'] if tail80['nonzero_numerator_blocks']>=4 else None,
                support=dict(consumer='one-site finite-P reweight',beads=beads,frames=len(x),
                             probability_r_gt_080_raw=tail80['value'],
                             probability_r_gt_080_nonzero_blocks=tail80['nonzero_numerator_blocks'],
                             kinetic='not inferred from three characteristic samples')))
        pairs = {}
        for first,second in [('quadratic','wide_quadratic'),('quartic','wide_quartic'),
                             ('quadratic','quartic'),('wide_quadratic','wide_quartic')]:
            v,b = estimates[first]; vv,bb = estimates[second]
            pairs[f'{second}_minus_{first}'] = {
                name:dict(change=float(vv[j]-v[j]),paired_block_se=float(c.se(bb-b)[j]))
                for j,name in enumerate(names)}
        records.append(dict(temperature_K=temperature,beads=beads,frames=len(x),
            original_atom_index=int(old['original_indices'][site]),
            path=str(path.relative_to(well.HERE.parents[1])),path_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            old_callback_control_max_abs_eV=float(control),variants=variants,paired_comparisons=pairs))
    result = dict(scope='one changed on-site term only; same original D mass, atom order, K, F0, host graph and archived joint paths',
        weight='exp[-beta mean_beads(V_variant-V_old)]',
        uncertainty='nine contiguous 89-frame ratio influence blocks; paired covariance retained; no P extrapolation or unseen-overlap guarantee',
        momentum='exact existing finite-P signed [111] displaced-path characteristic, not momentum density or tail',
        rare_tail='zero sampled >.80 counts are unobserved, not a probability or zero-error claim; common observable is null unless at least four blocks contribute',
        cases=records,observable_rows=common,characteristic_rows=curves)
    (out/'coupled.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(records,indent=2,allow_nan=False))

if __name__ == '__main__':
    main()
