"""Small, frozen-profile Pd/PdO beam calculation; run from any directory.

One-thread NumPy/SciPy only. No source data fitting or target evolution.
Outputs are per true incident D+ charge, normalized on the SAME foil at 10 keV.
The energy-integral implementation is checked below against csda_nodes.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
from functools import lru_cache
from pathlib import Path
import csv
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar, brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from g4_beam.forward import consumer as C, model as M
from g4_beam.transport import bca

HERE = Path(__file__).resolve().parent
E_LAB = np.arange(2.5, 10.01, .5)
E = np.linspace(.1, 13., 25801)
N_PD = 6.80e28                    # fixed Pd atom density; no fitted density
X_BULK = .01                     # declared sensitivity input, NOT measured local D/Pd
N_BULK = X_BULK*N_PD
# Crystalline PdO cell: two PdO units, a=.3043nm, c=.5336nm (source packet [B]).
N_PDO = 2/((.3043e-9)**2*.5336e-9)
QE = C.QE
SIGMA = M.reaction_sigma(E/2, .040, 'p', 'eq7')

# Exactly the consumer's Pd/D stopping evaluated once, then Bragg-linearly combined.
_LAYER0 = (0., 1e-6, N_PD, 0., 46., 106.42, N_PD)
_LAYER1 = (0., 1e-6, N_PD, N_PD, 46., 106.42, N_PD)
S_PD = np.array([C.stopping_eV_per_m(e, _LAYER0) for e in E])
S_D = np.array([C.stopping_eV_per_m(e, _LAYER1) for e in E])-S_PD
# Actual Pd+O composition, not a generic/effective atomic number. Pd and D
# retain ICRU/ZBL; O uses the repository's estimated LS/ZBL convention.
S_O_ATOM = (bca.se_ls(E*1000, 8)+bca.sn_zbl(E*1000, 8, 15.999))*1e-18
S_OX = S_PD*N_PDO/N_PD + N_PDO*S_O_ATOM


def cum(values):
    return cumulative_trapezoid(values, E, initial=0.)


class Column:
    """Constant-composition CSDA layer: dx=1000 dE/S, E in lab keV."""
    def __init__(self, stopping, nd):
        self.s = stopping
        self.nd = nd
        self.path = cum(1000/stopping)
        self.kernel = cum(SIGMA*1000/stopping)*nd*1e-31/QE

    def segment(self, incident, thickness_nm):
        start = np.interp(incident, E, self.path)
        remaining = np.maximum(0., start-thickness_nm*1e-9)
        exit_energy = np.interp(remaining, self.path, E)
        y = np.interp(incident, E, self.kernel)-np.interp(exit_energy, E, self.kernel)
        return exit_energy, y

    def thick(self, incident):
        return np.interp(incident, E, self.kernel)


@lru_cache(maxsize=None)
def pd_column(x):
    return Column(S_PD+x*S_D, x*N_PD)


@lru_cache(maxsize=None)
def ox_column(x, oxygen_scale=1.):
    return Column(S_PD*N_PDO/N_PD+oxygen_scale*N_PDO*S_O_ATOM+x*S_D*N_PDO/N_PD,
                  x*N_PDO)


BULK = pd_column(X_BULK)


def normalized(y):
    return y/y[-1]


def uniform(energies, u_keV=.040):
    kernel = cum(M.reaction_sigma(E/2, u_keV, 'p', 'eq7')*1000/BULK.s)*N_BULK*1e-31/QE
    return np.interp(energies, E, kernel)


def fit_u(y):
    """Equal log-yield weights at 2.5,3,...,9.5 lab; 10 lab is exact reference."""
    target = np.log(normalized(y)[:-1])
    def objective(u):
        return np.mean((np.log(normalized(uniform(E_LAB, u))[:-1])-target)**2)
    opt = minimize_scalar(objective, bounds=(0., 2.), method='bounded', options={'xatol': 1e-10})
    fit = normalized(uniform(E_LAB, opt.x))
    residual = np.log(normalized(y)/fit)
    return dict(U_app_eV=float(opt.x*1000), log_rms=float(np.sqrt(np.mean(residual[:-1]**2))),
                max_fractional_shape_residual=float(np.max(np.abs(np.expm1(residual[:-1])))))


def pd_profile(t_nm, ratio, feedback=True):
    # Frozen stopping uses bulk composition while changing only reaction density.
    layer = pd_column(X_BULK*ratio) if feedback else BULK
    exit_e, yl = layer.segment(E_LAB, t_nm)
    if not feedback:
        yl = yl*ratio
    return yl+BULK.thick(exit_e)


def oxide_profile(t_nm, ratio, buried=False, oxygen_scale=1.):
    # Exposed D-rich region is IN the cap: 0<t<30 nm. Buried region begins at
    # the 30 nm oxide/metal interface; in that case the entire cap is D-free.
    if buried:
        after_cap, _ = ox_column(0., oxygen_scale).segment(E_LAB, 30.)
        after_layer, yl = pd_column(X_BULK*ratio).segment(after_cap, t_nm)
        return yl+BULK.thick(after_layer)
    after_rich, yr = ox_column(N_BULK*ratio/N_PDO, oxygen_scale).segment(E_LAB, t_nm)
    after_cap, _ = ox_column(0., oxygen_scale).segment(after_rich, 30.-t_nm)
    return yr+BULK.thick(after_cap)


def stop_oxide(energy, layer):
    _, _, nh, nd, _, _, _ = layer
    # Same convention as S_OX but scalar for authoritative csda_nodes parity.
    pd = C.stopping_eV_per_m(energy, (0., 1e-6, nh, nd, 46., 106.42, nh))
    oxygen = nh*float(bca.se_ls(energy*1000, 8)+bca.sn_zbl(energy*1000, 8, 15.999))*1e-18
    return pd+oxygen


def consumer_yield(energies, layers, stopping_fn=C.stopping_eV_per_m, u=.040):
    result = []
    for e in energies:
        en, _, w, _ = C.csda_nodes(float(e), layers, stopping_fn)
        result.append(float(np.sum(w*M.reaction_sigma(en/2, u, 'p', 'eq7')))/QE)
    return np.array(result)


def write_csv(filename, rows):
    with open(HERE/filename, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    ts = [.05, .1, .3, 1., 3., 10., 30., 100.]
    rs = [1., 2., 5., 10., 30., 60., 100., 1e3, 1e5]
    rows = []
    curves = []
    def add_curve(label, y):
        fit = fit_u(y)
        fitted = normalized(uniform(E_LAB, fit['U_app_eV']/1000))
        for e, yy, rr, ff in zip(E_LAB, y, normalized(y), fitted):
            curves.append(dict(scenario=label, E_lab_keV=float(e), E_cm_keV=float(e/2),
                               reactions_per_true_C=float(yy), R_own_10keV=float(rr),
                               fitted_uniform_R=float(ff), fractional_shape_residual=float(rr/ff-1)))
        return fit
    for mode in ('Pd_frozen', 'Pd_D_feedback', 'in_30nm_PdO', 'under_30nm_PdO'):
        for t in ts:
            if mode == 'in_30nm_PdO' and t>30:
                continue
            for r in rs:
                if mode.startswith('Pd_'):
                    y = pd_profile(t, r, mode == 'Pd_D_feedback')
                else:
                    y = oxide_profile(t, r, mode == 'under_30nm_PdO')
                d_per_host = N_BULK*r/N_PDO if mode == 'in_30nm_PdO' else X_BULK*r
                rows.append(dict(geometry=mode, t_nm=t, r=r, D_per_host_in_rich_layer=d_per_host,
                                 capacity_label='hypothetical_inventory_not_measurement' if d_per_host<=1 else 'exceeds_single_site_inventory_benchmark',
                                 R_2p5_over_10=float(y[0]/y[-1]), **fit_u(y)))
    write_csv('map.csv', rows)
    baseline = uniform(E_LAB)
    sheet = M.reaction_sigma(E_LAB/2, .040, 'p', 'eq7')*1e19*1e-31/QE
    sheet_fit = add_curve('exposed_thin_sheet_1e19_m-2_frozen_stop_limit', sheet)
    add_curve('uniform_true40eV', baseline)
    add_curve('uniform_reported310eV_reference', uniform(E_LAB, .310))
    add_curve('uniform_reported600eV_reference', uniform(E_LAB, .600))
    for mode, t, r in [('Pd_frozen', 1., 100.), ('Pd_D_feedback', 1., 100.),
                       ('in_30nm_PdO', 1., 100.), ('under_30nm_PdO', 1., 100.)]:
        y = pd_profile(t, r, mode == 'Pd_D_feedback') if mode.startswith('Pd_') else oxide_profile(t, r, mode.startswith('under'))
        add_curve(mode+'_t1nm_r100', y)

    # Analytic positive-mixture identity: B+(r-1)L. At each energy normalized
    # ratios are a convex combination of normalized B and L with 10-keV weights.
    envelope = dict(sheet_fit=sheet_fit, sheet_over_uniform40_at_2p5=float(normalized(sheet)[0]/normalized(baseline)[0]))
    bare = uniform(E_LAB, 0.)
    envelope['sheet_over_same_model_bare_at_2p5'] = float(normalized(sheet)[0]/normalized(bare)[0])
    envelope['uniform40_over_same_model_bare_at_2p5'] = float(normalized(baseline)[0]/normalized(bare)[0])
    envelope['uniform600_over_same_model_bare_at_2p5'] = float(normalized(uniform(E_LAB,.600))[0]/normalized(bare)[0])
    envelope['source_comparison'] = 'Kasagi reports about50x its own adjusted bare curve at2.5lab; not identical to this commissioned Eq7/Bosch-Hale+ICRU/ZBL comparison, and not a raw-data refit'
    for u in (.310, .600):
        requested = normalized(uniform(E_LAB, u))
        envelope[str(int(u*1000))+'eV'] = dict(
            requested_over_uniform40_at_2p5=float(requested[0]/normalized(baseline)[0]),
            requested_over_sheet_at_2p5=float(requested[0]/normalized(sheet)[0]),
            frozen_Pd_reachable=False, t_max_nm=None, r_required=None,
            reason='Target lowest-point normalized yield exceeds exposed-sheet envelope')
    # Local energy derivative determines whether fixed-depth ratio h(z) decreases.
    ec = E/2
    a = M.BH['p']
    sp = a[1]+ec*(2*a[2]+ec*(3*a[3]+ec*4*a[4]))
    dlog_lab = .5*(sp/M.S_BH(ec, 'p')-.5/ec-.5/(ec+.040)+.5*M.BG/(ec+.040)**1.5)
    envelope['species_monotonic_flow_max_delta'] = {
        name: float(np.max(np.diff(s[E<=10]*dlog_lab[E<=10])))
        for name, s in [('Pd', S_PD), ('D', S_D), ('O', S_O_ATOM*N_PD)]}
    envelope['flow_check_meaning'] = 'finite-grid sign checks of each Bragg component s_j*d(log sigma)/dE over 0.1-10 lab; sufficient condition is analytic, sampled signs are not an interval proof'
    envelope['map_maxima'] = {mode: max((r for r in rows if r['geometry']==mode), key=lambda r:r['U_app_eV']) for mode in sorted({r['geometry'] for r in rows})}
    envelope['hypothetical_D_per_host_le1_map_maxima'] = {mode: max((r for r in rows if r['geometry']==mode and r['D_per_host_in_rich_layer']<=1), key=lambda r:r['U_app_eV']) for mode in sorted({r['geometry'] for r in rows})}
    envelope['oxide_stopping_sensitivity_t1nm_r30'] = [
        dict(oxygen_scale=scale, **fit_u(oxide_profile(1., 30., False, scale)))
        for scale in (.5, 1., 2.)]

    # Constant actual LAB offset applied equally at all scan energies AND monitor.
    # D+ energy gain delta_keV corresponds to potential phi=-1000*delta volts.
    offsets = []
    for base, target in [(.040, .340), (.310, .610), (.310, .600), (.040, .600)]:
        def shifted(delta):
            return uniform(E_LAB+delta, base)
        delta = brentq(lambda d:fit_u(shifted(d))['U_app_eV']-target*1000, 0., 2.)
        fit = add_curve(f'constant_lab_gain_{base*1000:g}_to_{target*1000:g}eV', shifted(delta))
        offsets.append(dict(U_true_eV=base*1000, target_U_app_eV=target*1000,
                            actual_lab_gain_eV=float(delta*1000), potential_V=float(-delta*1000),
                            field_V_per_m=float(delta*1000/30e-9), **fit))
    # Source held P_lab nominally fixed: I/I10=10/E. A purely resistive charge
    # would follow this current rather than a constant potential. This is a
    # kinematic WRONG-SIGN counterfactual, not a predicted negative charge.
    current_delta = brentq(lambda d:fit_u(uniform(E_LAB+d*10/E_LAB, .040))['U_app_eV']-340., 0., 2.)
    fit = add_curve('negative_potential_proportional_to_current_wrong_sign_control', uniform(E_LAB+current_delta*10/E_LAB, .040))
    current_control = dict(monitor_potential_V=float(-1000*current_delta), low_energy_potential_V=float(-4000*current_delta), **fit)
    positive_control = uniform(E_LAB-current_delta*10/E_LAB, .040)
    current_control['positive_ion_charge_same_magnitude_R2p5_relative_to_uniform40'] = float(normalized(positive_control)[0]/normalized(baseline)[0])

    # Even an unlimited exposed sheet plus an attractive potential must satisfy
    # the lowest point before a fit could succeed. This is an optimistic necessary
    # boundary, NOT a simultaneous fit or an actual oxide electrical mechanism.
    combined_bounds = []
    for target in (.310, .340, .600):
        wanted = normalized(uniform(E_LAB, target))[0]
        def sheet_low(delta):
            ss = M.reaction_sigma(np.array([2.5+delta, 10.+delta])/2, .040, 'p', 'eq7')
            return ss[0]/ss[1]
        delta = brentq(lambda d:sheet_low(d)-wanted, 0., 2.)
        combined_bounds.append(dict(target_U_app_eV=1000*target,
                                    optimistic_minimum_abs_negative_potential_V=float(1000*delta),
                                    field_V_per_m=float(1000*delta/30e-9)))
    # Positive voltage can fake enhancement only if it suppresses the monitor
    # much more strongly than the low-energy point. Price that loophole at the
    # lowest point; normal constant-power Ohmic charging has the opposite order.
    positive_monitor_bounds = []
    for base, target in [(.040, .340), (.310, .600)]:
        wanted = normalized(uniform(E_LAB, target))[0]
        low = uniform(np.array([2.5]), base)[0]
        delta = brentq(lambda d:low/uniform(np.array([10.-d]), base)[0]-wanted, 0., 5.)
        positive_monitor_bounds.append(dict(U_true_eV=base*1000, target_U_app_eV=target*1000,
                                           positive_monitor_only_potential_V=float(1000*delta)))

    # Actual paths, exits and yield-weighted depths are DIFFERENT quantities.
    depth_rows = []
    depth_grid = []
    for e in E_LAB:
        p0 = np.interp(e, E, BULK.path)
        en = E[E<=e]
        z = p0-np.interp(en, E, BULK.path)
        integrand = M.reaction_sigma(en/2, .040, 'p', 'eq7')*1000/np.interp(en, E, BULK.s)
        meanz = np.trapezoid(z*integrand, en)/np.trapezoid(integrand, en)
        after5 = np.interp(max(0.,p0-5e-9), BULK.path, E)
        fraction5 = 1-float(BULK.thick(after5)/BULK.thick(e))
        exits = []
        for scale in (.5, 1., 2.):
            col = ox_column(0., scale)
            cap_exit, _ = col.segment(e, 30.)
            fullpath = float(np.interp(e, E, col.path)*1e9)
            exits.append(float(cap_exit))
            depth_grid.append(dict(E_lab_keV=float(e), oxygen_stopping_scale=scale,
                                   oxide_CSDA_path_to_0p1keV_nm=fullpath,
                                   exit_30nm_cap_keV=float(cap_exit), buried_interface_accessible=bool(cap_exit>.100001)))
        depth_rows.append(dict(E_lab_keV=float(e), Pd_CSDA_path_to_0p1keV_nm=float(p0*1e9),
                               Pd_mean_reaction_depth_nm=float(meanz*1e9), Pd_yield_fraction_first5nm=fraction5,
                               nominal_30nm_PdO_exit_lab_keV=exits[1]))
    write_csv('depths.csv', depth_rows)
    write_csv('oxide_paths.csv', depth_grid)

    # Only narrow integrated consumer checks: homogeneous, exposed step,
    # oxide+buried step, same Eq7 baseline. No unrelated project suite.
    check_e = np.array([2.5, 5., 10.])
    cases = []
    uniform_layers = [(0., 1000e-9, N_PD, N_BULK, 46., 106.42, N_PD)]
    cases.append(('uniform', uniform(check_e), consumer_yield(check_e, uniform_layers)))
    layers = [(0., 1e-9, N_PD, 100*N_BULK, 46., 106.42, N_PD),
              (1e-9, 1000e-9, N_PD, N_BULK, 46., 106.42, N_PD)]
    cases.append(('surface_step_with_D_stopping', pd_profile(1., 100.)[[0,5,15]], consumer_yield(check_e, layers)))
    # A composition callback dispatches on the explicitly marked cap Z=8.
    oxide_layers = [(0., 30e-9, N_PDO, 0., 8., 15.999, N_PDO),
                    (30e-9, 31e-9, N_PD, 100*N_BULK, 46., 106.42, N_PD),
                    (31e-9, 1000e-9, N_PD, N_BULK, 46., 106.42, N_PD)]
    def mixed_stop(e, layer):
        return stop_oxide(e, layer) if layer[4]==8 else C.stopping_eV_per_m(e, layer)
    cases.append(('oxide_buried_step', oxide_profile(1.,100.,True)[[0,5,15]], consumer_yield(check_e, oxide_layers, mixed_stop)))
    in_cap_layers = [(0., 1e-9, N_PDO, 30*N_BULK, 8., 15.999, N_PDO),
                     (1e-9, 30e-9, N_PDO, 0., 8., 15.999, N_PDO),
                     (30e-9, 1000e-9, N_PD, N_BULK, 46., 106.42, N_PD)]
    cases.append(('oxide_exposed_step_exact_density_ratio', oxide_profile(1.,30.)[[0,5,15]],
                  consumer_yield(check_e, in_cap_layers, mixed_stop)))
    parity = {name: float(np.max(np.abs(actual/reference-1))) for name, actual, reference in cases}
    if max(parity.values()) > .002:
        raise AssertionError(('energy integral / csda_nodes parity', parity))
    # Exercise forward_yields normalization/current semantics explicitly; source
    # Kasagi beam is D+, not the unrelated default 85/15 molecular beam.
    actual = C.forward_yields(check_e/2, lambda e:(uniform_layers,0.,{}), Ue_keV=.040,
                              reaction='eq7', reference_E_cm_keV=5.,
                              measurement=dict(D1_fraction=1.,D2_fraction=0.,spread_FWHM=0.))
    parity['forward_yields_R'] = float(np.max(np.abs(actual['R']/normalized(uniform(check_e))-1)))
    if parity['forward_yields_R'] > .002:
        raise AssertionError(parity)

    # Physical current/field/leakage demands; primary gives 60uA conditioning at
    # 10keV, not a separately reported measurement-monitor current.
    area = np.pi*(.002)**2
    currents = np.array([60e-6, 240e-6])
    js = currents/area
    volts = abs(offsets[0]['potential_V'])
    fields = volts/30e-9
    charging = dict(conditioning_I_10keV_A=60e-6, spot_diameter_m=.004,
                    equal_monitor_power_current_assumption='If I10=60uA in measurement too: I(E)=60uA*10keV/E_lab',
                    conditional_J_10_and_2p5_A_per_m2=js.tolist(),
                    required_bulk_rho_ohm_m_for_abs_voltage=(fields/js).tolist(),
                    required_R_area_ohm_m2=(volts/js).tolist(),
                    required_vertical_resistance_ohm=(volts/currents).tolist(),
                    rho_comparison_note='Vertical Ohmic steady-state demand E/J for POSITIVE ion charging; a NEGATIVE sustaining electron supply is separately required',
                    dielectric_capacitance_per_area_F_per_m2_per_epsilon_r=8.8541878128e-12/30e-9,
                    required_sheet_charge_C_per_m2_per_epsilon_r=8.8541878128e-12*fields)
    write_csv('curves.csv', curves)
    summary = dict(model=dict(truth_U_eV=40., reaction="reaction_sigma(..., reaction='eq7'): Bosch-Hale D(d,p)T",
                              beam='monoenergetic D+, no molecular admixture or unstated spread',
                              nominal_E_lab_keV=E_LAB.tolist(), reference_E_lab_keV=10.,
                              fit='equal weights in log R at 2.5,3,...,9.5 lab; exact own-profile 10keV reference; U only, no amplitude',
                              density_inputs=dict(nPd_m3=N_PD, bulk_D_per_Pd=X_BULK, nPdO_formula_m3=N_PDO),
                              oxide_stopping='ICRU/ZBL Pd and D + repository LS/ZBL O, Bragg sum; O electronic estimate; oxygen total stopping 0.5/1/2 sensitivity',
                              cutoff_lab_keV=.1, integration_lab_grid_step_keV=float(E[1]-E[0]),
                              omissions=['straggling/backscatter/channeling', 'beam-induced time-dependent chemistry and D redistribution', 'measured local D inventory or film electrical data unavailable']),
                   exposed_envelope=envelope, offsets=offsets, current_dependent_negative_control=current_control,
                   combined_sheet_potential_lowest_point_bounds=combined_bounds,
                   positive_monitor_only_lowest_point_bounds=positive_monitor_bounds,
                   charging=charging, consumer_parity_max_relative_errors=parity,
                   zr_gate=dict(executed=False, reason='No defensible Pd branch reaches the reported contrast. Mathematical negative-voltage fits do not establish a beam-driven mechanism; no Zr calculation authorized by this gate.'),
                   decision='Fixed D-rich exposed profile cannot reach 310/600 in the frozen Pd envelope; mapped feedback and in-cap/buried cases are sensitivities, not all-profile exclusion. Constant negative collision potential fits require an undocumented wrong-sign supply and cannot survive conservative return to grounded buried metal.')
    with open(HERE/'results.json', 'w') as f:
        json.dump(summary, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps(dict(envelope=envelope, offsets=offsets, current_control=current_control,
                          parity=parity, charging=charging), indent=2))


if __name__ == '__main__':
    main()
