"""Build the human readout from actual saved site and pair artifacts, no replay."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from solve import OUT
from well import HERE


def read(path):
    return json.loads(path.read_text())


def main():
    rows = read(OUT/'single_site.json')
    transitions = read(OUT/'transitions_100.json')
    controls = read(OUT/'numerical_controls.json')
    source = read(OUT/'source_evidence.json')
    consumers = read(OUT/'consumer_evidence.json')
    pair = read(HERE/'pair/result.json')
    integration_evidence(rows)
    cuts = read(HERE/'inputs/cuts.json')['directions']
    anchor_lines = ['| r (Å) | [100] symmetry mean (eV) | [110] symmetry mean (eV) | +[111] (eV) | −[111] (eV) |', '|---|---:|---:|---:|---:|']
    for i, radius in enumerate((.12, .24, .48, .8)):
        energies = {key: [p['delta_energy_free_eV'] for p in cuts[key][i]['signed_points']] for key in cuts}
        anchor_lines.append(f"| {radius:.2f} | {np.mean(energies['100']):.6f} | {np.mean(energies['110']):.6f} | {energies['111'][1]:.6f} | {energies['111'][0]:.6f} |")
    lines = ['# TiD₂: a stiff, tetrahedral single-site well and a conditional pair comparison', '',
        'This capsule uses the supplied **constrained cubic, fully occupied fluorite** electronic calculation. It does not establish the phase, loading or surface of a titanium foil. H/D/T change only the nuclear mass in the same frozen-neighbour potential; the source static calculation used D. Energies below are eV (or explicitly meV), position Å, wave number Å⁻¹, and momentum p=ħk. The electronic energy is the force-consistent free energy per simulated cell; one-site energy rises are for one displaced nucleus.', '',
        '## What the source actually says', '',
        'Atom 8 is at the +tetrahedral site in conventional cubic x=[100], y=[010], z=[001] coordinates. +[111] moves away from one nearest Ti and toward the other three; −[111] moves toward that one Ti. At r=0.80 Å their electronic rises are respectively **0.584880 and 5.294128 eV**. This ninefold lopsidedness is retained, not averaged away. [100] softens as displacement grows; [110] is intermediate. These are computed frozen-ion rays, not measured three-dimensional spectra.', '',
        '\n'.join(anchor_lines), '',
        f"The true [100]/[110] twofold signed residuals are small: max energy differences {source['symmetry']['100']['max_signed_energy_residual_eV']:.3g}/{source['symmetry']['110']['max_signed_energy_residual_eV']:.3g} eV and rotated-force differences {source['symmetry']['100']['max_rotated_force_residual_eV_A']:.3g}/{source['symmetry']['110']['max_rotated_force_residual_eV_A']:.3g} eV/Å. Only those symmetry-equivalent rays are averaged. All 24 selected point bytes, displacement signs, energy rises and projected forces were checked against their native records. The mesh comparison is **4→6 at zero and −0.48 Å [111] only**, changing that rise by 0.002392 eV; it is not a whole-surface or pair error estimate.", '',
        'The regular T_d interpolant and matched synthetic exterior were declared before implementation in [model.md](model.md). It reproduces ray energies and slopes, includes the allowed odd xyz tetrahedral term, and is polynomial at the origin. Few rays do **not** identify a full angular surface. These are localized-well distributions conditional on that closure, not periodic diffusion among tetrahedral/interstitial sites. The oscillator basis is on all space, not a hard box.', '',
        '## Finite-model levels and the two [100] transitions', '',
        'The three-dimensional Hamiltonian is −ħ²∇²/(2m)+V. “Exact diagonalization” means the declared finite projected model, not an exact measured 3D potential. Twenty oscillator functions per Cartesian coordinate and order-40 Gauss-Hermite potential quadrature are used; the eight coordinate-parity sectors couple in tetrahedral pairs through xyz, rather than incorrectly imposing inversion symmetry. Full canonical multiplets are included.', '',
        '| isotope | 3D ground energy (meV) | first 3D gap (meV) | next 3D excitation clusters (meV; multiplicity) |',
        '|---|---:|---:|---|']
    for isotope in ('H', 'D', 'T'):
        row = next(r for r in rows if r['isotope'] == isotope)
        levels = row['levels']
        later = ', '.join(f"{r['excitation_eV']*1000:.3f} ({r['multiplicity']})" for r in levels[2:5])
        lines.append(f"| {isotope} | {levels[0]['energy_eV']*1000:.3f} | {levels[1]['excitation_eV']*1000:.3f} (triplet) | {later} |")
    lines += ['', 'Observed low-cluster multiplicities are 1,3,1,3,2; no inversion-parity labels are assigned. Both temperatures use these same levels with different canonical weights. Higher low levels and adjacent cluster gaps are in `output/single_site.json`.', '',
        '| isotope | 1D [100] E₁−E₀ (meV) | E₂−E₁ (meV) | χ=(E₂−E₁)−(E₁−E₀) (meV/quantum) |', '|---|---:|---:|---:|']
    for row in transitions:
        r = row['result']
        lines.append(f"| {row['isotope']} | {1000*r['first_gap_eV']:.3f} | {1000*r['second_gap_eV']:.3f} | {1000*r['shift_per_quantum_eV']:.3f} |")
    lines += ['', 'This is exactly the G13 convention: successive eigenlevel spacings of the **one-dimensional neighbour ray**, not the 3D triplet gap or a thermal stiffness. The PdD values 42.647/47.989 meV belong to its **[110]** ray; they are provenance, not Ti values. The 0.002→0.001 Å spacing comparison changes Ti gaps by at most 0.004 meV and leaves all isotope softening signs unchanged.', '',
        'For a carrier k and sideband p, define d=[ε(k+p)+ε(k−p)]/2−ε(k). Before damping, instability requires **d(d+2χn)<0**. Ti’s negative χ is compatible with positive d (positive local curvature at small p), with 0<d<2|χ|n in that model. The source-derived single-site sign does not locate a band, determine its bandwidth/coupling, give a growth rate, or establish a damped threshold. No envelope simulation was run.', '',
        '## Full distributions before widths', '',
        '![Spatial, momentum and signed distributions](distributions.png)', '',
        'Each `*_position.csv` gives angle-averaged density ρ̄(r) [Å⁻³], radial probability measure 4πr²ρ̄ [Å⁻¹], and integrated survival P(|x|>r). Each `*_momentum.csv` analogously gives n̄(k) [Å³], 4πk²n̄ [Å], and survival; it also gives physical SI p=ħk, density with the (ħ·10¹⁰)⁻³ Jacobian and radial measure with the (ħ·10¹⁰)⁻¹ Jacobian. Both position and momentum are derived from the same normalized canonical eigenstates, using Fourier phases (−i)^n, not mass-relabeling a position curve.', '',
        '**Own-peak semantics:** density and radial-measure cutoffs are independently the first outward post-peak crossing of 10⁻⁶ of that same curve’s own peak. These are different radii and different from “probability 10⁻⁶.” Integrated survival at each crossing is retained, not set equal to 10⁻⁶. The table below uses the density crossings. Synthetic support begins at r>0.8 Å in position only; there is no corresponding k-space measurement boundary.', '',
        '| isotope | K | per-axis MSD (Å²) | r at density 10⁻⁶ (Å) | P(r>r₆) | k at density 10⁻⁶ (Å⁻¹) | P(k>k₆) | P(r>0.8 Å) |', '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['isotope']} | {r['temperature_K']} | {r['mean_x2_A2']:.6f} | {r['position_density_1e6_A']:.4f} | {r['position_density_tail_at_1e6']:.3g} | {r['momentum_density_1e6_A_inverse']:.3f} | {r['momentum_density_tail_at_1e6']:.3g} | {r['probability_r_gt_080']:.3g} |")
    lines += ['', 'T_d makes all three canonical Cartesian MSDs equal (reported separately with their observed anisotropy ratio). It does not make the two [111] sides equivalent: `*_g9_plus111.csv` and `*_g9_minus111.csv` preserve the distinct signed on-axis cuts **and** actual normalized plane-integrated marginals. A line density ρ(s n) has units Å⁻³ and is not a 1D probability density. The marginal ∫ρ(s n+t a+u b)dtdu has units Å⁻¹. Their skew shapes can differ while the vector mean remains zero. Real stationary eigenstates have even momentum occupation under k→−k even though position is not inversion symmetric.', '',
        '### What the numerical and closure comparisons change', '']
    basis = [r for r in controls if r['comparison'].startswith('basis')]
    maxlevel = max(abs(v) for r in basis for v in r['low_cluster_excitation_delta_eV'])
    lines.append(f"Basis 16→20 changes the first six excitation clusters by at most {maxlevel*1000:.4g} meV. `output/numerical_controls.json` gives per-isotope/temperature MSD, kinetic, crossing and exterior-tail changes, and the absolute curve difference at each curve’s own 10⁻⁶ crossing. These are empirical numerical differences, not rigorous enclosures. The selected finite-model tolerances are 0.1 meV for low clusters, 1% for MSD/mean kinetic and 0.02 Å / 0.5 Å⁻¹ for density crossings; sub-threshold curve differences are shown directly, including any unresolved crossing.")
    for r in controls:
        if r['comparison'].startswith('basis'):
            continue
        d = r['deltas']
        lines.append(f"- H 300 K {r['comparison']}: ΔMSD={d['mean_x2_A2']:.3g} Å², Δr₆={d['position_density_1e6_A']:.3g} Å, Δk₆={d['momentum_density_1e6_A_inverse']:.3g} Å⁻¹, ΔP(r>0.8)={d['probability_r_gt_080']:.3g}.")
    worst = max(v['delta_over_threshold'] for r in basis for v in r['tail_comparison'].values())
    lines.append(f"The largest basis-difference/own-peak-threshold ratio at these crossings is {worst:.3g}; no unmeasured angular or electronic error is included in that ratio. The H 300 K density-grid/angular comparison separately prices interpolation/cubature. The ray-preserving angular penalty plus stiffer exterior is a **model sensitivity**, not a confidence interval: it changes unmeasured support while preserving all ray anchors and slopes. Where that change exceeds finite-model numerical differences, interpretation follows the closure sensitivity rather than extra solver digits.")
    lines.append('This comparison changes the H 300 K MSD by −3.9%, moves the spatial density crossing inward by 0.035 Å and the momentum crossing outward by 0.802 Å⁻¹, and reduces the already small synthetic-exterior probability from 1.01×10⁻⁷ to 2.03×10⁻⁸. The localization conclusion survives, but unique 3D tail digits are not identified by the rays. The main reported numbers remain the declared base closure, not an average of closures.')
    lines += ['', '## G8/G9 handoff actually exercised', '',
        '`output/*_g8_joint.json` contains the real G8 state dictionaries: `excitation_eV`, `kinetic_eV`, `probabilities`. Joint probabilities include the canonical level weight and normalized Fourier kinetic quadrature conditional on that eigenstate; excitation and kinetic energy are not assumed independent. Each isotope uses its own mass. Kinetic bins retain probability and conditional mean; angular correlations are deliberately isotropized for this consumer.', '',
        consumers['g8']['spectator_envelope'], '',
        'The actual `BoundSpectatorKernel.fold` read the D 300 K table with the already existing **TiD** profile as a numerical-interface exercise and returned a finite positive value; it was not relabeled as a TiD₂ rate or compared with foil charts. H/T state measures are supplied, but their unequal-mass collision physics was not invented inside this D-D reader. `output/consumer_evidence.json` records six normalizations, conditional kinetic checks and bound source hashes. Existing PdD tables are untouched.', '',
        'G9 has no crystallographic density-CSV API. Its applicable signed ballistic convention was exercised with the actual module constants: s=x·n, z_target=Z₀+s, Δt=s/v, n=+[111]/√3. The opposite normal reflects s and the density; `output/g9_mapping.csv` carries both mappings. This is a registered bulk crystal direction, **not an assigned Ti foil surface**.', '',
        '## Frozen nearest pair, beside PdD', '',
        '![Closing work and conditional contact](pair-comparison.png)', '',
        'The integrated pair files are consumed directly from [pair/README.md](pair/README.md) and `pair/result.json`; unchanged PdD evidence was reused, not rerun. Pair closing uses u per site and R=R₀−2u. The connected JSON force is −dW/du, distinct from the all-atom mixed-force response. All four measured pair rows, including the two shortest with unknown PAW-overlap behavior, remain in the main result. No dedicated pair mesh control exists. The contact comparison is the existing pure l=0 radial general-pair model with a common Coulomb core and **borrowed PdD electronic inner shape**, not a Ti nuclear rate or measured Ti inner potential.', '',
        '| separation R (Å) | native Ti closing cost (eV) | native connected mutual W (eV) |', '|---|---:|---:|']
    for r in pair['measured']['native_Ti_rows']:
        lines.append(f"| {r['R_A']:.6f} | {r['total_closing_work_eV']:.6f} | {r['connected_W_eV']:.6f} |")
    lines += ['', '| common join R (Å) | preparation | Ti log₁₀ρ(0) [Å⁻³] | Pd log₁₀ρ(0) [Å⁻³] | log₁₀(Ti/Pd) |', '|---|---|---:|---:|---:|']
    for r in pair['primary_comparison']:
        lines.append(f"| {r['join_A']:.6f} | {r['preparation']} | {r['Ti_log10_rho_A_minus3']:.5f} | {r['Pd_log10_rho_A_minus3']:.5f} | {r['log10_Ti_over_Pd']:.5f} |")
    lines += ['', 'Only the frozen longitudinal nearest-pair line is supplied. It is **not** a full coupling matrix, phonon band, relaxation response or lattice bath, and no PdD J is transferred to Ti.', '',
        '## Three dominant omissions', '',
        '1. The full angular and exterior energy surface and the actual foil phase/surface are unknown. The few signed rays support a conditional localized model, not an identified diffusion landscape.',
        '2. Neighbours and electrons do not respond dynamically: no coupled-lattice thermal state, electron-phonon damping, real-time response or change in the well during an ion’s arrival. Only the frozen nearest-pair closing path is present.',
        '3. Pair contact lacks a Ti electronic inner cut and uses a common borrowed PdD interior; the pair has no independent mesh control and the shortest two points have no PAW-overlap check. Contact values express this conditional outer-wall sensitivity, not a material yield.', '',
        '## Files and rerun', '',
        '- `inputs/`: selected native source bytes, signed ray table and mesh comparison; `output/source_evidence.json` binds their identities and signs.',
        '- `well.py`, `model.md`: declared regular interpolation and synthetic continuation.',
        '- `solve.py`, `output/single_site.json`, `output/*csv`, `output/*g8_joint.json`: isotope spectra, distributions and native reader states; saved eigenvectors allow density reuse.',
        '- `output/numerical_controls.json`, `output/*resolution.csv`, `output/*closure.csv`: deciding basis, quadrature, observable-grid and model-sensitivity comparisons.',
        '- `pair/`: actual G11 consumer, native-register evidence, reused PdD reference and conditional contact comparison.',
        '- `manifest.json`, `SHA256SUMS`: only this directory, accepted base `463910dcf13ce69d74b138a21d53183129a779de`; no lead pages, private paths or neighboring campaigns.', '',
        'The masses are 1.00782503223, 2.01410177812 and 3.01604928199 u (the existing G10 neutral-atom isotope-mass approximation to nuclear masses; electron mass is neglected). `output/integration_evidence.json` directly checks survival integrals, beyond-crossing support and the physical momentum Jacobian, and compares G8 conditional kinetic means with the Fourier-grid means.', '',
        'From the accepted repository root, with Python 3.12+, NumPy, SciPy and Matplotlib installed:', '',
        '```sh',
        'export OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=3 MKL_NUM_THREADS=3',
        'python3 g10_quantum_site/tid2-site/solve.py produce',
        'python3 g10_quantum_site/tid2-site/solve.py controls',
        'python3 g10_quantum_site/tid2-site/consume.py',
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 g10_quantum_site/tid2-site/pair/run.py',
        'python3 g10_quantum_site/tid2-site/report.py',
        'python3 g10_quantum_site/tid2-site/package.py',
        'sha256sum -c g10_quantum_site/tid2-site/SHA256SUMS',
        '```', '',
        'The source electronic calculations, PdD bulk evidence and envelope runs are not replayed by these commands. The authenticated draft release is unpublished; its download-byte receipt is external to the frozen scientific archive. No push, publication or merge is part of this delivery.']
    (HERE/'README.md').write_text('\n'.join(lines)+'\n')
    plot_distributions(rows)
    plot_pair(pair)


def integration_evidence(rows):
    evidence = []
    for row in rows:
        name = f"{row['isotope']}_{row['temperature_K']}K"
        joint = read(OUT/f'{name}_g8_joint.json')
        mean_kinetic = sum(sum(p*k for p, k in zip(s['probabilities'], s['kinetic_eV'])) for s in joint['states'])
        item = dict(case=name, g8_minus_Fourier_grid_mean_kinetic_eV=mean_kinetic-row['mean_kinetic_eV'])
        for label in ('position', 'momentum'):
            data = np.genfromtxt(OUT/f'{name}_{label}.csv', delimiter=',', names=True)
            x, density, measure = [data[column] for column in data.dtype.names[:3]]
            tail = data['tail_probability']
            integral = np.trapezoid(measure, x)
            peak = density.max()
            crossing = row['position_density_1e6_A' if label == 'position' else 'momentum_density_1e6_A_inverse']
            item[label] = dict(radial_normalization=float(integral), survival_at_zero=float(tail[0]),
                minimum_density_over_peak=float(density.min()/peak),
                maximum_density_beyond_crossing_over_peak=float(density[x > crossing].max()/peak),
                tail_last_coordinate=float(x[-1]), survival_at_last_coordinate=float(tail[-1]),
                maximum_survival_increase=float(np.diff(tail).max()))
            if label == 'momentum':
                item[label]['physical_p_radial_normalization'] = float(np.trapezoid(data['radial_measure_p_SI_inverse'], data['p_kg_m_s']))
        evidence.append(item)
    (OUT/'integration_evidence.json').write_text(json.dumps(evidence, indent=2)+'\n')


def plot_distributions(rows):
    fig, axes = plt.subplots(4, 3, figsize=(14, 16))
    for col, isotope in enumerate(('H', 'D', 'T')):
        for temperature, color in [(77, 'tab:blue'), (300, 'tab:red')]:
            name = f'{isotope}_{temperature}K'
            x = np.genfromtxt(OUT/f'{name}_position.csv', delimiter=',', names=True)
            k = np.genfromtxt(OUT/f'{name}_momentum.csv', delimiter=',', names=True)
            s = np.genfromtxt(OUT/f'{name}_g9_plus111.csv', delimiter=',', names=True)
            axes[0, col].semilogy(x['r_A'], x['density_A_minus3']/x['density_A_minus3'].max(), color=color, label=f'{temperature} K density')
            axes[0, col].semilogy(x['r_A'], x['tail_probability'], color=color, ls='--', label=f'{temperature} K P(r>R)')
            axes[1, col].semilogy(k['k_A_inverse'], k['density_A3']/k['density_A3'].max(), color=color, label=f'{temperature} K density')
            axes[1, col].semilogy(k['k_A_inverse'], k['tail_probability'], color=color, ls='--', label=f'{temperature} K P(k>K)')
            axes[2, col].plot(s['s_A'], s['marginal_A_inverse'], color=color, label=f'{temperature} K plane marginal')
            axes[2, col].plot(s['s_A'], s['marginal_A_inverse'][::-1], color=color, ls=':', label=f'{temperature} K reflected')
            axes[3, col].semilogy(s['s_A'], s['on_axis_line_density_A_minus3'], color=color, label=f'{temperature} K on-axis')
            axes[3, col].semilogy(s['s_A'], s['on_axis_line_density_A_minus3'][::-1], color=color, ls=':', label=f'{temperature} K reflected')
        for row in (0, 1):
            axes[row, col].axhline(1e-6, color='gray', lw=.7); axes[row, col].set_ylim(1e-10, 1.3)
        axes[0, col].axvline(.8, color='gray', ls=':'); axes[0, col].set_xlim(0, 1.)
        axes[1, col].set_xlim(0, 42)
        axes[2, col].set_xlim(-.5, .5)
        axes[0, col].set_title(isotope+' in the same TiD₂ site potential')
        axes[0, col].set_xlabel('r (Å)'); axes[1, col].set_xlabel('k (Å⁻¹), p=ħk'); axes[2, col].set_xlabel('signed s along +[111] (Å)')
        axes[3, col].set_xlim(-.8, .8); axes[3, col].set_ylim(1e-5, 1e2); axes[3, col].set_xlabel('signed s along +[111] (Å)')
        for ax in axes[:, col]:
            ax.legend(fontsize=7); ax.grid(alpha=.2)
    axes[0, 0].set_ylabel('own-peak density or integrated survival')
    axes[1, 0].set_ylabel('own-peak density or integrated survival')
    axes[2, 0].set_ylabel('normalized plane marginal (Å⁻¹)')
    axes[3, 0].set_ylabel('on-axis 3D density (Å⁻³), not marginal')
    fig.tight_layout(); fig.savefig(HERE/'distributions.png', dpi=160); plt.close(fig)


def plot_pair(pair):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for label, color in [('Ti', 'tab:orange'), ('Pd', 'tab:blue')]:
        axes[1].scatter([r['join_A'] for r in pair['primary_comparison']],
                        [r[f'{label}_log10_rho_A_minus3'] for r in pair['primary_comparison']], label=label, color=color)
    axes[1].set_xlabel('common electronic join R (Å)'); axes[1].set_ylabel('log₁₀ conditional contact density (Å⁻³)')
    axes[1].legend(); axes[1].grid(alpha=.2)
    for material, key, rkey, ekey in [('Ti', 'native_Ti_rows', 'R_A', 'total_closing_work_eV'),
                                     ('Pd', 'Pd_native_rows', 'D_D_separation_A', 'total_approach_free_energy_eV')]:
        rows = pair['measured'][key]
        axes[0].plot([r[rkey] for r in rows], [r[ekey] for r in rows], 'o-', label=material+' native frozen pair')
    axes[0].invert_xaxis(); axes[0].set_xlabel('pair separation R (Å)'); axes[0].set_ylabel('closing cost from own equilibrium (eV)')
    axes[0].legend(); axes[0].grid(alpha=.2)
    fig.suptitle('All native Ti pair points retained; contact uses a common borrowed Pd electronic interior')
    fig.tight_layout(); fig.savefig(HERE/'pair-comparison.png', dpi=160); plt.close(fig)


if __name__ == '__main__':
    main()
