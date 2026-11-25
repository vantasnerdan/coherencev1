"""Read saved trajectories/launches; render the integrated result without simulation."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import engine as E
import rate_mapping as R

HERE=Path(__file__).resolve().parent
PATTERNS=('converging_wave','reverse_spreading','collision','direct_optimized')
LABELS=('Waves','Reverse (control)','Collision (8 D)','Direct')


def fmt(x):return f'{x:.3g}'


def analyze():
    data=json.loads((HERE/'results.json').read_text())
    assert data['complete']
    launches=np.load(HERE/'launches.npz')
    lattices={n:E.Lattice(m) for n,m in ((32,2),(256,4),(2048,8))}
    target=[r for r in data['rows'] if r['scope']=='target']
    for row in target:
        key=f"n{row['n_D']}_e{row['energy_per_D_eV']:g}_{row['pattern']}"
        q,v=launches[key+'_q'],launches[key+'_v']
        diag=E.launch_energy_diagnostics(lattices[row['n_D']],q,v)
        eligible=diag['participation_positive_site_energy']>=row['n_D']/4 and diag['target_positive_site_energy_fraction']<.05
        for r in data['rows']+data['candidates']:
            if all(r[k]==row[k] for k in ('n_D','energy_per_D_eV','pattern','candidate')):
                r.update(diag,headline_eligible=bool(eligible))
    spectra=R.load_contact_spectra()
    data['conditional_readouts']=[dict(n_D=r['n_D'],energy_per_D_eV=r['energy_per_D_eV'],pattern=r['pattern'],readout=R.conditional_readout(r,spectra)) for r in target]
    data['scalings']=[dict(energy_per_D_eV=e,pattern=p,result=R.conditional_size_scaling([r for r in target if r['energy_per_D_eV']==e and r['pattern']==p and r['headline_eligible']]))
                      for e in (.0006,.0048,.025) for p in ('converging_wave','direct_optimized')]
    data['conditional_energy_thresholds']=[
        dict(n_D=2048,anchor_energy_per_D_eV=r['energy_per_D_eV'],anchor_pair_gain=r['pair_gain'],
             assumption='hold this finite-case gain fixed while increasing e; not demonstrated',
             energy_per_D_eV_for_100_eV=100/r['pair_gain'],energy_per_D_eV_for_1_keV=1000/r['pair_gain'])
        for r in target if r['n_D']==2048 and r['pattern']=='direct_optimized']
    if (HERE/'consumer_proof.json').exists():
        proof=json.loads((HERE/'consumer_proof.json').read_text())
        data['selected_step_refinement']=proof['single_limiting_step_refinement']
        data['refined_selected_row']=proof['refined_selected_row']
        data['refined_conditional_readout']=R.conditional_readout(proof['refined_selected_row'],spectra)
        data['numerical_status']='Finite-step comparisons; largest-energy peak materially changes under the single .5/.25-fs refinement. No convergence claim.'
        data['refinement_trace_columns']=data['trace_columns'][:10]+['most_displaced_site_local_V_eV','most_displaced_site_harmonic_V_eV']
        refined=proof['refined_selected_row']
        data['refined_threshold_arithmetic']=dict(
            assumption='same-launch .25-fs anchor efficiency or gain held fixed; not a prediction',
            n_D_for_100_eV=100/(refined['fraction_total_energy_at_peak']*.025),
            n_D_for_1_keV=1000/(refined['fraction_total_energy_at_peak']*.025),
            energy_per_D_eV_for_100_eV=100/refined['pair_gain'],
            energy_per_D_eV_for_1_keV=1000/refined['pair_gain'])
        for series in data['scalings']:
            if series['energy_per_D_eV']==.025:
                series['result']['resolution_limit']='The N2048 direct anchor changes materially at .25 fs; no converged high-energy size law is established.'
    data['launch_energy_eligibility']='Positive K_i+V_080 site-energy participation >=N/4 and target site-energy fraction<.05; signed interactions excluded.'
    data['candidate_counts_by_case']=[dict(n_D=n,energy_per_D_eV=e,count=sum(r['n_D']==n and r['energy_per_D_eV']==e for r in data['candidates'])) for n in (32,256,2048) for e in (.0006,.0048,.025)]
    (HERE/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    fields=list(data['rows'][0])
    for r in data['rows']:
        for k in r:
            if k not in fields:fields.append(k)
    with (HERE/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(data['rows'])
    def row(n,e,p):return next(r for r in target if r['n_D']==n and r['energy_per_D_eV']==e and r['pattern']==p)
    best=row(2048,.025,'direct_optimized')
    fine=data.get('refined_selected_row',best)
    reverse=[r for r in target if r['pattern']=='reverse_spreading']
    collisions=[r for r in target if r['pattern']=='collision']
    controls=data['controls']
    hlow,hhigh=[r for r in controls if r['pattern']=='matched_harmonic']
    loss=next(r for r in controls if r['pattern']=='loss_sensitivity')
    lines=['# Can a distributed bulk PdD launch focus onto one D–D pair?','',
           '**Finite-model result, unpublished.** The target is the actual D25→D63 [110] bond; every launch uses the same total N e. Wave/direct launches start with zero target velocity and zero displacement. Reverse recovery and localized collisions are shown, but are not silently promoted to collective preparation.','',
           f"For the same distributed direct launch at N=2048, 25 meV/D, the two resolutions give **{fine['peak_relative_kinetic_eV']:.2f}–{best['peak_relative_kinetic_eV']:.2f} eV**, both above the wave's {row(2048,.025,'converging_wave')['peak_relative_kinetic_eV']:.3f} eV. This two-result spread is **not an error bar** or converged peak. Only {100*fine['fraction_total_energy_at_peak']:.2f}–{100*best['fraction_total_energy_at_peak']:.2f}% of the 51.2-eV budget reaches the target; no common incoming state or saturation is established.",'',
           '## Same-budget comparison','',
           'At **dt=.5 fs**, each cell is **gain K_rel/e; contiguous half-peak hold (fs)**. † is initially preconcentrated/not headline-eligible; ‡ marks an edge-censored hold. The final direct anchor has the visible .25-fs sensitivity below. Absolute energies, target radial component, separation, crossings and same-trajectory whole-cell maximum are in [comparison.csv](comparison.csv).','',
           '| D count | meV/D | Waves | Reverse | Collision | Direct |','|---:|---:|---:|---:|---:|---:|']
    for n in (32,256,2048):
        for e in (.0006,.0048,.025):
            cells=[]
            for p in PATTERNS:
                r=row(n,e,p);mark=('†' if not r['headline_eligible'] else '')+('‡' if r['hold_left_censored'] or r['hold_right_censored'] else '')
                cells.append(f"{r['pair_gain']:.3g}; {r['half_peak_hold_fs']:.1f}{mark}")
            lines.append(f'| {n} | {e*1000:g} | '+' | '.join(cells)+' |')
    lines += ['',f"Direct launch positive-site-energy participation is {fmt(best['participation_positive_site_energy'])}/{best['n_D']}, target preload {100*best['target_positive_site_energy_fraction']:.3g}%. Collision participation is {fmt(collisions[-1]['participation_positive_site_energy'])}: its eight-site packet gets more energetic as N grows, not more extensive. Reverse states have positive-site participation {min(r['participation_positive_site_energy'] for r in reverse):.3g}–{max(r['participation_positive_site_energy'] for r in reverse):.3g} and target site-energy fraction {min(r['target_positive_site_energy_fraction'] for r in reverse)*100:.2g}–{max(r['target_positive_site_energy_fraction'] for r in reverse)*100:.2g}%; exact recovery is a reachability control wherever that launch is poorly spread.",'',
              '## What limits the result?','']
    exits=[r for r in target if r['first_site_or_bond_exit_fs'] is not None]
    lines.append(f"{len(exits)}/36 selected trajectories hit the .80-Å site or 1.287824-Å nearest-bond boundary; every quoted peak is strictly before any exit, and censored holds remain labelled. The same supported-prefix rule is used for every family. The model's extent is a box of assumed support, not measured angular coverage.")
    lines += ['',f"At N=2048/.6 meV the same-launch harmonic peak is {fmt(hlow['peak_relative_kinetic_eV'])} eV before its site-boundary exit at {fmt(hlow['first_site_or_bond_exit_fs'])} fs, versus {fmt(row(2048,.0006,'direct_optimized')['peak_relative_kinetic_eV'])} eV for the nonlinear trajectory. At 25 meV the harmonic trajectory exits at {fmt(hhigh['first_site_or_bond_exit_fs'])} fs, long before the nonlinear peak: its {fmt(hhigh['peak_relative_kinetic_eV'])}-eV supported prefix is **clipping, not a same-horizon nonlinear advantage**. The one fixed 4136-fs amplitude-loss sensitivity changes the coarse high-energy direct peak to {fmt(loss['peak_relative_kinetic_eV'])} eV ({100*(loss['peak_relative_kinetic_eV']/best['peak_relative_kinetic_eV']-1):+.1f}%). No damping sweep was used.",'']
    proof_path=HERE/'consumer_proof.json'
    if proof_path.exists():
        proof=json.loads(proof_path.read_text());ref=proof['single_limiting_step_refinement']
        lines.append(f"The same launch at .25 fs peaks at **{ref['fine_peak_eV']:.3g} eV** at {ref['fine_peak_time_fs']:g} fs, versus {ref['coarse_peak_eV']:.3g} eV at {ref['coarse_peak_time_fs']:g} fs at .5 fs: {100*ref['fractional_peak_change']:+.1f}%, so **the largest-energy peak is not numerically converged**. Holds are {ref['coarse_hold_fs']:g}/{ref['fine_hold_fs']:g} fs; energy drift is {ref['coarse_relative_energy_drift']:.2g}/{ref['fine_relative_energy_drift']:.2g}. Both exact same-launch rows and their pre-peak exits are retained. "+proof.get('wall_interpretation',''))
        lines += ['', '| Same launch, no reoptimization | .5 fs | .25 fs |','|---|---:|---:|',
            f"| Peak Krel = radial K (eV) | {best['peak_relative_kinetic_eV']:.10f} | {fine['peak_relative_kinetic_eV']:.10f} |",
            f"| Peak time / half-peak hold (fs) | {best['peak_time_fs']:g} / {best['half_peak_hold_fs']:g} | {fine['peak_time_fs']:g} / {fine['half_peak_hold_fs']:g} |",
            f"| Signed inward speed (Å/fs) | {best['signed_inward_speed_A_per_fs_at_peak']:+.8f} (inward) | {fine['signed_inward_speed_A_per_fs_at_peak']:+.8f} (receding) |",
            '| First support exit | none by 2400 fs | none by 2400 fs |','',
            'This is late-pulse/peak-selection sensitivity, **not changed clipping**: the fine trace near the coarse peak (±30 fs) reaches only '+fmt(ref['fine_peak_in_coarse_peak_window_eV'])+' eV. The two opposite signs must not be treated as the same incoming state.']
    lines += ['', '## Energy floors first; extrapolations second','',
        'Perfect-concentration floors for 100 eV / 1 keV are N ≥ **166,667 / 1,666,667** at .6 meV/D, **20,834 / 208,334** at 4.8 meV/D and **4,000 / 40,000** at 25 meV/D. At N=2048 the per-D floors are **48.83 / 488.28 meV**, before inefficiency. Neither target is simulated.','',
        '**Scaling is finite and preparation-dependent, not an asymptotic limit.** Direct gain slopes are between 32→256 and 256→2048 D at .5 fs. The 25-meV anchor is not converged; the table exposes both timestep-based arithmetic assumptions rather than claiming a fitted high-energy size law.','',
        '| anchor meV/D | two coarse finite-size exponents | 100 eV / 1 keV N if efficiency stays fixed | eV/D at N=2048 if gain stays fixed |','|---:|---:|---:|---:|']
    for e in (.0006,.0048,.025):
        ss=next(s['result'] for s in data['scalings'] if s['pattern']=='direct_optimized' and s['energy_per_D_eV']==e)
        fit=ss['power_law_scenario'];exponents=fit.get('local_exponents',[])
        eta=ss['constant_efficiency_scenario'];thresholds=eta['targets'] if eta else []
        gain=row(2048,e,'direct_optimized')['pair_gain']
        ntext=' / '.join(f"{t['n_D_if_constant_efficiency']:.3g}" for t in thresholds) if thresholds else 'unresolved'
        etext=f'{100/gain:.3g} / {1000/gain:.3g}'
        if e==.025 and 'refined_threshold_arithmetic' in data:
            a=data['refined_threshold_arithmetic']
            ntext+=f" (.5 fs); {a['n_D_for_100_eV']:.3g} / {a['n_D_for_1_keV']:.3g} (.25 fs)"
            etext+=f" (.5 fs); {a['energy_per_D_eV_for_100_eV']:.3g} / {a['energy_per_D_eV_for_1_keV']:.3g} (.25 fs)"
        lines.append(f"| {1000*e:g} | {', '.join(fmt(a) for a in exponents) or fit['status']} | {ntext} | {etext} |")
    lines += ['', 'Those N and e values are **assumption-only extrapolations**, not simulated reachable states; measured efficiency is already energy-dependent, and the largest-energy step sensitivity limits its numerical anchor. `results.json` also retains the finite-size power-law scenarios and their incompatible slopes/residuals. No three-point flattening or growth establishes saturation or its absence. The fixed 2400-fs horizon, nonlinear wall/defocusing and model exits remain distinct.','',
              '## Conditional quantum contact and G12','',
              'All 36 target peaks have a [retained per-level lookup](results.json), with both joins and both original flagged-row treatments. The illustrative assignment is **ΔE=classical K_rel above each quantum model’s own ground**, not equality to an eigenenergy or an inferred quantum state. Adjacent l=0 eigenlevels are separately reported; their endpoint span is not a rigorous wavepacket bound.']
    rd=next(r['readout'] for r in data['conditional_readouts'] if r['n_D']==2048 and r['energy_per_D_eV']==.025 and r['pattern']=='direct_optimized')
    fine_rd=data.get('refined_conditional_readout',rd)
    if rd['status']=='supported' and fine_rd['status']=='supported':
        spans=rd['across_join_endpoint_span'];lo=spans['lower'];hi=spans['upper']
        flo=fine_rd['across_join_endpoint_span']['lower'];fhi=fine_rd['across_join_endpoint_span']['upper']
        lines += ['', '| Base-10 logs, stationary-level assignment only | .5-fs peak: lower / upper endpoints | .25-fs peak: lower / upper endpoints |','|---|---:|---:|']
        for key,label in [('log10_contact_A_minus3','contact (Å⁻³)'),('log10_emitted_p_per_s','emitted p/s'),('log10_emitted_p_per_true_C','emitted p/true C'),('log10_delta_R_k0','PdD .25-keV k0'),('log10_delta_R_k1','PdD .25-keV k1')]:
            lines.append(f"| {label} | [{lo[key][0]:.2f},{lo[key][1]:.2f}] / [{hi[key][0]:.2f},{hi[key][1]:.2f}] | [{flo[key][0]:.2f},{flo[key][1]:.2f}] / [{fhi[key][0]:.2f},{fhi[key][1]:.2f}] |")
    else:lines += ['', 'A peak exceeds one or more retained contact spectra; its readout is unresolved, not zero.']
    lines += ['', 'These are **one continuously active normalized-pair** units, not detector counts, a foil rate or a peak×hold event yield. Preparation probabilities, radial populations and coherence are unknown. The retained excitation support ends at 7.906–7.969 eV; no 100-eV/1-keV contact or rate is extrapolated.','',
              '**Three largest omissions:** (1) the off-ray local/connected angular force surface is not measured; (2) dynamic Pd/electrons, finite-temperature disorder and physical damping are absent except the stated one sensitivity; (3) generating the prescribed phases and converting a classical focus into occupied quantum radial states is not supplied.','',
              'The closest connected alternate-H difference is −0.0047422856 eV, with signed slope difference +0.0413149205 eV/Å; this finite-cutoff comparison is not a PAW error bound.','',
              '![Finite-size gains and selected target histories](focusing.png)','',
              'Equations, selection, support, resource controls and portable commands: [method.md](method.md). Exact staged capsule bytes and whitelist: `capsule.json`, `ARCHIVE_FILES.txt`, `SHA256SUMS`; prior-input assumptions: `input_manifest.json`. No prior contact or .80-well calculation was rerun.']
    (HERE/'README.md').write_text('\n'.join(lines)+'\n')
    traces=np.load(HERE/'trajectories.npz')
    fig,axes=plt.subplots(2,3,figsize=(12,6.8),constrained_layout=True)
    colors=dict(zip(PATTERNS,('C0','C1','C2','C3')))
    for col,e in enumerate((.0006,.0048,.025)):
        for p,label in zip(PATTERNS,LABELS):
            rr=[row(n,e,p) for n in (32,256,2048)]
            axes[0,col].loglog([r['n_D'] for r in rr],[r['pair_gain'] for r in rr],'o-',color=colors[p],label=label)
            key=f'n2048_e{e:g}_{p}';tr=traces[key]
            valid=tr[:,0]<(row(2048,e,p)['first_site_or_bond_exit_fs'] or np.inf)
            if p!='reverse_spreading':
                axes[1,col].plot(tr[valid,0],tr[valid,1],color=colors[p],label=label)
        axes[0,col].set_title(f'{e*1000:g} meV/D')
        axes[0,col].set_xlabel('D count');axes[0,col].set_ylabel('Target gain Krel/e')
        axes[1,col].set_xlabel('Time (fs)');axes[1,col].set_ylabel('Target Krel (eV), N=2048')
        axes[0,col].grid(alpha=.2);axes[1,col].grid(alpha=.2)
    if 'n2048_e0.025_direct_optimized_dt025' in traces.files:
        fine_trace=traces['n2048_e0.025_direct_optimized_dt025']
        axes[1,2].plot(fine_trace[:,0],fine_trace[:,1],'--',color='black',label='Same direct, .25 fs')
        axes[0,2].plot([2048],[fine['pair_gain']],'*',color='black',ms=10)
        axes[1,2].legend(fontsize=7)
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    fig.savefig(HERE/'focusing.png',dpi=160);plt.close(fig)
    print(json.dumps({'target_rows':36,'largest_direct':best,'selected_exits':len(exits),'controls':controls},indent=2))

if __name__=='__main__':analyze()
