"""Write the G16 supplement in the existing G12 five-column ledger format."""
import argparse
import hashlib
import json
from pathlib import Path


PATH = 'g16_hot_electrons'


def rows(readout, readout_sha, join, join_sha):
    table = []
    native = f'`{PATH}/response.json`, `{PATH}/response.csv`'
    binding = f'response JSON SHA-256 `{readout_sha}`'
    join_src = f'`{PATH}/track-response.json`'
    join_binding = f'track-response JSON SHA-256 `{join_sha}`; G15 source `{join["track_source"]["sha256"]}`'
    base = 'Fixed-ion neutral rocksalt PdD; a=4.084 Å, PW500 PBE PAW, Fermi–Dirac kBT; no fit to Davis.'
    for r in readout['rows']:
        t = r['temperature_K']
        check = r['mesh_check']
        mesh = 'two meshes agree at the stated tolerance' if check['delta_K_mesh_agrees'] else 'mesh difference exceeds the stated tolerance'
        for m in (12, 16):
            v = r[f'mesh{m}']
            b = v['contrast']
            z = v['zero']
            setup = f'{base} Γ-centered {m}³ mesh, {t} K.'
            model = v['harmonic']['K_force_eV_A2']
            table.append((f'PdD Γ [100] stiffness: {model:.8g} eV/Å²; ΔK {b["K_change_eV_A2"]:+.7g} eV/Å² at {t} K, k{m}', native,
                          binding, f'Common 300 K zero, fixed Pd, ±0.01 Å. ΔK {mesh}. Collective, not isolated on-site K.', setup))
            table.append((f'PdD electronic excitation ΔU {b["electronic_excitation_eV_per_D_site"]:+.8g} eV/D site at {t} K, k{m}', native,
                          binding, 'Internal energy relative to 300 K per formula unit, not the free energy; retain mesh sensitivity.', setup))
            table.append((f'PdD pressure {z["pressure_GPa"]:+.8g} GPa; ΔP {b["pressure_change_GPa"]:+.8g} GPa at {t} K, k{m}', native,
                          binding, 'Pressure from stress at the undisplaced fixed volume; compare both meshes, not a relaxed lattice constant.', setup))
            d = b['centered_D_force_change_eV_A']
            table.append((f'PdD centered ΔF_D = ({d[0]:+.2g}, {d[1]:+.2g}, {d[2]:+.2g}) eV/Å at {t} K, k{m}', native,
                          binding, 'Raw residual at the inversion center, not a resolved directed force; the signed displacement response is in the native arrays.', setup))
            for axis, q in v['quartic'].items():
                fit = q['fit']
                table.append((f'PdD collective [{axis}] quartic: {fit["K4_force_fit_eV_A4"]:+.8g} eV/Å⁴ at {t} K, k{m}', native,
                              binding, 'Two-amplitude ±0.12/0.24 Å force fit with energy holdout in response JSON; no separate mesh control on K4.', setup))
    for label, result in join['modelled_peaks'].items():
        use = 'Sign above the observed numerical scale.' if result.get('sign_above_floor_on_both_meshes') else 'Do not select an interpolated sign as a material input.'
        table.append((f'PdD G15 {label.replace("_", " ")} peak: {result["excitation_eV_per_D_site"]:.8g} eV per model D site',
                      join_src, join_binding, f'Compare matched ΔU on both meshes; retain original bracket/status. {use}',
                      'Conditional carrier-model density × fixed PdD volume, paired with independent DFT; not fitted to Davis.'))
    return table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('readout', type=Path)
    p.add_argument('track_join', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    response = a.readout.read_bytes()
    track = a.track_join.read_bytes()
    h = hashlib.sha256(response).hexdigest()
    join = json.loads(track)
    if join['readout_sha256'] != h:
        raise ValueError('Track join and ledger refer to different native response bytes')
    with a.output.open('x') as out:
        out.write('# G16 proposed rows for `g12_residual/inputs_ledger.md`\n\n'
                  '| Input / adopted number or table | Source file | Revision / byte binding | Latest or superseded; what changes | Basis; fitted to these 48? |\n'
                  '|---|---|---|---|---|\n')
        for row in rows(json.loads(response), h, join, hashlib.sha256(track).hexdigest()):
            out.write('| ' + ' | '.join(row) + ' |\n')


if __name__ == '__main__':
    main()
