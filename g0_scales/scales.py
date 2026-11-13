"""Measured material inputs, one representative value per quantity.

Run: python3 scales.py   (writes scales.json)
Densities use four metal atoms per conventional FCC cell. Plasma energy is
hbar*sqrt(n e^2/(epsilon0 m)) for bare nuclear charge +e, not a screened mode.
"""
import json
from pathlib import Path
from scipy.constants import e, epsilon_0, hbar, m_p, physical_constants

ROWS = [
    dict(material="Ni/H", a_A=3.732, site="octahedral", loading=1.0, isotope="H",
         optic_meV=89.0, hop_barrier_eV=0.37,
         barrier_scope="dilute fcc Ni32H, computed O-to-T barrier",
         anharmonicity="first optic band has a shoulder to 130 meV",
         sources=dict(lattice_site_optic="Antonov et al., EPL 51, 140 (2000)",
                      hop="Di Stefano et al., PRB 92, 224301 (2015)")),
    dict(material="Ti/H (delta TiH2)", a_A=4.454, site="tetrahedral", loading=2.0, isotope="H",
         optic_meV=149.0, hop_barrier_eV=0.49,
         barrier_scope="TiH1.59 measured self-diffusion activation, 450-650 C",
         anharmonicity="TiH1.84 excitations described by a harmonic oscillator",
         sources=dict(lattice="Khidirov and Kayumova (2012), INIS 44129203",
                      site_optic="Hattori et al., PRB 106, 134309 (2022)",
                      hop="Stuhr et al., J. Less-Common Met. 172-174, 678 (1991)")),
    dict(material="Pd/D", a_A=4.084, site="octahedral", loading=1.0, isotope="D",
         optic_meV=37.0, hop_barrier_eV=None, barrier_scope=None,
         anharmonicity="shoulder to 65 meV; second band 72 meV against 2 x 37 meV",
         sources=dict(lattice_site_optic="Antonov et al., PRB 80, 134302 (2009)")),
]

rows = []
for raw in ROWS:
    mass = m_p if raw["isotope"] == "H" else physical_constants["deuteron mass"][0]
    n = 4*raw["loading"]/(raw["a_A"]*1e-10)**3
    rows.append(dict(**raw, nuclei_density_m3=n,
                     bare_ionic_plasma_eV=hbar*(n*e*e/(epsilon_0*mass))**0.5/e))
out = dict(materials=rows,
           missing="conduction-electron density, screening length, effective charge, PdD hop barrier")
Path(__file__).with_name("scales.json").write_text(json.dumps(out, indent=2) + "\n")
