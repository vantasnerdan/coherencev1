"""Actual same-background neutral quartet and inherited point-minus-Gaussian mapping."""
import gc
import math
import time
import numpy as np
from ase import Atom
from ase.units import Bohr
from substrate_framework.collision_bare_reference import PeriodicNuclearGreen
from substrate_framework.collision_observables import COULOMB_EV_A
from defect_calculator import scf, save_state, write_json, require
from defect_geometry import MASS_D


def probe_stem(role,radius):
    return role+'-probe-R'+format(radius,'.3f').replace('.','p')


def quartet(atoms,relaxation,spec,radii,output):
    output.mkdir(parents=True,exist_ok=False)
    started=time.process_time()
    target_index=relaxation['target_index']
    path=relaxation['path'];target=np.asarray(path['target_A']);direction=np.asarray(path['direction'])
    require(abs(np.linalg.norm(direction)-1)<1e-12,'Nonunit actual probe direction')
    states={}
    # Geometry is standard-H PAW; NEVER mix that energy into the Hae quartet.
    atoms.calc=None;gc.collect()
    states['occupied-host']=scf(atoms,spec,0,output,'occupied-host')
    original=atoms.copy();atoms.calc=None;gc.collect()
    for role in ('occupied','vacant'):
        host=original.copy()
        if role=='vacant':
            require(host[target_index].symbol=='H','Quartet target is not D/H')
            del host[target_index]
            states['vacant-host']=scf(host,spec,0,output,'vacant-host')
            host.calc=None;gc.collect()
        probe=host.copy()
        probe+=Atom('H',position=target-float(radii[0])*direction,mass=MASS_D)
        records=[]
        for radius in radii:
            radius=float(radius)
            require(radius>0,'Nonpositive pair radius')
            probe.positions[-1]=target-radius*direction
            row=scf(probe,spec,1,output,probe_stem(role,radius))
            row['R_A']=radius;row['projectile_A']=probe.positions[-1].tolist()
            records.append(row)
        states[role+'-probes']=records
        probe.calc=None;gc.collect()
    require(abs(states['occupied-host']['electrons']-states['vacant-host']['electrons']-1)<1e-7,'Target removal changed more than one electron')
    green=PeriodicNuclearGreen(original.cell.array)
    width=np.sqrt(5.)/Bohr  # inherited qualified Hae alpha1=10 for each unit nucleus
    rows=[]
    for occupied,vacant in zip(states['occupied-probes'],states['vacant-probes']):
        radius=occupied['R_A'];require(radius==vacant['R_A'],'Unmatched BO radii')
        require(occupied['electrons']==states['occupied-host']['electrons'] and vacant['electrons']==states['vacant-host']['electrons'],'Probe incorrectly added an electron')
        raw=math.fsum([occupied['energy']['internal_energy_eV'],-vacant['energy']['internal_energy_eV'],
                       -states['occupied-host']['energy']['internal_energy_eV'],states['vacant-host']['energy']['internal_energy_eV']])
        correction=green.point_minus_gaussian_eV(-radius*direction,width)
        total=raw+correction
        rows.append(dict(R_A=radius,BO_quartet_soft_eV=raw,point_minus_gaussian_eV=correction,
                         V_pair_total_eV=total,lowering_eV=COULOMB_EV_A/radius-total))
    result=dict(schema='pd-defect-neutral-quartet-v1',case=spec['case'],specification=spec,path=path,
        rows=rows,states=states,relaxation_raw_fmax_eV_A=relaxation['raw_fmax_eV_A'],
        process_cpu_s=time.process_time()-started,material_pair_export_eligible=False,
        energy_zero='Eprobe,o-Eprobe,v-Ehost,o+Ehost,v, frozen relaxed background; target-only neutral removal.',
        scope='Current-cell static point-regularized BO; representative C h/cell controls still required. No U0/rate/dynamic observable.')
    write_json(output/'quartet.json',result)
    return result
