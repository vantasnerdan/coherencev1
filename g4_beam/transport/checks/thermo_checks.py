"""Independent parameter-to-observable checks, with explicit thermodynamic limits."""
import numpy as np
from scipy.constants import Boltzmann, elementary_charge, gas_constant
from scipy.optimize import brentq

KB_EV = Boltzmann/elementary_charge


def dilute_sieverts(parameters, temperature, pressure_Pa):
    if not 0 < temperature < parameters['Tc'] or pressure_Pa <= 0:
        raise ValueError('Dilute Pd check needs positive pressure below its critical temperature')
    plateau_Pa = 101325*np.exp(parameters['dH']/(gas_constant*temperature)
                              -parameters['dS']/gas_constant)
    alpha = min(parameters['a300']*np.exp(-.11/KB_EV*(1/temperature-1/300)),
                parameters['xc']-.001)
    odds = alpha/(parameters['tha']-alpha)*np.sqrt(pressure_Pa/plateau_Pa)
    x = parameters['tha']*odds/(1+odds)
    if x >= alpha:
        raise ValueError('Requested pressure is outside the dilute branch')
    return dict(T_K=temperature, pressure_Pa=pressure_Pa,
                plateau_pressure_Pa=float(plateau_Pa), D_per_Pd=float(x),
                scope='closed-form inversion of specified alpha-Pd isotherm; not independent validation of its empirical vanthoff parameters')


def dilute_lag(parameters, thickness_m, temperature):
    diffusivity = parameters['D0']*np.exp(-parameters['Ea']/(KB_EV*temperature))
    lag = thickness_m**2/(6*diffusivity)
    n = np.arange(1, 10001, dtype=float)
    eigenmode_lag = 2*thickness_m**2/(np.pi**2*diffusivity)*np.sum((-1.)**(n+1)/n**2)
    if abs(eigenmode_lag/lag-1) > 1e-6:
        raise ValueError('Diffusion eigenmode moment disagrees with the lag formula')
    return dict(T_K=temperature, thickness_m=thickness_m,
                diffusivity_m2_s=float(diffusivity), lag_s=float(lag),
                independent_eigenmode_lag_s=float(eigenmode_lag),
                scope='dilute constant-D permeation time lag; not nonlinear hydride-front loading or a full equilibration time')


def film_lifetimes(model, energies_lab_keV=(.5, 1., 2.5, 6.5, 13.)):
    result = []
    initial = 2e-9
    terminal = .1e-9
    smoothing = .02e-9
    for name in ['TiO2', 'Al2O3']:
        material = model.OVERLAYERS[name]
        for energy in energies_lab_keV:
            current = float(np.interp(energy, model.IV_V, model.IV_I))*1e-6
            flux = current/elementary_charge*1.15/1e-4
            effective_yield = (.85*float(model.sputter_yield(energy*1000, name))
                               +.30*float(model.sputter_yield(energy*500, name)))/1.15
            velocity = flux*effective_yield/(material['n']*material['natoms'])
            simple = initial/velocity
            remove95 = (initial-terminal+smoothing*np.log(initial/terminal))/velocity
            row = dict(film=name, E_lab_keV=energy, thickness_m=initial,
                       incident_D_flux_m2_s=flux,
                       sputter_yield_atoms_per_incident_D=effective_yield,
                       nonregrowing_linear_removal_h=simple/3600,
                       nonregrowing_95percent_removal_h=remove95/3600)
            if name == 'TiO2':
                growth = 3.513e22*3e-7/np.sqrt(18*300)*1e4*.5/(2*material['n'])
                balance = lambda d: growth*np.exp(-d/.4e-9)-velocity*d/(d+smoothing)
                steady = brentq(balance, 0., 20e-9, xtol=1e-18)
                row.update(native_water_regrowth_assumptions=dict(T_K=300., pH2O_torr=3e-7, sticking=.5),
                           native_steady_thickness_nm=steady*1e9,
                           native_95percent_removal_possible=bool(steady < terminal))
            result.append(row)
    return dict(rows=result,
                scope='independent incident-particle/atom-inventory balance using specified sputter fits; film sputter-removal time is not a measured permeability lifetime; native regrowth may prevent removal')
