#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:01:39 2016

@author: u0092172, Sam Eyley

Python 3 script for the calculation of cellulose elemental analysis data
and empirical formulae. This script requires a large number of command line
arguments. Help can be obtained by -h or --help.

Refactored version with reduced code repetition.
"""

import click
import os
import sys
import lmfit as lmf
import numpy as np

# ============================================================================
# Constants and Configuration
# ============================================================================

# Atomic masses dictionary
ATOMIC_MASSES = {
    'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007,
    'S': 32.06, 'Si': 28.0855, 'F': 18.998, 'Cl': 35.45,
    'Br': 79.904, 'I': 126.904, 'Fe': 55.845, 'Cu': 63.546,
    'Na': 22.98976928
}

# Elements that can be present in modifications
MOD_ELEMENTS = ['C', 'H', 'O', 'N', 'S', 'Si', 'F', 'Cl', 'Br', 'I', 'Fe', 'Cu', 'Na']

# Cellulose composition
CELLULOSE = {'C': 6, 'H': 10, 'O': 5}

# Mapping of element symbols to click argument names
DETECTION_MAP = {
    'C': 'carb', 'H': 'hydr', 'O': 'oxyg', 'N': 'nitr', 'S': 'sulf',
    'Si': 'sili', 'F': 'fluo', 'Cl': 'chlor', 'Br': 'brom',
    'I': 'iodi', 'Fe': 'iron', 'Cu': 'copp', 'Na': 'sodium'
}

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

# ============================================================================
# Decorator Functions for Click Options
# ============================================================================

def add_modification_options(num_mods=3):
    """Decorator factory to add modification element options for all modifications"""
    def decorator(func):
        # Add options in reverse order (click requirement)
        for mod_num in range(num_mods, 0, -1):
            for elem in reversed(MOD_ELEMENTS):
                option_name = f'--{elem}{mod_num}'
                default = 1.0 if elem in ['H', 'O'] else 0.0
                required = (elem in ['H', 'O'] and mod_num == 1)
                help_text = f"Number of {elem} atoms in Mod {mod_num}"
                
                func = click.option(
                    option_name,
                    type=click.FLOAT,
                    default=float(default),
                    required=required,
                    help=help_text
                )(func)
        return func
    return decorator


def add_guess_options(num_mods=3):
    """Decorator factory to add DS guess and constraint options"""
    def decorator(func):
        for mod_num in range(num_mods, 0, -1):
            func = click.option(
                f'--G{mod_num}max',
                help=f"DS{mod_num} maximum",
                type=click.FLOAT,
                default=float(3)
            )(func)
            func = click.option(
                f'--G{mod_num}min',
                help=f"DS{mod_num} minimum",
                type=click.FLOAT,
                default=float(0)
            )(func)
            func = click.option(
                f'--G{mod_num}',
                help=f"Initial guess DS{mod_num}",
                type=click.FLOAT,
                default=float(0.5)
            )(func)
        return func
    return decorator


# ============================================================================
# Main CLI Command
# ============================================================================

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--humidity", help="Sample water content in percent", type=click.FLOAT, default=float(0))
@click.argument("carb", type=click.FLOAT)
@click.argument("hydr", type=click.FLOAT)
@click.argument("nitr", type=click.FLOAT)
@click.option("--sulf", help="Detected sulfur percent", type=click.FLOAT, default=float(0))
@click.option("--nmod", help="Number of modifications", required=True, type=click.INT, default=1)
@click.option("--oxyg", help="Detected oxygen percent", type=click.FLOAT, default=float(0))
@click.option("--sili", help="Detected silicon percent", type=click.FLOAT, default=float(0))
@click.option("--fluo", help="Detected fluorine percent", type=click.FLOAT, default=float(0))
@click.option("--chlor", help="Detected chlorine percent", type=click.FLOAT, default=float(0))
@click.option("--brom", help="Detected bromine percent", type=click.FLOAT, default=float(0))
@click.option("--iodi", help="Detected iodine percent", type=click.FLOAT, default=float(0))
@click.option("--iron", help="Detected iron percent", type=click.FLOAT, default=float(0))
@click.option("--copp", help="Detected copper percent", type=click.FLOAT, default=float(0))
@click.option("--sodium", help="Detected sodium percent", type=click.FLOAT, default=float(0))
@click.option("--varywater", help="Fit water content (Experimental, do NOT use!)", is_flag=True)
@click.option("--oxywater", help="Detected water percent during oxygen determination (if different for CHNS)", type=click.FLOAT, default=float(0))
@click.option("--chain_ratio", help="Ratio of surface to total chains (DSsurf). Set to 1 for DS or 6nm (0.375) 7.5nm (0.305) Elazzouzi-Hafraoui (0.237)", type=click.FLOAT, default=float(1))
@add_modification_options(num_mods=3)
@add_guess_options(num_mods=3)
@click.option("--watermax", help="water maximum (percent) only with vary water", type=click.FLOAT, default=float(1))
@click.option("--cweight", help="Carbon weighting during minimization", type=click.FLOAT, default=float(1e5))
@click.option("--hweight", help="Hydrogen weighting during minimization", type=click.FLOAT, default=float(1))
@click.option("--nweight", help="Nitrogen weighting during minimization", type=click.FLOAT, default=float(1e6))
@click.option("--sweight", help="Sulfur weighting during minimization", type=click.FLOAT, default=float(1e6))
@click.option("--oweight", help="Oxygen weighting during minimization", type=click.FLOAT, default=float(1e5))
@click.option("--method", help="Minimization routine (nelder or leastsq recommended) see lmfit documentation", default="nelder")
def cli(**kwargs):
    """This script will calculate the DS for modified cellulose samples based on (CHNS)
     elemental analysis data. Type help to see a list of required and optional arguments.

     CARB    Detected carbon concentration
     HYDR    Detected hydrogen concentration
     NITR    Detected nitrogen concentration"""

    # ========================================================================
    # Parse Input Data
    # ========================================================================
    
    # Convert detection percentages to fractions
    detected = {elem: kwargs[arg_name] / 100 for elem, arg_name in DETECTION_MAP.items()}
    
    # Parse modification compositions (3 modifications)
    modifications = []
    for mod_num in range(1, 4):
        mod = {elem: kwargs[f'{elem.lower()}{mod_num}'] for elem in MOD_ELEMENTS}
        modifications.append(mod)
    
    # Parse weights for residual calculation
    weights = {
        'C': kwargs['cweight'],
        'H': kwargs['hweight'],
        'O': kwargs['oweight'],
        'N': kwargs['nweight'],
        'S': kwargs['sweight'],
        'Si': 1, 'Cl': 1, 'Br': 1, 'F': 1, 'I': 1, 'Fe': 1, 'Cu': 1, 'Na': 1
    }
    
    # Extract other parameters
    humidity = kwargs['humidity']
    varywater = kwargs['varywater']
    oxywater = kwargs['oxywater']
    chain_ratio = kwargs['chain_ratio']
    nmod = kwargs['nmod']
    watermax = kwargs['watermax']
    method = kwargs['method']

    # ========================================================================
    # Define Residual Function
    # ========================================================================
    
    def residual(pars):
        """Calculate residuals between detected and calculated elemental fractions"""
        parvals = pars.valuesdict()
        DS = [parvals['DS1'], parvals['DS2'], parvals['DS3']]
        
        # Calculate modified composition for each element
        composition = {}
        for elem in MOD_ELEMENTS:
            # Start with cellulose composition (or 0 if not in cellulose)
            composition[elem] = CELLULOSE.get(elem, 0)
            
            # Subtract hydroxyl groups for H and O
            if elem in ['H', 'O']:
                composition[elem] -= sum(DS)
            
            # Add modification contributions
            for i, ds in enumerate(DS):
                composition[elem] += ds * modifications[i][elem]
        
        # Calculate total molecular masses
        total_mass = sum(composition[elem] * ATOMIC_MASSES[elem] for elem in MOD_ELEMENTS)
        
        # Handle water fractions
        water_frac = parvals['water']
        water_frac2 = oxywater / 100 if oxywater > 0 else water_frac
        
        water_H_mass = 2 * ATOMIC_MASSES['H']
        water_O_mass = 1 * ATOMIC_MASSES['O']
        water_tot_mass = water_H_mass + water_O_mass
        
        # Calculate mass fractions and residuals
        residuals = []
        for elem in MOD_ELEMENTS:
            elem_mass = composition[elem] * ATOMIC_MASSES[elem]
            
            # Calculate fraction based on element type
            if elem == 'H':
                frac = (elem_mass / total_mass) * (1 - water_frac) + \
                       water_frac * (water_H_mass / water_tot_mass)
            elif elem == 'O':
                frac = (elem_mass / total_mass) * (1 - water_frac2) + \
                       water_frac2 * (water_O_mass / water_tot_mass)
            else:
                frac = (elem_mass / total_mass) * (1 - water_frac)
            
            # Calculate residual if element was detected (not close to 0)
            if not np.allclose(detected[elem], 0):
                resid = ((detected[elem] - frac) ** 2) * weights[elem]
            else:
                resid = 0
            
            residuals.append(resid)
        
        return np.array(residuals)

    # ========================================================================
    # Setup and Run Optimization
    # ========================================================================
    
    # Create model from residual and parameterize
    model1 = lmf.Model(residual)
    
    # Set DS parameters based on number of modifications
    if nmod == 1:
        model1.set_param_hint('DS1', value=kwargs['g1'], min=kwargs['g1min'], max=kwargs['g1max'])
        model1.set_param_hint('DS2', value=0, vary=False)
        model1.set_param_hint('DS3', value=0, vary=False)
    elif nmod == 2:
        model1.set_param_hint('DS1', value=kwargs['g1'], min=kwargs['g1min'], max=kwargs['g1max'])
        model1.set_param_hint('DS2', value=kwargs['g2'], min=kwargs['g2min'], max=kwargs['g2max'])
        model1.set_param_hint('DS3', value=0, vary=False)
    else:
        model1.set_param_hint('DS1', value=kwargs['g1'], min=kwargs['g1min'], max=kwargs['g1max'])
        model1.set_param_hint('DS2', value=kwargs['g2'], min=kwargs['g2min'], max=kwargs['g2max'])
        model1.set_param_hint('DS3', value=kwargs['g3'], min=kwargs['g3min'], max=kwargs['g3max'])
    
    # Set water parameter
    if varywater:
        model1.set_param_hint('water', value=humidity/100, 
                            min=max(0, (humidity-0.5)/100),
                            max=max((humidity+2)/100, watermax/100))
    else:
        model1.set_param_hint('water', value=humidity/100, vary=False)
    
    pars = model1.make_params()
    
    # Run minimization
    result = lmf.minimize(residual, pars, method=method)
    
    # ========================================================================
    # Format and Output Results
    # ========================================================================
    
    report = lmf.fit_report(result)
    DS = [result.params['DS1'].value, result.params['DS2'].value, result.params['DS3'].value]
    water_result = result.params['water'].value
    
    # Calculate final composition
    final_composition = {}
    for elem in MOD_ELEMENTS:
        final_composition[elem] = CELLULOSE.get(elem, 0)
        if elem in ['H', 'O']:
            final_composition[elem] -= sum(DS)
        for i, ds in enumerate(DS):
            final_composition[elem] += ds * modifications[i][elem]
    
    # Calculate masses and fractions
    element_masses = {elem: final_composition[elem] * ATOMIC_MASSES[elem] 
                     for elem in MOD_ELEMENTS}
    total_mass = sum(element_masses.values())
    
    water_frac = water_result
    water_frac2 = oxywater / 100 if oxywater > 0 else water_frac
    water_H_mass = 2 * ATOMIC_MASSES['H']
    water_O_mass = 1 * ATOMIC_MASSES['O']
    water_tot_mass = water_H_mass + water_O_mass
    
    calculated_fractions = {}
    for elem in MOD_ELEMENTS:
        if elem == 'H':
            calculated_fractions[elem] = (element_masses[elem] / total_mass) * (1 - water_frac) + \
                                        water_frac * (water_H_mass / water_tot_mass)
        elif elem == 'O':
            calculated_fractions[elem] = (element_masses[elem] / total_mass) * (1 - water_frac2) + \
                                        water_frac2 * (water_O_mass / water_tot_mass)
        else:
            calculated_fractions[elem] = (element_masses[elem] / total_mass) * (1 - water_frac)
    
    # Calculate mass fraction of modification
    mod_masses = {elem: sum(DS[i] * modifications[i][elem] for i in range(3)) * ATOMIC_MASSES[elem]
                  for elem in MOD_ELEMENTS}
    total_mod_mass = sum(mod_masses.values())
    frac_mod = total_mod_mass / total_mass
    
    # Build empirical formula string
    formula_parts = []
    for elem in MOD_ELEMENTS:
        if final_composition[elem] > 0:
            formula_parts.append(f"{elem}{round(final_composition[elem], 3)}")
    empirical_formula = ' '.join(formula_parts)
    
    # Build modification description strings
    mod_descriptions = []
    for mod_num in range(1, 4):
        mod_parts = []
        for elem in MOD_ELEMENTS:
            val = modifications[mod_num-1][elem]
            if val != 0 or elem in ['H', 'O']:  # Always show H and O
                mod_parts.append(f"{elem}{round(val, 3)}")
        mod_descriptions.append(' '.join(mod_parts))
    
    # Build found/calculated element strings
    found_parts = []
    calc_parts = []
    for elem in MOD_ELEMENTS:
        found_parts.append(f"{elem}, {float(100*detected[elem])}")
        calc_parts.append(f"{elem}, {round(float(100*calculated_fractions[elem]), 3)}")
    
    # Build complete output string
    printout = (
        f'empirical formula: {empirical_formula}' + os.linesep +
        f'Chain ratio: {chain_ratio}' + os.linesep +
        f'Mod 1: {mod_descriptions[0]}' + os.linesep +
        f'Mod 2: {mod_descriptions[1]}' + os.linesep +
        f'Mod 3: {mod_descriptions[2]}' + os.linesep +
        f'DS1 = {round(DS[0], 3)}' + os.linesep +
        f'DS2 = {round(DS[1], 3)}' + os.linesep +
        f'DS3 = {round(DS[2], 3)}' + os.linesep +
        f'DS1(surf) = {round(DS[0]/chain_ratio, 3)}' + os.linesep +
        f'DS2(surf) = {round(DS[1]/chain_ratio, 3)}' + os.linesep +
        f'DS3(surf) = {round(DS[2]/chain_ratio, 3)}' + os.linesep +
        f'Found percent water = {humidity}' + os.linesep +
        f'Found percent water (oxygen) = {oxywater}' + os.linesep +
        f'Calc percent water = {water_result*100}' + os.linesep +
        f'Mass fraction modification = {round(frac_mod, 3)}' + os.linesep +
        f'Molecular mass modified AGU = {round(total_mass, 5)}' + os.linesep +
        'Found: ' + os.linesep +
        '; '.join(found_parts) + os.linesep +
        'Calculated: ' + os.linesep +
        '; '.join(calc_parts) + os.linesep + os.linesep +
        'Optimization data:' + os.linesep + os.linesep + str(report) +
        os.linesep + os.linesep + "Command line:" + os.linesep +
        os.linesep + str(sys.argv)
    )
    
    # Write to file and print
    with open('result.txt', 'w') as text_file:
        print(printout, file=text_file)
    print(printout)


if __name__ == '__main__':
    cli()
