#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:01:39 2016

@author: u0092172, Sam Eyley

Python 3 script for the calculation of cellulose elemental analysis data
and empirical formulae.  This script requires a large number of command line
arguments.  Help can be obtained by -h or --help.

Requires lmfit!  install with pip if using anaconda on windows, otherwise use
your favourite package manager.

For windows users:

python -m pip install lmfit

In order to get meaningful results, remember to constrain DS prior to fitting.

Example commandline:

./EA_calculator_lmfit.py 0.305 1.85 40.94 5.65 3.96 0.33 SE053 --O1 0 --H1 0
 --N1 3 --O2 4 --H2 1 --S2 1 --C3 7 --H3 7 --O3 3 --S3 1 --G1 0.5 --G1min 0.1
 --G1max 0.55 --G2 0.005 --G2min 0 --G2max 0.1 --G3 0.02 --G3min 0 --G3max 0.1
 
 Carbon, hydrogen, nitrogen and sulfur contents are mandatory, set to 0 if
 unknown or the script will not run!

"""

import argparse as ap
import os
import sys
import lmfit as lmf
import numpy as np

parser = ap.ArgumentParser(description='EA calculator for cellulose (CHNS minimum)')
parser.add_argument("chain_ratio", help=
                    "Chain ratio (DSsurf). Set to 1 for DS or 6nm (0.375) 7.5nm (0.305) Elazzouzi-Hafraoui (0.237)"
                    , type=np.float_)
parser.add_argument("humidity", help="Sample water content in percent", type=np.float_)
parser.add_argument("carb", help="Detected carbon percent",type=np.float_)
parser.add_argument("hydr", help="Detected hydrogen percent",type=np.float_)
parser.add_argument("nitr",help="Detected nitrogen percent",type=np.float_)
parser.add_argument("sulf",help="Detected sulfur percent",type=np.float_)
parser.add_argument("--oxyg",help="Detected oxygen percent",type=np.float_,default=np.float_(0))
parser.add_argument("--sili",help="Detected silicon percent",type=np.float_,default=np.float_(0))
parser.add_argument("--fluo",help="Detected fluorine percent",type=np.float_,default=np.float_(0))
parser.add_argument("--chlor",help="Detected chlorine percent",type=np.float_,default=np.float_(0))
parser.add_argument("--brom",help="Detected bromine percent",type=np.float_,default=np.float_(0))
parser.add_argument("--iodi",help="Detected iodine percent",type=np.float_,default=np.float_(0))
parser.add_argument("--iron",help="Detected iron percent",type=np.float_,default=np.float_(0))
parser.add_argument("--copp",help="Detected copper percent",type=np.float_,default=np.float_(0))
parser.add_argument("--sodium",help="Detected sodium percent",type=np.float_,default=np.float_(0))
parser.add_argument("--varywater",help="Fit water content (Experimental, do NOT use!)",action='store_true')
parser.add_argument("--oxywater",help="Detected water percent during oxygen determination (if different for CHNS)",type=np.float_,default=np.float_(0))


parser.add_argument("fname", help="Output filename (without extension)")
parser.add_argument("--nmod",help="Number of modifications",type=np.int_,default=1)

"""Define Modification"""
parser.add_argument("--C1",help="Number of carbons in Mod 1",type=np.float_,default=np.float_(0))
parser.add_argument("--H1",type=np.float_,default=np.float_(1))
parser.add_argument("--O1",type=np.float_,default=np.float_(1))
parser.add_argument("--N1",type=np.float_,default=np.float_(0))
parser.add_argument("--S1",type=np.float_,default=np.float_(0))
parser.add_argument("--Si1",type=np.float_,default=np.float_(0))
parser.add_argument("--F1",type=np.float_,default=np.float_(0))
parser.add_argument("--Cl1",type=np.float_,default=np.float_(0))
parser.add_argument("--Br1",type=np.float_,default=np.float_(0))
parser.add_argument("--I1",type=np.float_,default=np.float_(0))
parser.add_argument("--Fe1",type=np.float_,default=np.float_(0))
parser.add_argument("--Cu1",type=np.float_,default=np.float_(0))
parser.add_argument("--Na1",type=np.float_,default=np.float_(0))

parser.add_argument("--C2",help="Number of carbons in Mod 2",type=np.float_,default=np.float_(0))
parser.add_argument("--H2",type=np.float_,default=np.float_(1))
parser.add_argument("--O2",type=np.float_,default=np.float_(1))
parser.add_argument("--N2",type=np.float_,default=np.float_(0))
parser.add_argument("--S2",type=np.float_,default=np.float_(0))
parser.add_argument("--Si2",type=np.float_,default=np.float_(0))
parser.add_argument("--F2",type=np.float_,default=np.float_(0))
parser.add_argument("--Cl2",type=np.float_,default=np.float_(0))
parser.add_argument("--Br2",type=np.float_,default=np.float_(0))
parser.add_argument("--I2",type=np.float_,default=np.float_(0))
parser.add_argument("--Fe2",type=np.float_,default=np.float_(0))
parser.add_argument("--Cu2",type=np.float_,default=np.float_(0))
parser.add_argument("--Na2",type=np.float_,default=np.float_(0))

parser.add_argument("--C3",help="Number of carbons in Mod 3",type=np.float_,default=np.float_(0))
parser.add_argument("--H3",type=np.float_,default=np.float_(1))
parser.add_argument("--O3",type=np.float_,default=np.float_(1))
parser.add_argument("--N3",type=np.float_,default=np.float_(0))
parser.add_argument("--S3",type=np.float_,default=np.float_(0))
parser.add_argument("--Si3",type=np.float_,default=np.float_(0))
parser.add_argument("--F3",type=np.float_,default=np.float_(0))
parser.add_argument("--Cl3",type=np.float_,default=np.float_(0))
parser.add_argument("--Br3",type=np.float_,default=np.float_(0))
parser.add_argument("--I3",type=np.float_,default=np.float_(0))
parser.add_argument("--Fe3",type=np.float_,default=np.float_(0))
parser.add_argument("--Cu3",type=np.float_,default=np.float_(0))
parser.add_argument("--Na3",type=np.float_,default=np.float_(0))

"""Minimizer Constraints"""
parser.add_argument("--G1",help="Initial guess DS1",type=np.float_,default=np.float_(0.5))
parser.add_argument("--G1min",help="DS1 minimum",type=np.float_,default=np.float_(0))
parser.add_argument("--G1max",help="DS1 maximum",type=np.float_,default=np.float_(3))

parser.add_argument("--G2",help="Initial guess DS2",type=np.float_,default=np.float_(0.5))
parser.add_argument("--G2min",help="DS2 minimum",type=np.float_,default=np.float_(0))
parser.add_argument("--G2max",help="DS2 maximum",type=np.float_,default=np.float_(3))

parser.add_argument("--G3",help="Initial guess DS3",type=np.float_,default=np.float_(0.5))
parser.add_argument("--G3min",help="DS3 minimum",type=np.float_,default=np.float_(0))
parser.add_argument("--G3max",help="DS3 maximum",type=np.float_,default=np.float_(3))

parser.add_argument("--watermax",help="water maximum (percent) only with vary water",type=np.float_,default=np.float_(1))

parser.add_argument("--cweight",help="Carbon weighting during minimization",type=np.float_,default=np.float_(1e5))
parser.add_argument("--hweight",help="Hydrogen weighting during minimization",type=np.float_,default=np.float_(1))
parser.add_argument("--nweight",help="Nitrogen weighting during minimization",type=np.float_,default=np.float_(1e4))
parser.add_argument("--sweight",help="Sulfur weighting during minimization",type=np.float_,default=np.float_(1e4))
parser.add_argument("--oweight",help="Oxygen weighting during minimization",type=np.float_,default=np.float_(1e5))

parser.add_argument("--method", help="Minimization routine (nelder or leastsq recommended) see lmfit documentation",type=str,default="nelder")

args = parser.parse_args()

"""Detected Mass Concentrations"""
C_det = args.carb / 100
H_det = args.hydr / 100
O_det = args.oxyg / 100
N_det = args.nitr / 100
S_det = args.sulf / 100
Si_det = args.sili / 100
Cl_det = args.chlor / 100
Br_det = args.brom / 100
F_det = args.fluo / 100
I_det =  args.iodi / 100
Fe_det = args.iron / 100
Cu_det = args.copp / 100
Na_det = args.sodium / 100

"""Atomic masses"""
Cmass = np.float_(12.011)
Hmass = np.float_(1.008)
Omass = np.float_(15.999)
Nmass = np.float_(14.007)
Smass = np.float_(32.06)
Simass = np.float_(28.0855)
Fmass = np.float_(18.998)
Clmass = np.float_(35.45)
Brmass = np.float_(79.904)
Imass = np.float_(126.904)
Femass = np.float_(55.845)
Cumass = np.float_(63.546)
Namass = np.float_(22.98976928)

"""Define cellulose"""

c_cel = np.float_(6)
o_cel = np.float_(5)
h_cel = np.float_(10)


"""Define residual funtion to minimize"""
def residual (pars):
    parvals = pars.valuesdict()
    DS1 = parvals['DS1']
    DS2 = parvals['DS2']
    DS3 = parvals['DS3']

    H_mod = h_cel - DS1 - DS2 - DS3 + (DS1*args.H1) + (DS2*args.H2) + (DS3*args.H3)
    O_mod = o_cel - DS1 - DS2 - DS3 + (DS1*args.O1) + (DS2*args.O2) + (DS3*args.O3)
    
    C_mod = c_cel + (DS1*args.C1) + (DS2*args.C2) + (DS3*args.C3)
    N_mod = (DS1*args.N1) + (DS2*args.N2) + (DS3*args.N3)
    S_mod = (DS1*args.S1) + (DS2*args.S2) + (DS3*args.S3)
    Si_mod = (DS1*args.Si1) + (DS2*args.Si2) + (DS3*args.Si3)
    Cl_mod = (DS1*args.Cl1) + (DS2*args.Cl2) + (DS3*args.Cl3)
    Br_mod = (DS1*args.Br1) + (DS2*args.Br2) + (DS3*args.Br3)
    F_mod = (DS1*args.F1) + (DS2*args.F2) + (DS3*args.F3)
    I_mod = (DS1*args.I1) + (DS2*args.I2) + (DS3*args.I3)
    Fe_mod = (DS1*args.Fe1) + (DS2*args.Fe2) + (DS3*args.Fe3)
    Cu_mod = (DS1*args.Cu1) + (DS2*args.Cu2) + (DS3*args.Cu3)
    Na_mod =  (DS1*args.Na1) + (DS2*args.Na2) + (DS3*args.Na3)
    
    """Define empirical formula (include chain ratio here if required
    otherwise equal to X_mod"""

    C_emp = C_mod
    H_emp = H_mod
    O_emp = O_mod
    N_emp = N_mod
    S_emp = S_mod
    Si_emp = Si_mod
    Cl_emp = Cl_mod
    Br_emp = Br_mod
    F_emp = F_mod
    I_emp = I_mod
    Fe_emp = Fe_mod
    Cu_emp = Cu_mod
    Na_emp = Na_mod
    
    """Time to get massive"""
    
    C_tot = C_emp * Cmass
    H_tot = H_emp * Hmass
    O_tot = O_emp * Omass
    N_tot = N_emp * Nmass
    S_tot = S_emp * Smass
    Si_tot = Si_emp * Simass
    Cl_tot = Cl_emp * Clmass
    Br_tot = Br_emp * Brmass
    F_tot = F_emp * Fmass
    I_tot = I_emp * Imass
    Fe_tot = Fe_emp * Femass
    Cu_tot = Cu_emp * Cumass
    Na_tot = Na_emp * Namass
    
    M_tot = C_tot + H_tot + O_tot + N_tot + S_tot + Si_tot + Cl_tot + Br_tot + F_tot + I_tot + Fe_tot + Cu_tot + Na_tot
    
    water_frac = parvals['water']
    
    if args.oxywater > 0:
        water_frac2 = args.oxywater/100
        args.varywater = False
    else:
        water_frac2 = water_frac
    
    water_H = 2 * Hmass
    water_O = 1 * Omass
    water_tot = water_H + water_O

    C_frac = (C_tot / M_tot)*(1-water_frac)
    H_frac = (H_tot / M_tot)*(1-water_frac) + (water_frac*(water_H/water_tot))
    O_frac = (O_tot / M_tot)*(1-water_frac2) + (water_frac2*(water_O/water_tot))
    N_frac = (N_tot / M_tot)*(1-water_frac)
    S_frac = (S_tot / M_tot)*(1-water_frac)
    Si_frac = (Si_tot / M_tot)*(1-water_frac)
    Cl_frac = (Cl_tot / M_tot)*(1-water_frac)
    Br_frac = (Br_tot / M_tot)*(1-water_frac)
    F_frac = (F_tot / M_tot)*(1-water_frac)
    I_frac = (I_tot / M_tot)*(1-water_frac)
    Fe_frac = (Fe_tot / M_tot)*(1-water_frac)
    Cu_frac = (Cu_tot / M_tot)*(1-water_frac)
    Na_frac = (Na_tot / M_tot)*(1-water_frac)
    
    """Life gets complicated!"""
 
    if np.allclose(C_det,0) == False:
        residC = ((C_det-C_frac))**2
    else:
        residC = 0
    if np.allclose(H_det,0) == False:
        residH = ((H_det-H_frac))**2
    else:
        residH = 0
    if np.allclose(O_det,0) == False:
        residO = ((O_det-O_frac))**2
    else:
        residO = 0
    if np.allclose(N_det,0) == False:
        residN = ((N_det-N_frac))**2
    else:
        residN = 0
    if np.allclose(S_det,0) == False:
        residS = ((S_det-S_frac))**2
    else:
        residS = 0
    if np.allclose(Si_det,0) == False:
        residSi = ((Si_det-Si_frac))**2
    else:
        residSi = 0
    if np.allclose(Cl_det,0) == False:
        residCl = ((Cl_det-Cl_frac))**2
    else:
        residCl = 0
    if np.allclose(Br_det,0) == False:
        residBr = ((Br_det-Br_frac))**2
    else:
        residBr = 0
    if np.allclose(F_det,0) == False:
        residF = ((F_det-F_frac))**2
    else:
        residF = 0
    if np.allclose(I_det,0) == False:
        residI = ((I_det-I_frac))**2
    else: 
        residI = 0
    if np.allclose(Fe_det,0) == False:
        residFe = ((Fe_det-Fe_frac))**2
    else: 
        residFe = 0
    if np.allclose(Cu_det,0) == False:
        residCu = ((Cu_det-Cu_frac))**2
    else: 
        residCu = 0
    if np.allclose(Na_det,0) == False:
        residNa = ((Na_det-Na_frac))**2
    else: 
        residNa = 0     
    return np.array((1*residC*args.cweight,1*residH*args.hweight,1*residO*args.oweight,1*residN*args.nweight,
                     1*residS*args.sweight,1*residSi,1*residCl,1*residBr,1*residF,1*residI,1*residFe,
                     1*residCu,1*residNa))


"""Create model from residual and parameterize"""
model1 = lmf.Model(residual)
if args.nmod == 1:
    model1.set_param_hint('DS1',value=args.G1,min=args.G1min,max=args.G1max)
    model1.set_param_hint('DS2',value=0,vary=False)
    model1.set_param_hint('DS3',value=0,vary=False)
elif args.nmod == 2:
    model1.set_param_hint('DS1',value=args.G1,min=args.G1min,max=args.G1max)
    model1.set_param_hint('DS2',value=args.G2,min=args.G2min,max=args.G2max)
    model1.set_param_hint('DS3',value=0,vary=False)
else:
    model1.set_param_hint('DS1',value=args.G1,min=args.G1min,max=args.G1max)
    model1.set_param_hint('DS2',value=args.G2,min=args.G2min,max=args.G2max)
    model1.set_param_hint('DS3',value=args.G3,min=args.G3min,max=args.G3max)
if args.varywater == True:
    model1.set_param_hint('water',value=args.humidity/100,min=max(0,(args.humidity-0.5)/100),
                          max=max((args.humidity+2)/100,args.watermax/100))
else:
    model1.set_param_hint('water',value=args.humidity/100,vary=False)
pars = model1.make_params()

"""Here's the magic line. leastsq is LM algorithm, others are available."""
result = lmf.minimize(residual,pars,method=args.method)

"""Format results"""
report = lmf.fit_report(result)
res_dict = result.params.valuesdict()
DS = list(res_dict.values())


H_mod = h_cel - DS[0] - DS[1] - DS[2] + (DS[0]*args.H1) + (DS[1]*args.H2) + (DS[2]*args.H3)
O_mod = o_cel - DS[0] - DS[1] - DS[2] + (DS[0]*args.O1) + (DS[1]*args.O2) + (DS[2]*args.O3)

C_mod = c_cel + (DS[0]*args.C1) + (DS[1]*args.C2) + (DS[2]*args.C3)
N_mod = (DS[0]*args.N1) + (DS[1]*args.N2) + (DS[2]*args.N3)
S_mod = (DS[0]*args.S1) + (DS[1]*args.S2) + (DS[2]*args.S3)
Si_mod = (DS[0]*args.Si1) + (DS[1]*args.Si2) + (DS[2]*args.Si3)
Cl_mod = (DS[0]*args.Cl1) + (DS[1]*args.Cl2) + (DS[2]*args.Cl3)
Br_mod = (DS[0]*args.Br1) + (DS[1]*args.Br2) + (DS[2]*args.Br3)
F_mod = (DS[0]*args.F1) + (DS[1]*args.F2) + (DS[2]*args.F3)
I_mod = (DS[0]*args.I1) + (DS[1]*args.I2) + (DS[2]*args.I3)

Fe_mod = (DS[0]*args.Fe1) + (DS[1]*args.Fe2) + (DS[2]*args.Fe3)
Cu_mod = (DS[0]*args.Cu1) + (DS[1]*args.Cu2) + (DS[2]*args.Cu3)
Na_mod = (DS[0]*args.Na1) + (DS[1]*args.Na2) + (DS[2]*args.Na3)


"""Define empirical formula (include chain ratio here if required
otherwise equal to X_mod"""
C_emp = C_mod
H_emp = H_mod
O_emp = O_mod
N_emp = N_mod
S_emp = S_mod
Si_emp = Si_mod
Cl_emp = Cl_mod
Br_emp = Br_mod
F_emp = F_mod
I_emp = I_mod
Fe_emp = Fe_mod
Cu_emp = Cu_mod
Na_emp = Na_mod    

"""Time to get massive"""

C_tot = C_emp * Cmass
H_tot = H_emp * Hmass
O_tot = O_emp * Omass
N_tot = N_emp * Nmass
S_tot = S_emp * Smass
Si_tot = Si_emp * Simass
Cl_tot = Cl_emp * Clmass
Br_tot = Br_emp * Brmass
F_tot = F_emp * Fmass
I_tot = I_emp * Imass
Fe_tot = Fe_emp * Femass
Cu_tot = Cu_emp * Cumass
Na_tot = Na_emp * Namass

M_tot = C_tot + H_tot + O_tot + N_tot + S_tot + Si_tot + Cl_tot + Br_tot + F_tot + I_tot + Fe_tot + Cu_tot + Na_tot

water_frac = DS[3] #This takes the result from the fit

if args.oxywater >0:
    water_frac2 = args.oxywater/100
else:
    water_frac2 = water_frac

water_H = 2 * Hmass
water_O = 1 * Omass
water_tot = water_H + water_O

C_frac = (C_tot / M_tot)*(1-water_frac)
H_frac = (H_tot / M_tot)*(1-water_frac) + (water_frac*(water_H/water_tot))
O_frac = (O_tot / M_tot)*(1-water_frac2) + (water_frac2*(water_O/water_tot))
N_frac = (N_tot / M_tot)*(1-water_frac)
S_frac = (S_tot / M_tot)*(1-water_frac)
Si_frac = (Si_tot / M_tot)*(1-water_frac)
Cl_frac = (Cl_tot / M_tot)*(1-water_frac)
Br_frac = (Br_tot / M_tot)*(1-water_frac)
F_frac = (F_tot / M_tot)*(1-water_frac)
I_frac = (I_tot / M_tot)*(1-water_frac)
Fe_frac = (Fe_tot / M_tot)*(1-water_frac)
Cu_frac = (Cu_tot / M_tot)*(1-water_frac)
Na_frac = (Na_tot / M_tot)*(1-water_frac) 

"""Extra for calc of mass percent mod"""

H_mod1 = (DS[0]*args.H1) + (DS[1]*args.H2) + (DS[2]*args.H3)
O_mod1 = (DS[0]*args.O1) + (DS[1]*args.O2) + (DS[2]*args.O3)
C_mod1 = (DS[0]*args.C1) + (DS[1]*args.C2) + (DS[2]*args.C3)

C_tot1 = C_mod1 * Cmass
H_tot1 = H_mod1 * Hmass
O_tot1 = O_mod1 * Omass

M_tot1 = C_tot1 + H_tot1 + O_tot1 + N_tot + S_tot + Si_tot + Cl_tot + Br_tot + F_tot + I_tot + Fe_tot + Cu_tot + Na_tot

Frac_mod = M_tot1 / M_tot

with open(args.fname+'_result.txt','w') as text_file:
    print('Empirical formula: C'+str(round(C_emp,3))+" H"+str(round(H_emp,3))+' O'+str(round(O_emp,3))+' N'+str(round(N_emp,3))
          +' S'+str(round(S_emp,3))+' F'+str(round(F_emp,3))+' Cl'+str(round(Cl_emp,3))+' Br'+str(round(Br_emp,3))+' I'+str(round(I_emp,3))
          +' Fe'+str(round(Fe_emp,3))+' Cu'+str(round(Cu_emp,3))+' Na'+str(round(Na_emp,3))+os.linesep
          +'Chain ratio: '+str(args.chain_ratio)+os.linesep
          +'Mod 1: C'+str(round(args.C1,3))+' H'+str(round(args.H1,3))+' O'+str(round(args.O1,3))+' N'+str(round(args.N1,3))+' S'+str(round(args.S1,3))+' Si'+str(round(args.Si1,3))+' F'+str(round(args.F1,3))+' Cl'+str(round(args.Cl1,3))+' Br'+str(round(args.Br1,3))+' I'+str(round(args.I1,3))+' Fe'+str(round(args.Fe1,3))+' Cu'+str(round(args.Cu1,3))+' Na'+str(round(args.Na1,3))+os.linesep
          +'Mod 2: C'+str(round(args.C2,3))+' H'+str(round(args.H2,3))+' O'+str(round(args.O2,3))+' N'+str(round(args.N2,3))+' S'+str(round(args.S2,3))+' Si'+str(round(args.Si2,3))+' F'+str(round(args.F2,3))+' Cl'+str(round(args.Cl2,3))+' Br'+str(round(args.Br2,3))+' I'+str(round(args.I2,3))+' Fe'+str(round(args.Fe2,3))+' Cu'+str(round(args.Cu2,3))+' Na'+str(round(args.Na2,3))+os.linesep
          +'Mod 3: C'+str(round(args.C3,3))+' H'+str(round(args.H3,3))+' O'+str(round(args.O3,3))+' N'+str(round(args.N3,3))+' S'+str(round(args.S3,3))+' Si'+str(round(args.Si3,3))+' F'+str(round(args.F3,3))+' Cl'+str(round(args.Cl3,3))+' Br'+str(round(args.Br3,3))+' I'+str(round(args.I3,3))+' Fe'+str(round(args.Fe3,3))+' Cu'+str(round(args.Cu3,3))+' Na'+str(round(args.Na3,3))+os.linesep
          +'DS1 = '+str(round(DS[0],3))+os.linesep
          +'DS2 = '+str(round(DS[1],3))+os.linesep
          +'DS3 = '+str(round(DS[2],3))+os.linesep
          +'DS1(surf) = '+str(round(DS[0]/args.chain_ratio,3))+os.linesep
          +'DS2(surf) = '+str(round(DS[1]/args.chain_ratio,3))+os.linesep
          +'DS3(surf) = '+str(round(DS[2]/args.chain_ratio,3))+os.linesep
          +'Found percent water = '+str(args.humidity)+os.linesep
          +'Found percent water (oxygen) = '+str(args.oxywater)+os.linesep
          +'Calc percent water = '+str(DS[3]*100)+os.linesep
          +'Mass fraction modification = '+str(round(Frac_mod,3))+os.linesep
          +'Molecular mass modified AGU = '+str(round(M_tot,5))+os.linesep
          +'Found: '+os.linesep
          +'C, '+str(np.float_(100*C_det))+'; H, '+str(np.float_(100*H_det))+'; O, '+str(np.float_(100*O_det))+'; N, '+str(np.float_(100*N_det))
          +'; S, '+str(np.float_(100*S_det))+'; Si, '+str(np.float_(100*Si_det))+'; F, '+str(np.float_(100*F_det))+'; Cl, '+str(np.float_(100*Cl_det))
          +'; Br, '+str(np.float_(100*Br_det))+'; I, '+str(np.float_(100*I_det))+'; Fe, '+str(np.float_(100*Fe_det))+'; Cu, '+str(np.float_(100*Cu_det))+'; Na, '+str(np.float_(100*Na_det))+os.linesep
          +'Calculated: '+os.linesep
          +'C, '+str(round(np.float_(100*C_frac),3))+'; H, '+str(round(np.float_(100*H_frac),3))+'; O, '+str(round(np.float_(100*O_frac),3))+'; N, '+str(round(np.float_(100*N_frac),3))
          +'; S, '+str(round(np.float_(100*S_frac),3))+'; Si, '+str(round(np.float_(100*Si_frac),3))+'; F, '+str(round(np.float_(100*F_frac),3))+'; Cl, '+str(round(np.float_(100*Cl_frac),3))
          +'; Br, '+str(round(np.float_(100*Br_frac),3))+'; I, '+str(round(np.float_(100*I_frac),3))+'; Fe, '+str(round(np.float_(100*Fe_frac),3))+'; Cu, '+str(round(np.float_(100*Cu_frac),3))+'; Na, '+str(round(np.float_(100*Na_frac),3))
          +os.linesep+os.linesep+'Optimization data:'+os.linesep+os.linesep+str(report)
          +os.linesep+os.linesep+"Command line:"+os.linesep
          +os.linesep+str(sys.argv),file=text_file)
print('Empirical formula: C'+str(round(C_emp,3))+" H"+str(round(H_emp,3))+' O'+str(round(O_emp,3))+' N'+str(round(N_emp,3))
      +' S'+str(round(S_emp,3))+' F'+str(round(F_emp,3))+' Cl'+str(round(Cl_emp,3))+' Br'+str(round(Br_emp,3))+' I'+str(round(I_emp,3))
      +' Fe'+str(round(Fe_emp,3))+' Cu'+str(round(Cu_emp,3))+' Na'+str(round(Na_emp,3))+os.linesep
      +'Chain ratio: '+str(args.chain_ratio)+os.linesep
      +'Mod 1: C'+str(round(args.C1,3))+' H'+str(round(args.H1,3))+' O'+str(round(args.O1,3))+' N'+str(round(args.N1,3))+' S'+str(round(args.S1,3))+' Si'+str(round(args.Si1,3))+' F'+str(round(args.F1,3))+' Cl'+str(round(args.Cl1,3))+' Br'+str(round(args.Br1,3))+' I'+str(round(args.I1,3))+' Fe'+str(round(args.Fe1,3))+' Cu'+str(round(args.Cu1,3))+' Na'+str(round(args.Na1,3))+os.linesep
      +'Mod 2: C'+str(round(args.C2,3))+' H'+str(round(args.H2,3))+' O'+str(round(args.O2,3))+' N'+str(round(args.N2,3))+' S'+str(round(args.S2,3))+' Si'+str(round(args.Si2,3))+' F'+str(round(args.F2,3))+' Cl'+str(round(args.Cl2,3))+' Br'+str(round(args.Br2,3))+' I'+str(round(args.I2,3))+' Fe'+str(round(args.Fe2,3))+' Cu'+str(round(args.Cu2,3))+' Na'+str(round(args.Na2,3))+os.linesep
      +'Mod 3: C'+str(round(args.C3,3))+' H'+str(round(args.H3,3))+' O'+str(round(args.O3,3))+' N'+str(round(args.N3,3))+' S'+str(round(args.S3,3))+' Si'+str(round(args.Si3,3))+' F'+str(round(args.F3,3))+' Cl'+str(round(args.Cl3,3))+' Br'+str(round(args.Br3,3))+' I'+str(round(args.I3,3))+' Fe'+str(round(args.Fe3,3))+' Cu'+str(round(args.Cu3,3))+' Na'+str(round(args.Na3,3))+os.linesep
      +'DS1 = '+str(round(DS[0],3))+os.linesep
      +'DS2 = '+str(round(DS[1],3))+os.linesep
      +'DS3 = '+str(round(DS[2],3))+os.linesep
      +'DS1(surf) = '+str(round(DS[0]/args.chain_ratio,3))+os.linesep
      +'DS2(surf) = '+str(round(DS[1]/args.chain_ratio,3))+os.linesep
      +'DS3(surf) = '+str(round(DS[2]/args.chain_ratio,3))+os.linesep
      +'Found percent water = '+str(args.humidity)+os.linesep
      +'Found percent water (oxygen) = '+str(args.oxywater)+os.linesep
      +'Calc percent water = '+str(DS[3]*100)+os.linesep
      +'Mass fraction modification = '+str(round(Frac_mod,3))+os.linesep
      +'Molecular mass modified AGU = '+str(round(M_tot,5))+os.linesep
      +'Found: '+os.linesep
      +'C, '+str(np.float_(100*C_det))+'; H, '+str(np.float_(100*H_det))+'; O, '+str(np.float_(100*O_det))+'; N, '+str(np.float_(100*N_det))
      +'; S, '+str(np.float_(100*S_det))+'; Si, '+str(np.float_(100*Si_det))+'; F, '+str(np.float_(100*F_det))+'; Cl, '+str(np.float_(100*Cl_det))
      +'; Br, '+str(np.float_(100*Br_det))+'; I, '+str(np.float_(100*I_det))+'; Fe, '+str(np.float_(100*Fe_det))+'; Cu, '+str(np.float_(100*Cu_det))+'; Na, '+str(np.float_(100*Na_det))+os.linesep
      +'Calculated: '+os.linesep
      +'C, '+str(round(np.float_(100*C_frac),3))+'; H, '+str(round(np.float_(100*H_frac),3))+'; O, '+str(round(np.float_(100*O_frac),3))+'; N, '+str(round(np.float_(100*N_frac),3))
      +'; S, '+str(round(np.float_(100*S_frac),3))+'; Si, '+str(round(np.float_(100*Si_frac),3))+'; F, '+str(round(np.float_(100*F_frac),3))+'; Cl, '+str(round(np.float_(100*Cl_frac),3))
      +'; Br, '+str(round(np.float_(100*Br_frac),3))+'; I, '+str(round(np.float_(100*I_frac),3))+'; Fe, '+str(round(np.float_(100*Fe_frac),3))+'; Cu, '+str(round(np.float_(100*Cu_frac),3))+'; Na, '+str(round(np.float_(100*Na_frac),3))
      +os.linesep+os.linesep+'Optimization data:'+os.linesep+os.linesep+str(report)
      +os.linesep+os.linesep+"Command line:"+os.linesep
      +os.linesep+str(sys.argv))