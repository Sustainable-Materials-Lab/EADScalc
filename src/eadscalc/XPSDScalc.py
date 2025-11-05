#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:01:39 2016

@author: u0092172, Sam Eyley

Python 3 script for the calculation of cellulose elemental analysis data
and empirical formulae.  This script requires a large number of command line
arguments.  Help can be obtained by -h or --help.

In order to get meaningful results, remember to constrain DS prior to fitting.

Example commandline:

./EA_calculator_lmfit.py 0.305 1.85 40.94 5.65 3.96 0.33 SE053 --O1 0 --H1 0
 --N1 3 --O2 4 --H2 1 --S2 1 --C3 7 --H3 7 --O3 3 --S3 1 --G1 0.5 --G1min 0.1
 --G1max 0.55 --G2 0.005 --G2min 0 --G2max 0.1 --G3 0.02 --G3min 0 --G3max 0.1
 
 Carbon, hydrogen, nitrogen and sulfur contents are mandatory, set to 0 if
 unknown or the script will not run.

"""

import argparse as ap
import os
import sys
import lmfit as lmf
import numpy as np

parser = ap.ArgumentParser(description='DS calculator for cellulose from XPS')
parser.add_argument("carb", help="Detected carbon percent",type=float)
parser.add_argument("nitr",help="Detected nitrogen percent",type=float)
parser.add_argument("sulf",help="Detected sulfur percent",type=float)

parser.add_argument("--fname", help="Output filename (without extension)",default="result")
parser.add_argument("--nmod",help="Number of modifications",required=True,type=int,default=1)

parser.add_argument("--oxyg",help="Detected oxygen percent",type=float,default=float(0))
parser.add_argument("--sili",help="Detected silicon percent",type=float,default=float(0))
parser.add_argument("--fluo",help="Detected fluorine percent",type=float,default=float(0))
parser.add_argument("--chlor",help="Detected chlorine percent",type=float,default=float(0))
parser.add_argument("--brom",help="Detected bromine percent",type=float,default=float(0))
parser.add_argument("--iodi",help="Detected iodine percent",type=float,default=float(0))
parser.add_argument("--iron",help="Detected iron percent",type=float,default=float(0))
parser.add_argument("--copp",help="Detected copper percent",type=float,default=float(0))
parser.add_argument("--sodium",help="Detected sodium percent",type=float,default=float(0))


"""Define Modification"""
parser.add_argument("--C1",help="Number of carbons in Mod 1",type=float,default=float(0))
parser.add_argument("--H1",type=float,default=float(1),required=True)
parser.add_argument("--O1",type=float,default=float(1),required=True)
parser.add_argument("--N1",type=float,default=float(0))
parser.add_argument("--S1",type=float,default=float(0))
parser.add_argument("--Si1",type=float,default=float(0))
parser.add_argument("--F1",type=float,default=float(0))
parser.add_argument("--Cl1",type=float,default=float(0))
parser.add_argument("--Br1",type=float,default=float(0))
parser.add_argument("--I1",type=float,default=float(0))
parser.add_argument("--Fe1",type=float,default=float(0))
parser.add_argument("--Cu1",type=float,default=float(0))
parser.add_argument("--Na1",type=float,default=float(0))

parser.add_argument("--C2",help="Number of carbons in Mod 2",type=float,default=float(0))
parser.add_argument("--H2",type=float,default=float(1))
parser.add_argument("--O2",type=float,default=float(1))
parser.add_argument("--N2",type=float,default=float(0))
parser.add_argument("--S2",type=float,default=float(0))
parser.add_argument("--Si2",type=float,default=float(0))
parser.add_argument("--F2",type=float,default=float(0))
parser.add_argument("--Cl2",type=float,default=float(0))
parser.add_argument("--Br2",type=float,default=float(0))
parser.add_argument("--I2",type=float,default=float(0))
parser.add_argument("--Fe2",type=float,default=float(0))
parser.add_argument("--Cu2",type=float,default=float(0))
parser.add_argument("--Na2",type=float,default=float(0))

parser.add_argument("--C3",help="Number of carbons in Mod 3",type=float,default=float(0))
parser.add_argument("--H3",type=float,default=float(1))
parser.add_argument("--O3",type=float,default=float(1))
parser.add_argument("--N3",type=float,default=float(0))
parser.add_argument("--S3",type=float,default=float(0))
parser.add_argument("--Si3",type=float,default=float(0))
parser.add_argument("--F3",type=float,default=float(0))
parser.add_argument("--Cl3",type=float,default=float(0))
parser.add_argument("--Br3",type=float,default=float(0))
parser.add_argument("--I3",type=float,default=float(0))
parser.add_argument("--Fe3",type=float,default=float(0))
parser.add_argument("--Cu3",type=float,default=float(0))
parser.add_argument("--Na3",type=float,default=float(0))

"""Minimizer Constraints"""
parser.add_argument("--G1",help="Initial guess DS1",type=float,default=float(0.5))
parser.add_argument("--G1min",help="DS1 minimum",type=float,default=float(0))
parser.add_argument("--G1max",help="DS1 maximum",type=float,default=float(3))

parser.add_argument("--G2",help="Initial guess DS2",type=float,default=float(0.5))
parser.add_argument("--G2min",help="DS2 minimum",type=float,default=float(0))
parser.add_argument("--G2max",help="DS2 maximum",type=float,default=float(3))

parser.add_argument("--G3",help="Initial guess DS3",type=float,default=float(0.5))
parser.add_argument("--G3min",help="DS3 minimum",type=float,default=float(0))
parser.add_argument("--G3max",help="DS3 maximum",type=float,default=float(3))

parser.add_argument("--cweight",help="Carbon weighting during minimization",type=float,default=float(1))
parser.add_argument("--nweight",help="Nitrogen weighting during minimization",type=float,default=float(1))
parser.add_argument("--sweight",help="Sulfur weighting during minimization",type=float,default=float(1))
parser.add_argument("--oweight",help="Oxygen weighting during minimization",type=float,default=float(1))

parser.add_argument("--method", help="Minimization routine (nelder or leastsq recommended) see lmfit documentation",type=str,default="nelder")

args = parser.parse_args()

"""Detected At. Concentrations"""
C_det = args.carb / 100
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

"""Define cellulose"""

c_cel = float(6)
o_cel = float(5)
h_cel = float(10)


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
    
    
    
    """Time to get massive"""
        
    M_mod = C_mod + O_mod + N_mod + S_mod + Si_mod + Cl_mod + Br_mod + F_mod + I_mod + Fe_mod + Cu_mod + Na_mod
    

    C_frac = (C_mod / M_mod)
    O_frac = (O_mod / M_mod)
    N_frac = (N_mod / M_mod)
    S_frac = (S_mod / M_mod)
    Si_frac = (Si_mod / M_mod)
    Cl_frac = (Cl_mod / M_mod)
    Br_frac = (Br_mod / M_mod)
    F_frac = (F_mod / M_mod)
    I_frac = (I_mod / M_mod)
    Fe_frac = (Fe_mod / M_mod)
    Cu_frac = (Cu_mod / M_mod)
    Na_frac = (Na_mod / M_mod)
    
    """Life gets complicated!"""
 
    if np.allclose(C_det,0) == False:
        residC = ((C_det-C_frac))**2
    else:
        residC = 0
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
    return np.array((1*residC*args.cweight,1*residO*args.oweight,1e7*residN*args.nweight,
                     1*residS*args.sweight,1*residSi,1e5*residCl,1e5*residBr,1*residF,1*residI,1*residFe,
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
pars = model1.make_params()

"""Here's the magic line. leastsq is LM algorithm, others are available."""
result = lmf.minimize(residual,pars,method=args.method)

"""Format results"""
report = lmf.fit_report(result)
DS = list(result.params.valuesdict().values())

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



M_mod = C_mod + O_mod + N_mod + S_mod + Si_mod + Cl_mod + Br_mod + F_mod + I_mod + Fe_mod + Cu_mod + Na_mod


C_frac = (C_mod / M_mod)
O_frac = (O_mod / M_mod)
N_frac = (N_mod / M_mod)
S_frac = (S_mod / M_mod)
Si_frac = (Si_mod / M_mod)
Cl_frac = (Cl_mod / M_mod)
Br_frac = (Br_mod / M_mod)
F_frac = (F_mod / M_mod)
I_frac = (I_mod / M_mod)
Fe_frac = (Fe_mod / M_mod)
Cu_frac = (Cu_mod / M_mod)
Na_frac = (Na_mod / M_mod)

"""Extra for calc of mass percent mod"""

printout = ('empirical formula: C'+str(round(C_mod,3))+' O'+str(round(O_mod,3))+' N'+str(round(N_mod,3))
      +' S'+str(round(S_mod,3))+' F'+str(round(F_mod,3))+' Cl'+str(round(Cl_mod,3))+' Br'+str(round(Br_mod,3))+' I'+str(round(I_mod,3))
      +' Fe'+str(round(Fe_mod,3))+' Cu'+str(round(Cu_mod,3))+' Na'+str(round(Na_mod,3))+os.linesep
      +'Mod 1: C'+str(round(args.C1,3))+' H'+str(round(args.H1,3))+' O'+str(round(args.O1,3))+' N'+str(round(args.N1,3))+' S'+str(round(args.S1,3))+' Si'+str(round(args.Si1,3))+' F'+str(round(args.F1,3))+' Cl'+str(round(args.Cl1,3))+' Br'+str(round(args.Br1,3))+' I'+str(round(args.I1,3))+' Fe'+str(round(args.Fe1,3))+' Cu'+str(round(args.Cu1,3))+' Na'+str(round(args.Na1,3))+os.linesep
      +'Mod 2: C'+str(round(args.C2,3))+' H'+str(round(args.H2,3))+' O'+str(round(args.O2,3))+' N'+str(round(args.N2,3))+' S'+str(round(args.S2,3))+' Si'+str(round(args.Si2,3))+' F'+str(round(args.F2,3))+' Cl'+str(round(args.Cl2,3))+' Br'+str(round(args.Br2,3))+' I'+str(round(args.I2,3))+' Fe'+str(round(args.Fe2,3))+' Cu'+str(round(args.Cu2,3))+' Na'+str(round(args.Na2,3))+os.linesep
      +'Mod 3: C'+str(round(args.C3,3))+' H'+str(round(args.H3,3))+' O'+str(round(args.O3,3))+' N'+str(round(args.N3,3))+' S'+str(round(args.S3,3))+' Si'+str(round(args.Si3,3))+' F'+str(round(args.F3,3))+' Cl'+str(round(args.Cl3,3))+' Br'+str(round(args.Br3,3))+' I'+str(round(args.I3,3))+' Fe'+str(round(args.Fe3,3))+' Cu'+str(round(args.Cu3,3))+' Na'+str(round(args.Na3,3))+os.linesep
      +'DS1 = '+str(round(DS[0],3))+os.linesep
      +'DS2 = '+str(round(DS[1],3))+os.linesep
      +'DS3 = '+str(round(DS[2],3))+os.linesep
      +'Molecular mass modified AGU = '+str(round(M_mod,5))+os.linesep
      +'Found: '+os.linesep
      +'C, '+str(float(100*C_det))+'; O, '+str(float(100*O_det))+'; N, '+str(float(100*N_det))
      +'; S, '+str(float(100*S_det))+'; Si, '+str(float(100*Si_det))+'; F, '+str(float(100*F_det))+'; Cl, '+str(float(100*Cl_det))
      +'; Br, '+str(float(100*Br_det))+'; I, '+str(float(100*I_det))+'; Fe, '+str(float(100*Fe_det))+'; Cu, '+str(float(100*Cu_det))+'; Na, '+str(float(100*Na_det))+os.linesep
      +'Calculated: '+os.linesep
      +'C, '+str(round(float(100*C_frac),3))+'; O, '+str(round(float(100*O_frac),3))+'; N, '+str(round(float(100*N_frac),3))
      +'; S, '+str(round(float(100*S_frac),3))+'; Si, '+str(round(float(100*Si_frac),3))+'; F, '+str(round(float(100*F_frac),3))+'; Cl, '+str(round(float(100*Cl_frac),3))
      +'; Br, '+str(round(float(100*Br_frac),3))+'; I, '+str(round(float(100*I_frac),3))+'; Fe, '+str(round(float(100*Fe_frac),3))+'; Cu, '+str(round(float(100*Cu_frac),3))+'; Na, '+str(round(float(100*Na_frac),3))
      +os.linesep+os.linesep+'Optimization data:'+os.linesep+os.linesep+str(report)
      +os.linesep+os.linesep+"Command line:"+os.linesep
      +os.linesep+str(sys.argv))

with open(args.fname+'_result.txt','w') as text_file:
    print(printout,file=text_file)
print(printout)
