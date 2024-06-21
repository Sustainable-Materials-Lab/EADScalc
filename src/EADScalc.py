#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:01:39 2016

@author: u0092172, Sam Eyley

Python 3 script for the calculation of cellulose elemental analysis data
and empirical formulae.  This script requires a large number of command line
arguments.  Help can be obtained by -h or --help.

"""

import click
import os
import sys
import lmfit as lmf
import numpy as np

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--humidity", help="Sample water content in percent",type=click.FLOAT,default=float(0))
@click.argument("carb",type=click.FLOAT)
@click.argument("hydr",type=click.FLOAT)
@click.argument("nitr",type=click.FLOAT)
@click.option("--sulf", help="Detected sulfur percent",type=click.FLOAT,default=float(0))

@click.option("--nmod",help="Number of modifications",required=True,type=click.INT,default=1)

@click.option("--oxyg",help="Detected oxygen percent",type=click.FLOAT,default=float(0))
@click.option("--sili",help="Detected silicon percent",type=click.FLOAT,default=float(0))
@click.option("--fluo",help="Detected fluorine percent",type=click.FLOAT,default=float(0))
@click.option("--chlor",help="Detected chlorine percent",type=click.FLOAT,default=float(0))
@click.option("--brom",help="Detected bromine percent",type=click.FLOAT,default=float(0))
@click.option("--iodi",help="Detected iodine percent",type=click.FLOAT,default=float(0))
@click.option("--iron",help="Detected iron percent",type=click.FLOAT,default=float(0))
@click.option("--copp",help="Detected copper percent",type=click.FLOAT,default=float(0))
@click.option("--sodium",help="Detected sodium percent",type=click.FLOAT,default=float(0))
@click.option("--varywater",help="Fit water content (Experimental, do NOT use!)",is_flag=True)
@click.option("--oxywater",help="Detected water percent during oxygen determination (if different for CHNS)",type=click.FLOAT,default=float(0))
@click.option("--chain_ratio", help="Ratio of surface to total chains (DSsurf). Set to 1 for DS or 6nm (0.375) 7.5nm (0.305) Elazzouzi-Hafraoui (0.237)",
              type=click.FLOAT,default=float(1))
@click.option("--C1",help="Number of carbons in Mod 1",type=click.FLOAT,default=float(0))
@click.option("--H1",type=click.FLOAT,default=float(1),required=True)
@click.option("--O1",type=click.FLOAT,default=float(1),required=True)
@click.option("--N1",type=click.FLOAT,default=float(0))
@click.option("--S1",type=click.FLOAT,default=float(0))
@click.option("--Si1",type=click.FLOAT,default=float(0))
@click.option("--F1",type=click.FLOAT,default=float(0))
@click.option("--Cl1",type=click.FLOAT,default=float(0))
@click.option("--Br1",type=click.FLOAT,default=float(0))
@click.option("--I1",type=click.FLOAT,default=float(0))
@click.option("--Fe1",type=click.FLOAT,default=float(0))
@click.option("--Cu1",type=click.FLOAT,default=float(0))
@click.option("--Na1",type=click.FLOAT,default=float(0))

@click.option("--C2",help="Number of carbons in Mod 2",type=click.FLOAT,default=float(0))
@click.option("--H2",type=click.FLOAT,default=float(1))
@click.option("--O2",type=click.FLOAT,default=float(1))
@click.option("--N2",type=click.FLOAT,default=float(0))
@click.option("--S2",type=click.FLOAT,default=float(0))
@click.option("--Si2",type=click.FLOAT,default=float(0))
@click.option("--F2",type=click.FLOAT,default=float(0))
@click.option("--Cl2",type=click.FLOAT,default=float(0))
@click.option("--Br2",type=click.FLOAT,default=float(0))
@click.option("--I2",type=click.FLOAT,default=float(0))
@click.option("--Fe2",type=click.FLOAT,default=float(0))
@click.option("--Cu2",type=click.FLOAT,default=float(0))
@click.option("--Na2",type=click.FLOAT,default=float(0))

@click.option("--C3",help="Number of carbons in Mod 3",type=click.FLOAT,default=float(0))
@click.option("--H3",type=click.FLOAT,default=float(1))
@click.option("--O3",type=click.FLOAT,default=float(1))
@click.option("--N3",type=click.FLOAT,default=float(0))
@click.option("--S3",type=click.FLOAT,default=float(0))
@click.option("--Si3",type=click.FLOAT,default=float(0))
@click.option("--F3",type=click.FLOAT,default=float(0))
@click.option("--Cl3",type=click.FLOAT,default=float(0))
@click.option("--Br3",type=click.FLOAT,default=float(0))
@click.option("--I3",type=click.FLOAT,default=float(0))
@click.option("--Fe3",type=click.FLOAT,default=float(0))
@click.option("--Cu3",type=click.FLOAT,default=float(0))
@click.option("--Na3",type=click.FLOAT,default=float(0))

@click.option("--G1",help="Initial guess DS1",type=click.FLOAT,default=float(0.5))
@click.option("--G1min",help="DS1 minimum",type=click.FLOAT,default=float(0))
@click.option("--G1max",help="DS1 maximum",type=click.FLOAT,default=float(3))

@click.option("--G2",help="Initial guess DS2",type=click.FLOAT,default=float(0.5))
@click.option("--G2min",help="DS2 minimum",type=click.FLOAT,default=float(0))
@click.option("--G2max",help="DS2 maximum",type=click.FLOAT,default=float(3))

@click.option("--G3",help="Initial guess DS3",type=click.FLOAT,default=float(0.5))
@click.option("--G3min",help="DS3 minimum",type=click.FLOAT,default=float(0))
@click.option("--G3max",help="DS3 maximum",type=click.FLOAT,default=float(3))

@click.option("--watermax",help="water maximum (percent) only with vary water",type=click.FLOAT,default=float(1))

@click.option("--cweight",help="Carbon weighting during minimization",type=click.FLOAT,default=float(1e5))
@click.option("--hweight",help="Hydrogen weighting during minimization",type=click.FLOAT,default=float(1))
@click.option("--nweight",help="Nitrogen weighting during minimization",type=click.FLOAT,default=float(1e6))
@click.option("--sweight",help="Sulfur weighting during minimization",type=click.FLOAT,default=float(1e6))
@click.option("--oweight",help="Oxygen weighting during minimization",type=click.FLOAT,default=float(1e5))

@click.option("--method", help="Minimization routine (nelder or leastsq recommended) see lmfit documentation",default="nelder")

def cli(humidity,carb,hydr,nitr,sulf,nmod,oxyg,sili,fluo,
        chlor,brom,iodi,iron,copp,sodium,varywater,oxywater,
        chain_ratio,C1,H1,O1,N1,S1,Si1,F1,Cl1,Br1,I1,Fe1,Cu1,Na1,
        C2,H2,O2,N2,S2,Si2,F2,Cl2,Br2,I2,Fe2,Cu2,Na2,
        C3,H3,O3,N3,S3,Si3,F3,Cl3,Br3,I3,Fe3,Cu3,Na3,
        G1,G1min,G1max,G2,G2min,G2max,G3,G3min,G3max,watermax,
        cweight,hweight,nweight,sweight,oweight,method):
    """This script will calculate the DS for modified cellulose samples based on (CHNS)
     elemental analysis data. Type help to see a list of required and optional arguments.
     
     CARB    Detected carbon concentration
     HYDR    Detected hydrogen concentration
     NITR    Detected nitrogen concentration"""
    
    C_det = carb / 100
    H_det = hydr / 100
    O_det = oxyg / 100
    N_det = nitr / 100
    S_det = sulf / 100
    Si_det = sili / 100
    Cl_det = chlor / 100
    Br_det = brom / 100
    F_det = fluo / 100
    I_det =  iodi / 100
    Fe_det = iron / 100
    Cu_det = copp / 100
    Na_det = sodium / 100

    """Atomic masses"""
    Cmass = float(12.011)
    Hmass = float(1.008)
    Omass = float(15.999)
    Nmass = float(14.007)
    Smass = float(32.06)
    Simass = float(28.0855)
    Fmass = float(18.998)
    Clmass = float(35.45)
    Brmass = float(79.904)
    Imass = float(126.904)
    Femass = float(55.845)
    Cumass = float(63.546)
    Namass = float(22.98976928)

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

        H_mod = h_cel - DS1 - DS2 - DS3 + (DS1*H1) + (DS2*H2) + (DS3*H3)
        O_mod = o_cel - DS1 - DS2 - DS3 + (DS1*O1) + (DS2*O2) + (DS3*O3)
        
        C_mod = c_cel + (DS1*C1) + (DS2*C2) + (DS3*C3)
        N_mod = (DS1*N1) + (DS2*N2) + (DS3*N3)
        S_mod = (DS1*S1) + (DS2*S2) + (DS3*S3)
        Si_mod = (DS1*Si1) + (DS2*Si2) + (DS3*Si3)
        Cl_mod = (DS1*Cl1) + (DS2*Cl2) + (DS3*Cl3)
        Br_mod = (DS1*Br1) + (DS2*Br2) + (DS3*Br3)
        F_mod = (DS1*F1) + (DS2*F2) + (DS3*F3)
        I_mod = (DS1*I1) + (DS2*I2) + (DS3*I3)
        Fe_mod = (DS1*Fe1) + (DS2*Fe2) + (DS3*Fe3)
        Cu_mod = (DS1*Cu1) + (DS2*Cu2) + (DS3*Cu3)
        Na_mod =  (DS1*Na1) + (DS2*Na2) + (DS3*Na3)
        
        
        
        """Time to get massive"""
        
        C_tot = C_mod * Cmass
        H_tot = H_mod * Hmass
        O_tot = O_mod * Omass
        N_tot = N_mod * Nmass
        S_tot = S_mod * Smass
        Si_tot = Si_mod * Simass
        Cl_tot = Cl_mod * Clmass
        Br_tot = Br_mod * Brmass
        F_tot = F_mod * Fmass
        I_tot = I_mod * Imass
        Fe_tot = Fe_mod * Femass
        Cu_tot = Cu_mod * Cumass
        Na_tot = Na_mod * Namass
        
        M_tot = C_tot + H_tot + O_tot + N_tot + S_tot + Si_tot + Cl_tot + Br_tot + F_tot + I_tot + Fe_tot + Cu_tot + Na_tot
        
        water_frac = parvals['water']
        
        if oxywater > 0:
            water_frac2 = oxywater/100
            varywater = False
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
        return np.array((1*residC*cweight,1*residH*hweight,1*residO*oweight,1*residN*nweight,
                        1*residS*sweight,1*residSi,1*residCl,1*residBr,1*residF,1*residI,1*residFe,
                        1*residCu,1*residNa))


    """Create model from residual and parameterize"""
    model1 = lmf.Model(residual)
    if nmod == 1:
        model1.set_param_hint('DS1',value=G1,min=G1min,max=G1max)
        model1.set_param_hint('DS2',value=0,vary=False)
        model1.set_param_hint('DS3',value=0,vary=False)
    elif nmod == 2:
        model1.set_param_hint('DS1',value=G1,min=G1min,max=G1max)
        model1.set_param_hint('DS2',value=G2,min=G2min,max=G2max)
        model1.set_param_hint('DS3',value=0,vary=False)
    else:
        model1.set_param_hint('DS1',value=G1,min=G1min,max=G1max)
        model1.set_param_hint('DS2',value=G2,min=G2min,max=G2max)
        model1.set_param_hint('DS3',value=G3,min=G3min,max=G3max)
    if varywater == True:
        model1.set_param_hint('water',value=humidity/100,min=max(0,(humidity-0.5)/100),
                            max=max((humidity+2)/100,watermax/100))
    else:
        model1.set_param_hint('water',value=humidity/100,vary=False)
    pars = model1.make_params()

    """Here's the magic line. leastsq is LM algorithm, others are available."""
    result = lmf.minimize(residual,pars,method=method)

    """Format results"""
    report = lmf.fit_report(result)
    DS = list(result.params.valuesdict().values())

    H_mod = h_cel - DS[0] - DS[1] - DS[2] + (DS[0]*H1) + (DS[1]*H2) + (DS[2]*H3)
    O_mod = o_cel - DS[0] - DS[1] - DS[2] + (DS[0]*O1) + (DS[1]*O2) + (DS[2]*O3)

    C_mod = c_cel + (DS[0]*C1) + (DS[1]*C2) + (DS[2]*C3)
    N_mod = (DS[0]*N1) + (DS[1]*N2) + (DS[2]*N3)
    S_mod = (DS[0]*S1) + (DS[1]*S2) + (DS[2]*S3)
    Si_mod = (DS[0]*Si1) + (DS[1]*Si2) + (DS[2]*Si3)
    Cl_mod = (DS[0]*Cl1) + (DS[1]*Cl2) + (DS[2]*Cl3)
    Br_mod = (DS[0]*Br1) + (DS[1]*Br2) + (DS[2]*Br3)
    F_mod = (DS[0]*F1) + (DS[1]*F2) + (DS[2]*F3)
    I_mod = (DS[0]*I1) + (DS[1]*I2) + (DS[2]*I3)

    Fe_mod = (DS[0]*Fe1) + (DS[1]*Fe2) + (DS[2]*Fe3)
    Cu_mod = (DS[0]*Cu1) + (DS[1]*Cu2) + (DS[2]*Cu3)
    Na_mod = (DS[0]*Na1) + (DS[1]*Na2) + (DS[2]*Na3)   

    """Time to get massive"""

    C_tot = C_mod * Cmass
    H_tot = H_mod * Hmass
    O_tot = O_mod * Omass
    N_tot = N_mod * Nmass
    S_tot = S_mod * Smass
    Si_tot = Si_mod * Simass
    Cl_tot = Cl_mod * Clmass
    Br_tot = Br_mod * Brmass
    F_tot = F_mod * Fmass
    I_tot = I_mod * Imass
    Fe_tot = Fe_mod * Femass
    Cu_tot = Cu_mod * Cumass
    Na_tot = Na_mod * Namass

    M_tot = C_tot + H_tot + O_tot + N_tot + S_tot + Si_tot + Cl_tot + Br_tot + F_tot + I_tot + Fe_tot + Cu_tot + Na_tot

    water_frac = DS[3] #This takes the result from the fit

    if oxywater >0:
        water_frac2 = oxywater/100
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

    H_mod1 = (DS[0]*H1) + (DS[1]*H2) + (DS[2]*H3)
    O_mod1 = (DS[0]*O1) + (DS[1]*O2) + (DS[2]*O3)
    C_mod1 = (DS[0]*C1) + (DS[1]*C2) + (DS[2]*C3)

    C_tot1 = C_mod1 * Cmass
    H_tot1 = H_mod1 * Hmass
    O_tot1 = O_mod1 * Omass

    M_tot1 = C_tot1 + H_tot1 + O_tot1 + N_tot + S_tot + Si_tot + Cl_tot + Br_tot + F_tot + I_tot + Fe_tot + Cu_tot + Na_tot

    Frac_mod = M_tot1 / M_tot

    printout = ('empirical formula: C'+str(round(C_mod,3))+" H"+str(round(H_mod,3))+' O'+str(round(O_mod,3))+' N'+str(round(N_mod,3))
        +' S'+str(round(S_mod,3))+' F'+str(round(F_mod,3))+' Cl'+str(round(Cl_mod,3))+' Br'+str(round(Br_mod,3))+' I'+str(round(I_mod,3))
        +' Fe'+str(round(Fe_mod,3))+' Cu'+str(round(Cu_mod,3))+' Na'+str(round(Na_mod,3))+os.linesep
        +'Chain ratio: '+str(chain_ratio)+os.linesep
        +'Mod 1: C'+str(round(C1,3))+' H'+str(round(H1,3))+' O'+str(round(O1,3))+' N'+str(round(N1,3))+' S'+str(round(S1,3))+' Si'+str(round(Si1,3))+' F'+str(round(F1,3))+' Cl'+str(round(Cl1,3))+' Br'+str(round(Br1,3))+' I'+str(round(I1,3))+' Fe'+str(round(Fe1,3))+' Cu'+str(round(Cu1,3))+' Na'+str(round(Na1,3))+os.linesep
        +'Mod 2: C'+str(round(C2,3))+' H'+str(round(H2,3))+' O'+str(round(O2,3))+' N'+str(round(N2,3))+' S'+str(round(S2,3))+' Si'+str(round(Si2,3))+' F'+str(round(F2,3))+' Cl'+str(round(Cl2,3))+' Br'+str(round(Br2,3))+' I'+str(round(I2,3))+' Fe'+str(round(Fe2,3))+' Cu'+str(round(Cu2,3))+' Na'+str(round(Na2,3))+os.linesep
        +'Mod 3: C'+str(round(C3,3))+' H'+str(round(H3,3))+' O'+str(round(O3,3))+' N'+str(round(N3,3))+' S'+str(round(S3,3))+' Si'+str(round(Si3,3))+' F'+str(round(F3,3))+' Cl'+str(round(Cl3,3))+' Br'+str(round(Br3,3))+' I'+str(round(I3,3))+' Fe'+str(round(Fe3,3))+' Cu'+str(round(Cu3,3))+' Na'+str(round(Na3,3))+os.linesep
        +'DS1 = '+str(round(DS[0],3))+os.linesep
        +'DS2 = '+str(round(DS[1],3))+os.linesep
        +'DS3 = '+str(round(DS[2],3))+os.linesep
        +'DS1(surf) = '+str(round(DS[0]/chain_ratio,3))+os.linesep
        +'DS2(surf) = '+str(round(DS[1]/chain_ratio,3))+os.linesep
        +'DS3(surf) = '+str(round(DS[2]/chain_ratio,3))+os.linesep
        +'Found percent water = '+str(humidity)+os.linesep
        +'Found percent water (oxygen) = '+str(oxywater)+os.linesep
        +'Calc percent water = '+str(DS[3]*100)+os.linesep
        +'Mass fraction modification = '+str(round(Frac_mod,3))+os.linesep
        +'Molecular mass modified AGU = '+str(round(M_tot,5))+os.linesep
        +'Found: '+os.linesep
        +'C, '+str(float(100*C_det))+'; H, '+str(float(100*H_det))+'; O, '+str(float(100*O_det))+'; N, '+str(float(100*N_det))
        +'; S, '+str(float(100*S_det))+'; Si, '+str(float(100*Si_det))+'; F, '+str(float(100*F_det))+'; Cl, '+str(float(100*Cl_det))
        +'; Br, '+str(float(100*Br_det))+'; I, '+str(float(100*I_det))+'; Fe, '+str(float(100*Fe_det))+'; Cu, '+str(float(100*Cu_det))+'; Na, '+str(float(100*Na_det))+os.linesep
        +'Calculated: '+os.linesep
        +'C, '+str(round(float(100*C_frac),3))+'; H, '+str(round(float(100*H_frac),3))+'; O, '+str(round(float(100*O_frac),3))+'; N, '+str(round(float(100*N_frac),3))
        +'; S, '+str(round(float(100*S_frac),3))+'; Si, '+str(round(float(100*Si_frac),3))+'; F, '+str(round(float(100*F_frac),3))+'; Cl, '+str(round(float(100*Cl_frac),3))
        +'; Br, '+str(round(float(100*Br_frac),3))+'; I, '+str(round(float(100*I_frac),3))+'; Fe, '+str(round(float(100*Fe_frac),3))+'; Cu, '+str(round(float(100*Cu_frac),3))+'; Na, '+str(round(float(100*Na_frac),3))
        +os.linesep+os.linesep+'Optimization data:'+os.linesep+os.linesep+str(report)
        +os.linesep+os.linesep+"Command line:"+os.linesep
        +os.linesep+str(sys.argv))

    with open('result.txt','w') as text_file:
        print(printout,file=text_file)
    print(printout)

if __name__ == '__main__':
    cli()