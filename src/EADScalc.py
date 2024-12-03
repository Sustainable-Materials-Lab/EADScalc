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
        chain_ratio,c1,h1,o1,n1,s1,si1,f1,cl1,br1,i1,fe1,cu1,na1,
        c2,h2,o2,n2,s2,si2,f2,cl2,br2,i2,fe2,cu2,na2,
        c3,h3,o3,n3,s3,si3,f3,cl3,br3,i3,fe3,cu3,na3,
        g1,g1min,g1max,g2,g2min,g2max,g3,g3min,g3max,watermax,
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

    # Atomic masses
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

    # Define cellulose

    c_cel = float(6)
    o_cel = float(5)
    h_cel = float(10)


    # Define residual funtion to minimize
    def residual (pars):
        parvals = pars.valuesdict()
        DS1 = parvals['DS1']
        DS2 = parvals['DS2']
        DS3 = parvals['DS3']

        H_mod = h_cel - DS1 - DS2 - DS3 + (DS1*h1) + (DS2*h2) + (DS3*h3)
        O_mod = o_cel - DS1 - DS2 - DS3 + (DS1*o1) + (DS2*o2) + (DS3*o3)
        
        C_mod = c_cel + (DS1*c1) + (DS2*c2) + (DS3*c3)
        N_mod = (DS1*n1) + (DS2*n2) + (DS3*n3)
        S_mod = (DS1*s1) + (DS2*s2) + (DS3*s3)
        Si_mod = (DS1*si1) + (DS2*si2) + (DS3*si3)
        Cl_mod = (DS1*cl1) + (DS2*cl2) + (DS3*cl3)
        Br_mod = (DS1*br1) + (DS2*br2) + (DS3*br3)
        F_mod = (DS1*f1) + (DS2*f2) + (DS3*f3)
        I_mod = (DS1*i1) + (DS2*i2) + (DS3*i3)
        Fe_mod = (DS1*fe1) + (DS2*fe2) + (DS3*fe3)
        Cu_mod = (DS1*cu1) + (DS2*cu2) + (DS3*cu3)
        Na_mod =  (DS1*na1) + (DS2*na2) + (DS3*na3)
        
        
        
        # Time to get massive
        
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
        
        # Life gets complicated!
    
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


    # Create model from residual and parameterize
    model1 = lmf.Model(residual)
    if nmod == 1:
        model1.set_param_hint('DS1',value=g1,min=g1min,max=g1max)
        model1.set_param_hint('DS2',value=0,vary=False)
        model1.set_param_hint('DS3',value=0,vary=False)
    elif nmod == 2:
        model1.set_param_hint('DS1',value=g1,min=g1min,max=g1max)
        model1.set_param_hint('DS2',value=g2,min=g2min,max=g2max)
        model1.set_param_hint('DS3',value=0,vary=False)
    else:
        model1.set_param_hint('DS1',value=g1,min=g1min,max=g1max)
        model1.set_param_hint('DS2',value=g2,min=g2min,max=g2max)
        model1.set_param_hint('DS3',value=g3,min=g3min,max=g3max)
    if varywater == True:
        model1.set_param_hint('water',value=humidity/100,min=max(0,(humidity-0.5)/100),
                            max=max((humidity+2)/100,watermax/100))
    else:
        model1.set_param_hint('water',value=humidity/100,vary=False)
    pars = model1.make_params()

    # Here's the magic line. leastsq is LM algorithm, others are available.
    result = lmf.minimize(residual,pars,method=method)

    # Format results
    report = lmf.fit_report(result)
    DS = list(result.params.valuesdict().values())

    H_mod = h_cel - DS[0] - DS[1] - DS[2] + (DS[0]*h1) + (DS[1]*h2) + (DS[2]*h3)
    O_mod = o_cel - DS[0] - DS[1] - DS[2] + (DS[0]*o1) + (DS[1]*o2) + (DS[2]*o3)

    C_mod = c_cel + (DS[0]*c1) + (DS[1]*c2) + (DS[2]*c3)
    N_mod = (DS[0]*n1) + (DS[1]*n2) + (DS[2]*n3)
    S_mod = (DS[0]*s1) + (DS[1]*s2) + (DS[2]*s3)
    Si_mod = (DS[0]*si1) + (DS[1]*si2) + (DS[2]*si3)
    Cl_mod = (DS[0]*cl1) + (DS[1]*cl2) + (DS[2]*cl3)
    Br_mod = (DS[0]*br1) + (DS[1]*br2) + (DS[2]*br3)
    F_mod = (DS[0]*f1) + (DS[1]*f2) + (DS[2]*f3)
    I_mod = (DS[0]*i1) + (DS[1]*i2) + (DS[2]*i3)

    Fe_mod = (DS[0]*fe1) + (DS[1]*fe2) + (DS[2]*fe3)
    Cu_mod = (DS[0]*cu1) + (DS[1]*cu2) + (DS[2]*cu3)
    Na_mod = (DS[0]*na1) + (DS[1]*na2) + (DS[2]*na3)   

    # Time to get massive

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

    # Extra for calc of mass percent mod

    H_mod1 = (DS[0]*h1) + (DS[1]*h2) + (DS[2]*h3)
    O_mod1 = (DS[0]*o1) + (DS[1]*o2) + (DS[2]*o3)
    C_mod1 = (DS[0]*c1) + (DS[1]*c2) + (DS[2]*c3)

    C_tot1 = C_mod1 * Cmass
    H_tot1 = H_mod1 * Hmass
    O_tot1 = O_mod1 * Omass

    M_tot1 = C_tot1 + H_tot1 + O_tot1 + N_tot + S_tot + Si_tot + Cl_tot + Br_tot + F_tot + I_tot + Fe_tot + Cu_tot + Na_tot

    Frac_mod = M_tot1 / M_tot

    printout = ('empirical formula: C'+str(round(C_mod,3))+" H"+str(round(H_mod,3))+' O'+str(round(O_mod,3))+' N'+str(round(N_mod,3))
        +' S'+str(round(S_mod,3))+' F'+str(round(F_mod,3))+' Cl'+str(round(Cl_mod,3))+' Br'+str(round(Br_mod,3))+' I'+str(round(I_mod,3))
        +' Fe'+str(round(Fe_mod,3))+' Cu'+str(round(Cu_mod,3))+' Na'+str(round(Na_mod,3))+os.linesep
        +'Chain ratio: '+str(chain_ratio)+os.linesep
        +'Mod 1: C'+str(round(c1,3))+' H'+str(round(h1,3))+' O'+str(round(o1,3))+' N'+str(round(n1,3))+' S'+str(round(s1,3))+' Si'+str(round(si1,3))+' F'+str(round(f1,3))+' Cl'+str(round(cl1,3))+' Br'+str(round(br1,3))+' I'+str(round(i1,3))+' Fe'+str(round(fe1,3))+' Cu'+str(round(cu1,3))+' Na'+str(round(na1,3))+os.linesep
        +'Mod 2: C'+str(round(c2,3))+' H'+str(round(h2,3))+' O'+str(round(o2,3))+' N'+str(round(n2,3))+' S'+str(round(s2,3))+' Si'+str(round(si2,3))+' F'+str(round(f2,3))+' Cl'+str(round(cl2,3))+' Br'+str(round(br2,3))+' I'+str(round(i2,3))+' Fe'+str(round(fe2,3))+' Cu'+str(round(cu2,3))+' Na'+str(round(na2,3))+os.linesep
        +'Mod 3: C'+str(round(c3,3))+' H'+str(round(h3,3))+' O'+str(round(o3,3))+' N'+str(round(n3,3))+' S'+str(round(s3,3))+' Si'+str(round(si3,3))+' F'+str(round(f3,3))+' Cl'+str(round(cl3,3))+' Br'+str(round(br3,3))+' I'+str(round(i3,3))+' Fe'+str(round(fe3,3))+' Cu'+str(round(cu3,3))+' Na'+str(round(na3,3))+os.linesep
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