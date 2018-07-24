import numpy as np
import Maiolino_relation_functions as mrf
import metallicity_emcee_functions as mef
import metallicity_emcee_plots as mc_plots

#define constants
solar_x = 8.69
avg_E_bv = 0.21 #this is derived from Bridge et. al. for the HPS OII emitters

disp_dict = {'R23':0.03771, 'O32':0.16025, 'O3Hb':0.06760, 'NeO2':0.14452, 'O2Hb':0.10521}

#define an observations of fluxes and their errors 
#in order [OII], [OIII]5007, [Hb]4861
# [OIII]4363 - line not needed since ratio fixed with 5007 - accounted for in model equation. 

flux = [243.8, 67.0, 24.8] #HPS030638+000015
S_N = [23.5, 5.5, 6.8]
flux_err = np.divide(flux, S_N)
print flux_err

OII  = flux[0]
OIII = flux[1]
Hb   = flux[2]

OII_e  = flux_err[0]
OIII_e = flux_err[1]
Hb_e   = flux_err[2]

#emcee stuff
thetaGuess = [8.48547306, 0.12]
ndim, nwalkers = 2, 100
args = (OII, OIII, Hb, OII_e, OIII_e, Hb_e)

