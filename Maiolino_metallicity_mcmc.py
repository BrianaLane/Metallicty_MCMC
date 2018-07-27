import numpy as np
import Maiolino_relation_functions as mrf
import metallicity_emcee_functions as mef
import metallicity_emcee_plots as mc_plots
import pandas as pd

#***********************************#
# Load line fluxes and their errors #
#***********************************#

hps = pd.read_csv('HPS_cat_table.dat')

col_dict = {'Obj_names':'HPS_name',
			'OII_flux':'OII_Flux',
			'OIII_flux':'OIII1_Flux',
			'Hb_flux':'Hb_Flux',
			'NeIII_flux':'NeIII1_Flux',
			'OII_SN':'OII_SN',
			'OIII_SN':'OIII1_SN',
			'Hb_SN':'Hb_SN',
			'NeIII_SN':'NeIII1_SN'}

hps = hps.dropna(subset=['OII_Flux', 'OIII1_Flux', 'Hb_Flux'])
print 'Found '+str(len(hps))+' objects with [OII], [OIII], and Hb'

#************************#
# Define emcee variables #
#************************#
avg_E_bv = hps['E(B-V)'].mean() #this is derived from Bridge et. al. for the HPS OII emitters

thetaGuess = [8.45, avg_E_bv]
ndim, nwalkers = 2, 100
nchains = (200, 500) #200 for burn in and 500 after

hps = mef.find_metallicity(hps, col_dict, thetaGuess, ndim, nwalkers, nchains, show_plots=False, save_plots=False)