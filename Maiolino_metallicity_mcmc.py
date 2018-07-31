import numpy as np
import pandas as pd
import Maiolino_relation_functions as mrf
import metallicity_emcee_functions as mef
import plot_results as pr

#***********************************#
# Load line fluxes and their errors #
#***********************************#

#import pandas dataframe with line fluxes and S/N for objects 
hps = pd.read_csv('HPS_cat_table.dat')

#Define the column names with the needed information
col_dict = {'Obj_names':'HPS_name',
			'OII_flux':'OII_Flux',
			'OIII_flux':'OIII1_Flux',
			'Hb_flux':'Hb_Flux',
			'NeIII_flux':'NeIII1_Flux',
			'OII_SN':'OII_SN',
			'OIII_SN':'OIII1_SN',
			'Hb_SN':'Hb_SN',
			'NeIII_SN':'NeIII1_SN'}

#Only keep rows that have the minimum needed line fluxes: [OIII]5007, [OII]3727, Hb
hps = hps.dropna(subset=[col_dict['OII_flux'], col_dict['OIII_flux'], col_dict['Hb_flux']])
print 'Found '+str(len(hps))+' objects with [OII], [OIII], and Hb'

#************************#
# Define emcee variables #
#************************#
avg_E_bv = hps['E(B-V)'].mean() #this is derived from Bridge et. al. for the HPS OII emitters

thetaGuess = [8.45, avg_E_bv] #guess values for metallicity and E(B-V)
ndim, nwalkers = 2, 100 #number of dimensions (len(thetaGuess)), and number of walkers 
nchains = (200, 500) #number of MC chains (200 for burn in and 500 after)

#***************************#
# Run MCMC and save results #
#***************************#
save_path = 'met_mcmc_results'
hps = mef.find_metallicity(hps, col_dict, thetaGuess, ndim, nwalkers, nchains, show_plots=False, save_plots=True, plot_path=save_path)
#hps.to_csv('HPS_cat_table.dat')
pr.plot_mass_z_relation(hps, plot_path=save_path)

