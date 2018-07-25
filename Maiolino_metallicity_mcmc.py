import numpy as np
import Maiolino_relation_functions as mrf
import metallicity_emcee_functions as mef
import metallicity_emcee_plots as mc_plots

#******************#
# define constants #
#******************#
avg_E_bv = 0.21 #this is derived from Bridge et. al. for the HPS OII emitters

#***************************************************#
# define an observations of fluxes and their errors #
#***************************************************#
#in order [OII], [OIII]5007, [Hb]4861
# [OIII]4363 - line not needed since ratio fixed with 5007 - accounted for in model equation. 
object_name = 'HPS_test'

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

#*************#
# emcee stuff #
#*************#
thetaGuess = [8.45, avg_E_bv]
ndim, nwalkers = 2, 100
nchains = (200, 500) #200 for burn in and 500 after
args = (OII, OIII, Hb, OII_e, OIII_e, Hb_e)

flat_samples, x_mcmc, E_mcmc = mef.solve_model_emcee(mef.lnprob, ndim, nwalkers, nchains, thetaGuess, args)

mc_plots.build_corner_plot(object_name, flat_samples, x_mcmc, E_mcmc, show=True, save=False)
mc_plots.plot_best_solution(object_name, flat_samples, x_mcmc, E_mcmc, show=True, save=False)
mc_plots.plot_result_ratios(object_name, flat_samples, x_mcmc, E_mcmc, args, show=True, save=False)

