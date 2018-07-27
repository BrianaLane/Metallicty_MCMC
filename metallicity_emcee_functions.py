import numpy as np
import emcee
import math
import Maiolino_relation_functions as mrf
import metallicity_emcee_plots as mc_plots

solar_x = 8.69
disp_dict = {'R23':0.03771, 'O32':0.16025, 'O3Hb':0.06760, 'NeO2':0.14452, 'O2Hb':0.10521}

#define the log likelihood function
def lnlike(theta, OII, OIII, Hb, NeIII, OII_e, OIII_e, Hb_e, NeIII_e):
	x = theta[0] - solar_x
	E_bv = theta[1]

	RO32_mod = mrf.RO32_model(x)
	O3Hb_mod  = mrf.RO3Hb_model(x)

	RO32_obs = mrf.RO32_ratio(OIII, OII, E_bv)
	O3Hb_obs  = mrf.RO3Hb_ratio(OIII, Hb, E_bv)

	RO32_mod_var = (10**disp_dict['O32'])**2
	O3Hb_mod_var  = (10**disp_dict['O3Hb'])**2

	RO32_obs_var = mrf.RO32_ratio_err(OIII, OII, OIII_e, OII_e, E_bv)
	O3Hb_obs_var = mrf.RO3Hb_ratio_err(OIII, Hb, OIII_e, Hb_e, E_bv)

	if np.isnan(NeIII):
		return -0.5*((((O3Hb_obs-O3Hb_mod)**2)/np.sqrt(O3Hb_obs_var+O3Hb_mod_var))
					+(((RO32_obs-RO32_mod)**2)/np.sqrt(RO32_obs_var+RO32_mod_var)))

	else:
		NeO2_mod = mrf.RNeO2_model(x)
		NeO2_obs = mrf.RNeO2_ratio(NeIII, OII, E_bv)
		NeO2_mod_var  = (10**disp_dict['NeO2'])**2
		NeO2_obs_var = mrf.RNeO2_ratio_err(NeIII, OII, NeIII_e, OII_e, E_bv)

		return -0.5*((((O3Hb_obs-O3Hb_mod)**2)/np.sqrt(O3Hb_obs_var+O3Hb_mod_var))
					+(((RO32_obs-RO32_mod)**2)/np.sqrt(RO32_obs_var+RO32_mod_var))
					+(((NeO2_obs-NeO2_mod)**2)/np.sqrt(NeO2_obs_var+NeO2_mod_var)))

#define the log prior function
def lnprior(theta):
	x = theta[0] - solar_x
	E_bv = theta[1]

	#between metallicity 6.5 and 10
	#range of E_bv came from Joanna et. al. 2008
	if (-2.19 < x < 1.31) and (0.01 <= E_bv < 0.60):
		return 0.0

	return -np.inf

#define log postierior to sovle with emcee
def lnprob(theta, OII, OIII, Hb, NeIII, OII_e, OIII_e, Hb_e, NeIII_e):
	lp = lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(theta, OII, OIII, Hb, NeIII, OII_e, OIII_e, Hb_e, NeIII_e)

#set up and solve the model with MCMC 
def solve_model_emcee(lnprob, ndim, nwalkers, nchains, thetaGuess, args):
	pos = [thetaGuess + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)

	print "Burning in ..."
	pos, prob, state = sampler.run_mcmc(pos, nchains[0])

	# Reset the chain to remove the burn-in samples.
	sampler.reset()

	# Starting from the final position in the burn-in chain, sample for 1000
	# steps. (rstate0 is the state of the internal random number generator)
	print "Running MCMC ..."
	pos, prob, state = sampler.run_mcmc(pos, nchains[1], rstate0=state)

	flat_samples = sampler.flatchain
	samples = sampler.chain[:, :, :].reshape((-1, ndim))

	x_mcmc, E_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
						zip(*np.percentile(samples, [16, 50, 84], axis=0)))

	return flat_samples, x_mcmc, E_mcmc

def find_metallicity(pandas_df, col_dict, theta_guess, ndim, nwalkers, nchains, show_plots=False, save_plots=True):

	x_lis = []
	E_lis = []
	for i in range(len(pandas_df)):
		object_name = pandas_df[col_dict['Obj_names']].iloc[i]
		print i, object_name

		flux = [pandas_df[col_dict['OII_flux']].iloc[i], pandas_df[col_dict['OIII_flux']].iloc[i], 
				pandas_df[col_dict['Hb_flux']].iloc[i], pandas_df[col_dict['NeIII_flux']].iloc[i]] 
		S_N = [pandas_df[col_dict['OII_SN']].iloc[i], pandas_df[col_dict['OIII_SN']].iloc[i], 
				pandas_df[col_dict['Hb_SN']].iloc[i], pandas_df[col_dict['NeIII_SN']].iloc[i]]
		flux_err = np.divide(flux, S_N)

		OII   = flux[0]
		OIII  = flux[1]
		Hb    = flux[2]
		NeIII = flux[3]

		OII_e   = flux_err[0]
		OIII_e  = flux_err[1]
		Hb_e    = flux_err[2]
		NeIII_e = flux_err[3]

		args = (OII, OIII, Hb, NeIII, OII_e, OIII_e, Hb_e, NeIII_e)
		if not np.isnan(NeIII):
			print 'Using NeIII'

		#*********************#
		# Run emcee and plots #
		#*********************#
		flat_samples, x_mcmc, E_mcmc = solve_model_emcee(lnprob, ndim, nwalkers, nchains, theta_guess, args)
		x_lis.append(x_mcmc)
		E_lis.append(E_mcmc)
		print x_mcmc[0], E_mcmc[0]

		if show_plots or save_plots:
			#mc_plots.build_corner_plot(object_name, flat_samples, x_mcmc, E_mcmc, show=show_plots, save=save_plots)
			#mc_plots.plot_best_solution(object_name, flat_samples, x_mcmc, E_mcmc, show=show_plots, save=save_plots)
			mc_plots.plot_result_ratios(object_name, flat_samples, x_mcmc, E_mcmc, args, show=show_plots, save=save_plots)

		print '\n'

	pandas_df['Metal'] = [x[0] for x in x_lis]
	pandas_df['Metal_Uerr'] = [x[1] for x in x_lis]
	pandas_df['Metal_Lerr'] = [x[2] for x in x_lis]
	pandas_df['E(B-V)_mc'] = [e[0] for e in E_lis]
	pandas_df['E(B-V)_mc_Uerr'] = [e[1] for e in E_lis]
	pandas_df['E(B-V)_mc_Lerr'] = [e[2] for e in E_lis]

	return pandas_df


