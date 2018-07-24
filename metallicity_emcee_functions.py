import numpy as np
import emcee
import math
import Maiolino_relation_functions as mrf

#define the log likelihood function
def lnlike(theta, OII, OIII, Hb, OII_e, OIII_e, Hb_e):
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

	return -0.5*((((O3Hb_obs-O3Hb_mod)**2)/np.sqrt(O3Hb_obs_var+O3Hb_mod_var))+(((RO32_obs-RO32_mod)**2)/np.sqrt(RO32_obs_var+RO32_mod_var)))

#define the log prior function
def lnprior(theta):
	x = theta[0] - solar_x
	E_bv = theta[1]

	if (-2.19 < x < 1.31) and (0.01 <= E_bv < 0.60):
		return 0.0

	return -np.inf

#define log postierior to sovle with emcee
def lnprob(theta, OII, OIII, Hb, OII_e, OIII_e, Hb_e):
	lp = lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(theta, OII, OIII, Hb, OII_e, OIII_e, Hb_e)

#set up and solve the model with MCMC 
def solve_model_emcee(lnprob, ndim, nwalkers, thetaGuess, args):
	pos = thetaGuess + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)

	print "Burning in ..."
	pos, prob, state = sampler.run_mcmc(pos, 200)

	# Reset the chain to remove the burn-in samples.
	sampler.reset()

	# Starting from the final position in the burn-in chain, sample for 1000
	# steps. (rstate0 is the state of the internal random number generator)
	print "Running MCMC ..."
	pos, prob, state = sampler.run_mcmc(pos, 500, rstate0=state)

	flat_samples = sampler.flatchain

	samples = sampler.chain[:, :, :].reshape((-1, ndim))

	x_mcmc, E_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
						zip(*np.percentile(samples, [16, 50, 84], axis=0)))

	return flat_samples, x_mcmc, E_mcmc