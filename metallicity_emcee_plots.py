import numpy as np
import os.path as op
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import corner
import Maiolino_relation_functions as mrf

solar_x = 8.69
disp_dict = {'R23':0.03771, 'O32':0.16025, 'O3Hb':0.06760, 'NeO2':0.14452, 'O2Hb':0.10521}


def build_corner_plot(object_name, flatchain, mc_x, mc_E_bv, path, show=False, save=True):
	fig = corner.corner(flatchain, labels=["$X$", "$E(b-v)$"], truths=[mc_x[0], mc_E_bv[0]])
	if show:
		plt.show()
	if save:
		fig.savefig(op.join(path, object_name+"_cornerplot.png"))

	plt.close()

def plot_best_solution(object_name, flatchain, mc_x, mc_E_bv, path, show=False, save=True):
	x_samples = flatchain[:,0]
	E_bv_samples = flatchain[:,1]

	#make the confidience interval plots seen in the maiolino paper
	cov = np.cov(E_bv_samples, x_samples)
	lambda_, v = np.linalg.eig(cov)
	lambda_ = np.sqrt(lambda_)

	ax = plt.subplot(111)
	for j in xrange(1, 4):
		ell = Ellipse(xy=(mc_E_bv[0], mc_x[0]),
						width=lambda_[0]*j*2, height=lambda_[1]*j*2,
						angle=np.rad2deg(np.arccos(v[0, 0])))
		ell.set_facecolor('none')
		ell.set_edgecolor('red')
		ax.add_artist(ell)

	plt.scatter(mc_E_bv[0], mc_x[0])
	plt.xlabel('E(b-v)')
	plt.ylabel('X')
	plt.xlim(0.0, 0.60)
	plt.ylim(7.0, 10.0)
	plt.text(0.5, 7.3, 'X: '+str(round(mc_x[0],2))+'\n'+'E(B-V): '+str(round(mc_E_bv[0],2)), 
		horizontalalignment='center', verticalalignment='center', bbox=dict(facecolor='red', alpha=0.3)) 
	#plt.text(x, y, s, bbox=dict(facecolor='red', alpha=0.5))

	if show:
		plt.show()
	if save:
		plt.savefig(op.join(path, object_name+"_best_result.png"))

	plt.close()

def plot_result_ratios(object_name, flatchain, mc_x, mc_E_bv, args, path, show=False, save=True):
	met = np.arange(7.5, 9.3, 0.1)
	met_norm = np.subtract(met, solar_x)
	model_mc_x = mc_x[0] - solar_x

	OIII    = args[0]
	OII     = args[1]
	Hb      = args[2]
	NeIII   = args[3]
	OIII_e  = args[4]
	OII_e   = args[5]
	Hb_e    = args[6]
	NeIII_e = args[7]

	x_samples = flatchain[:,0]
	E_bv_samples = flatchain[:,1]

	ratio_dict = {'R23':[mrf.R23_model, 
						mrf.R23_ratio(OIII, OII, Hb, E_bv_samples), 
						mrf.R23_ratio(OIII, OII, Hb, mc_E_bv[0]), 
						mrf.R23_ratio(OIII, OII, Hb, 0), 
						mrf.R23_ratio_err(OIII, OII, Hb, OIII_e, OII_e, Hb_e, 0)], 
				'O32':[mrf.RO32_model, 
						mrf.RO32_ratio(OIII, OII, E_bv_samples),
						mrf.RO32_ratio(OIII, OII, mc_E_bv[0]),
						mrf.RO32_ratio(OIII, OII, 0), 
						mrf.RO32_ratio_err(OIII, OII, OIII_e, OII_e, 0)],
				'O3Hb':[mrf.RO3Hb_model, 
						mrf.RO3Hb_ratio(OIII, Hb, E_bv_samples),
						mrf.RO3Hb_ratio(OIII, Hb, mc_E_bv[0]),
						mrf.RO3Hb_ratio(OIII, Hb, 0), 
						mrf.RO3Hb_ratio_err(OIII, Hb, OIII_e, Hb_e, 0)]}

	if not np.isnan(NeIII):

		ratio_dict['NeO2'] = [mrf.RNeO2_model, 
							mrf.RNeO2_ratio(NeIII, OII, E_bv_samples),
							mrf.RNeO2_ratio(NeIII, OII, mc_E_bv[0]),
							mrf.RNeO2_ratio(NeIII, OII, 0), 
							mrf.RNeO2_ratio_err(NeIII, OII, NeIII_e, OII_e, 0)]

	f, ax = plt.subplots(1,len(ratio_dict), sharex=True, figsize=(18, 6))
	f.suptitle('Results')

	for i, ratio in enumerate(ratio_dict):
		#plot the model curves with disperson 
		model_log_ratio = np.log10(ratio_dict[ratio][0](met_norm))
		model_log_disp  = disp_dict[ratio]

		ax[i].plot(met, model_log_ratio, color='black')
		ax[i].plot(met, model_log_ratio+model_log_disp, ls=':', color='black')
		ax[i].plot(met, model_log_ratio-model_log_disp, ls=':', color='black')

		#plot the values derived from the model and the observed values
		mc_model_ratio     = np.log10(ratio_dict[ratio][0](model_mc_x))
		mc_obs_ratio       = np.log10(ratio_dict[ratio][3])
		mc_obs_ratio_deRed = np.log10(ratio_dict[ratio][2])
		mc_obs_ratio_err   = ratio_dict[ratio][4]
		mc_log_ratio_err   = (np.log10((10**mc_obs_ratio)+mc_obs_ratio_err)-np.log10((10**mc_obs_ratio)-mc_obs_ratio_err))/2

		ax[i].errorbar(mc_x[0], mc_obs_ratio, xerr=(mc_x[1]+mc_x[2])/2, yerr= mc_log_ratio_err, capsize=7, capthick=2, color='green')
		ax[i].scatter(mc_x[0], mc_obs_ratio_deRed , marker='x', s=100, color='blue')

		#compute the covariance matrix of the ratios and the metallicties 
		#this is based on the values for X and E_bv
		ratio_samples = np.log10(ratio_dict[ratio][1])

		cov = np.cov(x_samples, ratio_samples)
		lambda_, v = np.linalg.eig(cov)
		lambda_ = np.sqrt(lambda_)

		ell = Ellipse(xy=(mc_x[0], mc_obs_ratio_deRed),width=lambda_[0]*2, height=lambda_[1]*2, angle=np.rad2deg(np.arccos(v[0, 0])))
		ell.set_facecolor('none')
		ell.set_edgecolor('red')
		ax[i].add_artist(ell)

		#plot the axis labels
		ax[i].set_xlabel('12+log(O/H)')
		ax[i].set_ylabel('log '+ ratio)

	if show:
		plt.show()
	if save:
		plt.savefig(op.join(path, object_name+"_result_ratios.png"))
		
	plt.close()

