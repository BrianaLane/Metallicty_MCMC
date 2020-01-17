import numpy as np
import pandas as pd
import os.path as op
#import plot_results as pr
import metallicity_emcee_class as met_mc

#***********************************#
# Load line fluxes and their errors #
#***********************************#

#import pandas dataframe with line fluxes and errors for objects 
#hps = pd.read_csv('/Volumes/Briana_mac3/HPS/HPS_OII_Fluxes.csv')
hps = pd.read_csv('/Users/Briana/Documents/Grad_School/HPS/HPS_dataframes/OII_norm_fluxes.csv')
print hps.columns

# hps=pd.DataFrame(index = np.arange(50), columns = hps.columns)
# hps.drop(['Unnamed: 0', 'ID', '[Ne]3870', '[OIII]_HPS', '[OIII]	_e_HPS'], axis=1, inplace=True)

# hps['HPS_name'] = np.arange(50)
# hps['Source'] = 'OII's
# hps['[OII]3727']   = 0.0
# hps['[OIII]5007']  = 0.0
# hps['[Hb]4861']    = 0.0
# hps['[NeIII]3870'] = 0.0
# hps['[Ha]6562']    = 0.0
# hps['[NII]6583']   = 0.0

# hps['Hb_absorption']  = 0.0

# hps['[OII]3727_e']   = 1.0
# hps['[OIII]5007_e']  = 1.0
# hps['[Hb]4861_e']    = 1.0
# hps['[NeIII]3870_e'] = 1.0
# hps['[Ha]6562_e']    = 1.0
# hps['[NII]6583_e']   = 1.0

#hps = pd.DataFrame(columns=hps.columns)

show_plot = False
save_plot = True
make_dataframe = True

omit_redlines = True
omit_Ne = True

model='maiolino'

#outpath = 'met_mcmc_results_OII/plots'
outpath = '/Users/Briana/Documents/Grad_School/HPS/Paper_Plots'
#met_df = '/Volumes/Briana_mac3/HPS/hps_metallicty_measures.csv'
#met_df = op.join(outpath, 'hps_metallicty_measures_norm_OII_test_w_NeIII.csv')
met_df = '/Users/Briana/Documents/Grad_School/HPS/HPS_dataframes/hps_metallicty_measures_norm_OII.csv'
#ss_out = '/Users/Briana/Documents/Grad_School/HPS/HPS_dataframes/sample_stack_'+model+'_OIII_w_RL.txt'

columns = ['HPS_name', 'metallicity', 'E(B-V)', 'OIII_intrinsic', 'met_err', 'E(B-V)_err', 'OIII_int_err']
if make_dataframe:
	df = pd.DataFrame(columns=columns)

else:
	df = pd.read_csv(met_df)

#assign the dtypes to be float for the line fluxes
#and objects for the errors since they will be a list of 2 values
dtypes = {k: object for k in columns[:1]}
dtypes.update({k: float for k in columns[2:4]})
dtypes.update({k: object for k in columns[5:]})
df = df.astype(dtypes)

x_lis = []
E_lis = []
sample_lis = []
#for i in range(len(hps)):
for i in [15]:
	object_name = hps['HPS_name'].iloc[i]
	print i, object_name


	OII   = hps['[OII]3727'].iloc[i]
	OIII  = hps['[OIII]5007'].iloc[i]
	Hb    = hps['[Hb]4861'].iloc[i]
	NeIII = hps['[NeIII]3870'].iloc[i]
	Ha    = hps['[Ha]6562'].iloc[i]
	NII   = hps['[NII]6583'].iloc[i]

	Hb_ab   = hps['Hb_absorption'].iloc[i]
	Hb_corr = Hb+Hb_ab

	OII_e   = hps['[OII]3727_e'].iloc[i]
	OIII_e  = hps['[OIII]5007_e'].iloc[i]
	Hb_e    = hps['[Hb]4861_e'].iloc[i]
	NeIII_e = hps['[NeIII]3870_e'].iloc[i]
	Ha_e    = hps['[Ha]6562_e'].iloc[i]
	NII_e   = hps['[NII]6583_e'].iloc[i]

	#NII = 0.001
	#NII_e = 10.0

	Hb_ab_e   = Hb_ab*0.3 #30% error on measurement
	Hb_corr_e = (np.sqrt((Hb_e**2) + (Hb_ab_e**2)))

	# print 'NII: ', NII, 'NII_e: ', NII_e, 'Ha: ', Ha, 'Ha_e: ', Ha_e
	# print 'NII_e: ', NII_e
	print 'OIII flux: ', OIII
	print 'OII: ', OII, 'OII_e: ', OII_e, 'OIII: ', OIII, 'OIII_e: ', OIII_e
	print 'Hb', Hb_corr, 'Hb_e', Hb_corr_e
	print 'N2 ratio: ', np.log10(NII/Ha)
	print 'O3O2 ratio: ', np.log10(OIII/OII)
	print 'R23: ', np.log10((((4/3)*OIII)+OII)/Hb_corr)

	if np.isnan(NeIII) or (omit_Ne == True) or (NeIII < 0.0):
		use_NeIII = False
		if np.isnan(Ha) or (omit_redlines == True):
			use_redlines = False
			data =  np.array([OII, Hb_corr, OIII])
			error = np.array([OII_e, Hb_corr_e, OIII_e])
		else:
			use_redlines = True
			print 'Using red lines'
			data =  np.array([OII, Hb_corr, OIII, Ha, NII])
			error = np.array([OII_e, Hb_corr_e, OIII_e, Ha_e, NII_e])

	else:
		print 'Using NeIII'
		use_NeIII = True
		if np.isnan(Ha) or (omit_redlines == True):
			use_redlines = False
			data =  np.array([OII, NeIII, Hb_corr, OIII])
			error = np.array([OII_e, NeIII_e, Hb_corr_e, OIII_e])
		else:
			use_redlines = True
			print 'Using red lines'
			data =  np.array([OII, NeIII, Hb_corr, OIII, Ha, NII])
			error = np.array([OII_e, NeIII_e, Hb_corr_e, OIII_e, Ha_e, NII_e])


	print 'DATA: ', data

	#***********#
	# Run emcee #
	#***********#
	Fitter = met_mc.FitterEmcee(object_name, data, error, use_NeIII=use_NeIII, use_redlines=use_redlines, model=model)
	Fitter.fit()
	results = Fitter.results
	samples = Fitter.samples
	sample_lis.append(samples)
	Fitter.plot_corner_results(outpath, show=show_plot, save=save_plot)
	#Fitter.plot_best_solution(outpath, show=show_plot, save=save_plot)
	#Fitter.plot_result_ratios(outpath, show=show_plot, save=save_plot)

	#**************#
	# Save Results #
	#**************#
	print results

	df.at[i, 'HPS_name'] = object_name

	df.at[i, 'metallicity']    = results[0][0]
	df.at[i, 'E(B-V)']         = results[1][0]
	df.at[i, 'OIII_intrinsic'] = results[2][0]

	df.at[i,'met_err']      = [results[0][1], results[0][2]]
	df.at[i,'E(B-V)_err']   = [results[1][1], results[1][2]]
	df.at[i,'OIII_int_err'] = [results[2][1], results[2][2]]

#sample_stack = np.hstack(sample_lis)
#np.savetxt(ss_out, sample_stack)

#df.to_csv(met_df)
#print df.head(5)
#print df['metallicity']

		#**************#
		# Plot Results #
		#**************#
		# if show_plots or save_plots:
		# 	if not os.path.isdir(os.path.join(plot_path, 'corner_plots')):
		# 		os.mkdir(os.path.join(plot_path, 'corner_plots' ))
		# 		os.mkdir(os.path.join(plot_path, 'best_solution_plots' ))
		# 		os.mkdir(os.path.join(plot_path, 'result_ratio_plots' ))
		# 	mc_plots.build_corner_plot(object_name, flat_samples, x_mcmc, E_mcmc, os.path.join(plot_path, 'corner_plots'), show=show_plots, save=save_plots)
		# 	mc_plots.plot_best_solution(object_name, flat_samples, x_mcmc, E_mcmc, os.path.join(plot_path, 'best_solution_plots'), show=show_plots, save=save_plots)
		# 	mc_plots.plot_result_ratios(object_name, flat_samples, x_mcmc, E_mcmc, args, os.path.join(plot_path, 'result_ratio_plots'), show=show_plots, save=save_plots)

