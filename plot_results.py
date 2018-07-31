import matplotlib.pyplot as plt
from astropy.table import Table, Column, join, unique
import pandas as pd
import numpy as np
import os

def plot_mass_z_relation(pandas_df, plot_path='./'):

	#table for green pea sample 
	gp = Table.read('/Volumes/Briana_mac3/Green_Peas/green_pea_table.ascii', format='ascii') #Cardamone 2009
	gp_mass = gp['Stellar_mass'].data
	gp_metals = gp['metal'].data
	gp_mask = np.where(gp_metals !='\xe2\x8b\xaf')
	gp_metals = gp_metals[gp_mask]
	gp_metals = gp_metals.astype(np.float)
	gp_mass   = gp_mass[gp_mask]

	#table for blueberry sample 
	bb = Table.read('/Volumes/Briana_mac3/Green_Peas/blueberry_table.ascii', format='ascii') #Yang 2017
	bb_mass = bb['log(M/Msun)']
	bb_metals = bb['12+log(O/H)'].data
	bb_mask = np.where(bb_metals > 0)
	bb_metals = bb_metals[bb_mask]
	bb_mass = bb_mass[bb_mask]

	#table for blue compact dwarf sample 
	bcd = pd.read_csv('/Volumes/Briana_mac3/Green_Peas/blue_compact_dwarfs_table.ascii', delim_whitespace=True) #Lian 2016
	bcd_mass = bcd['logM'].values
	bcd_metals = bcd['MZ'].values

	#plot curve from SDSS SF galaxies 
	mass_range = np.linspace(7.4,10.5,40)
	mzr_a  = [7.8,8.0,8.1,8.3,8.4,8.5,8.7,8.75] 
	mass_a = [7.4,7.7,8.0,8.4,8.65,9.0,9.9,10.5]

	pol_fit = np.polyfit(mass_a, mzr_a, 2)
	p = np.poly1d(pol_fit)
	mzr_range = p(mass_range)

	#define data from pandas dataframe
	mass = pandas_df['logM'].values
	met_lis = pandas_df['Metal'].values
	m_u_lis = pandas_df['Metal_Uerr'].values
	m_l_lis = pandas_df['Metal_Lerr'].values

	met_err = np.divide(np.add(m_u_lis, m_l_lis),2)

	plt.scatter(bcd_mass, bcd_metals, facecolors='none', edgecolors='lightblue', label='blue compact dwarfs (Lian+ 2016)')
	plt.scatter(bb_mass, bb_metals, facecolors='none', edgecolors='blue', label='blueberries (Yang+ 2017)')
	plt.scatter(gp_mass, gp_metals, facecolors='none', edgecolors='green', label='green peas (Cardamone 2009)')
	plt.plot(mass_range, mzr_range, color='grey', label='SDSS SF galaxies (Andrews+ 2013)')
	plt.errorbar(mass, met_lis, yerr=(m_u_lis, m_l_lis), capsize=7, capthick=2, color='lightcoral', ls='None')
	plt.scatter(mass, met_lis, color='red', label='this work')
	plt.xlabel(r'$log(M_{*}/M_{\odot})$')
	plt.ylabel(r'$12+log(O/H)$')

	plt.legend()
	#plt.show()
	plt.savefig(os.path.join(plot_path, 'MCMC_mass_metallicity.png'))