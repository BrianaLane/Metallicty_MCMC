#*********************#
# Reddening Functions #
#*********************#

#define Calzetti law 2000 
#wl in angstroms
def calzetti(wl):
	wl = wl/10000.
	if wl <= 6300:
		return 2.695*(-2.156+(1.509/wl)-(0.189/(wl**2))+(0.011/(wl**3)))+4.05
	else:
		return 2.659*(-1.857+(1.040/wl))+4.05

def reddening(wl, E_bv):
	return 10**(0.4*E_bv*calzetti(wl))

#***********************#
# Ratio Model Functions #
#***********************#
def RO32_model(x):
	r_O32 = -0.2839+(-1.3881*x)+(-0.3172*(x**2))
	return 10**r_O32

def R23_model(x):
	r_23 = 0.7462+(-0.7149*x)+(-0.9401*(x**2))+(-0.6154*(x**3))+(-0.2524*(x**4))
	return 10**r_23

def RO3Hb_model(x):
	o3hb = 0.1549+(-1.5031*x)+(-0.9790*(x**2))+(-0.0297*(x**3))
	return 10**o3hb

def RNeO2_model(x):
	neo2 = -1.2608+(-1.0861*x)+(-0.1470*(x**2))
	return 10**neo2

def RO2Hb_model(x):
	o2hb = 0.5603+(0.0450*x)+(-1.8017*(x**2))+(-1.8434*(x**3))+(-0.6549*(x**4))
	return 10**o2hb

#**************************#
# Observed Ratio Functions #
#**************************#
def RO32_ratio(OIII, OII, E_bv):
	return (OIII*reddening(5007, E_bv))/(OII*reddening(3727, E_bv))
    
def R23_ratio(OIII, OII, Hb, E_bv):
	return ((OII*reddening(3727, E_bv))+((OIII*reddening(5007, E_bv))*(1+(1/2.98))))/(Hb*reddening(4861, E_bv))
    
def RO3Hb_ratio(OIII, Hb, E_bv):
	return (OIII*reddening(5007, E_bv))/(Hb*reddening(4861, E_bv))

def RNeO2_ratio(NeIII, OII, E_bv):
	return (NeIII*reddening(3869, E_bv))/(OII*reddening(3727, E_bv))

def RO2Hb_ratio(OII, Hb, E_bv):
	return (OII*reddening(3727, E_bv))/(Hb*reddening(4861, E_bv))

#********************************#
# Observed Ratio Error Functions #
#********************************#

def RO32_ratio_err(OIII, OII, OIII_e, OII_e, E_bv):
	e_5007 = reddening(5007, E_bv)
	e_3727 = reddening(3727, E_bv)
	return (((e_5007*OIII*OII_e)/((OII**2)*e_3727))**2)+(((e_5007*OIII_e)/(OII*e_3727))**2)

def R23_ratio_err(OIII, OII, Hb, OIII_e, OII_e, Hb_e, E_bv):
	e_5007 = reddening(5007, E_bv)
	e_3727 = reddening(3727, E_bv)
	e_4861 = reddening(4861, E_bv)
	return (((e_3727*OII_e)/(Hb*e_4861))**2)+((((1+(1/2.98))*e_5007*OIII_e)/(Hb*e_4861))**2)+((((1+(1/2.98))*e_5007*OIII*Hb_e)/((Hb**2)*e_4861))**2)

def RO3Hb_ratio_err(OIII, Hb, OIII_e, Hb_e, E_bv):
	e_5007 = reddening(5007, E_bv)
	e_4861 = reddening(4861, E_bv)
	return (((e_5007*OIII*Hb_e)/((Hb**2)*e_4861))**2)+(((e_5007*OIII_e)/(Hb*e_4861))**2)

def RNeO2_ratio_err(NeIII, OII, NeIII_e, OII_e, E_bv):
	e_3727 = reddening(3727, E_bv)
	e_3869 = reddening(3869, E_bv)
	return (((e_3869*NeIII*OII_e)/((OII**2)*e_3727))**2)+(((e_3869*NeIII_e)/(OII*e_3727))**2)

def RO2Hb_ratio_err(OII, Hb, OII_e, Hb_e, E_bv):
	e_3727 = reddening(3727, E_bv)
	e_4861 = reddening(4861, E_bv)
	return (((e_3727*OII*Hb_e)/((Hb**2)*e_4861))**2)+(((e_3727*OII_e)/(Hb*e_4861))**2)