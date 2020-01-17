# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 13:17:25 2018

@author: gregz, bindahl
"""

import numpy as np
import emcee
import corner
import time
import os.path as op
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
import Maiolino_relation_functions as mrf


class FitterEmcee:
    def __init__(self, obj_name, data, error, use_NeIII=False, use_redlines=False, model='maiolino', nsteps=(200,1000), nwalkers=100):
        ''' Initialize variables '''
        self.solar_x = 8.69
        self.obj_name = str(obj_name)
        self.model  = model #maiolino, curti 
        self.use_redlines = use_redlines

        if (self.model == 'curti') and (use_NeIII == True):
        #if using curti models need to remove NeIII because there is no realiton for it 
            self.data  = np.delete(data, [1])
            self.error = np.delete(error, [1])
            self.use_NeIII = False
        else:
            self.data = data * 1.
            self.error = error * 1.
            self.use_NeIII = use_NeIII

        if self.use_redlines:
            self.lims = [[7., 10.], [-0.05, 1.], [0., 1e6], [0., 1e6]]
            if self.use_NeIII:
                self.calz_fac = self.calzettilaw(np.array([3727., 3869., 4861., 5007., 6562., 6583.]))
            else:
                self.calz_fac = self.calzettilaw(np.array([3727., 4861., 5007., 6562., 6583.]))
                
        else:
            self.lims = [[7., 10.], [-0.05, 1.], [0., 1e6]]
            if self.use_NeIII:
                self.calz_fac = self.calzettilaw(np.array([3727., 3869., 4861., 5007.]))
            else:
                self.calz_fac = self.calzettilaw(np.array([3727., 4861., 5007.]))

        self.ndim = len(self.lims)
        self.nsteps = nsteps
        self.nwalkers = nwalkers

        self.truths = None
        self.samples = None
        self.flatchain = None
        self.results = None 

        #self.model_errs = {'R23':0.03771, 'O3O2':0.16025, 'O3Hb':0.06760, 'Ne3O2':0.14452, 'O2Hb':0.10521}
        # if self.model == 'maiolino':
        #     self.model_errs = {'R23':0.1, 'O3O2':0.1, 'Ne3O2':0.1} #assume all models == 0.1 , ~10% error in model predictions
        # elif self.model == 'curti':
        #     self.model_errs = {'R23':0.12, 'O3O2':0.14} 

        self.sigma_m = 0.3  # ~30% error in model predictions

    def calzettilaw(self, wave):
        ''' Calzetti et al. (2000) dust attenuation curve, k(wave)

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength

        Returns
        -------
        k : numpy array (1 dim)
            A(wave) = R_V * k(wave)
        '''
        invwv = 1/(wave/1e4)
        sel1 = np.nonzero(wave < 0.63e4)[0]
        sel2 = np.nonzero(np.logical_and(wave >= 0.63e4, wave < 2.2e4))[0]
        k1 = np.zeros(sel1.shape)
        k2 = np.zeros(sel2.shape)
        k1 = (2.659 * (-2.156 + 1.509 * invwv[sel1] - 0.198 * invwv[sel1]**2 +
              0.011 * invwv[sel1]**3) + 4.05)
        k2 = 2.659*(-1.857 + 1.040*invwv[sel2]) + 4.05
        k = np.zeros(wave.shape)
        k[sel1] = k1
        k[sel2] = k2
        return k

    # def get_model_err(self, thetas):
    #     O3O2_err  = self.model_errs['O3O2']
    #     R23_err   = self.model_errs['R23']

    #     O2_err = 1. / 10**O3O2_err
    #     HB_err = (4./3. + O2_err) / 10**R23_err

    #     if self.use_NeIII:
    #         Ne3O2_err = self.model_errs['Ne3O2']
    #         Ne3_err = 10**Ne3O2_err * O2_err

    #         model_err = np.array([O2_err, Ne3_err, HB_err, 1.])
    #     else:
    #         model_err = np.array([O2_err, HB_err, 1.])

    #     return model_err

    def get_model_fluxes(self, thetas):
        ''' Get Fluxes from model '''
        x = thetas[0] - self.solar_x
        Abv = self.calz_fac * thetas[1]

        if self.model == 'maiolino':
            O3O2 = -.2839 - 1.3881 * x - 0.3172 * x**2
            R23 = 0.7462 - 0.7149 * x - 0.9401 * x**2 - 0.6154 * x**3 - 0.2524 * x**4
            N2   = -0.7732 + 1.2357 * x - 0.2811 * x**2 - 0.7201 * x**3 - 0.3330 * x**4
            O3N2 = 0.4520 - 2.6096 * x - 0.7170 * x**2 + 0.1347 * x**3

        elif self.model == 'curti':
            O3O2 = -0.691 - 2.944 * x - 1.308 * x**2
            R23  = 0.527 - 1.569 * x - 1.652 * x**2 - 0.421 * x**3
            N2   = -0.489 + 1.513 * x - 2.554 * x**2 - 5.293 * x**2 - 2.867 * x**4

        OII = 1. / 10**O3O2
        HB = (4./3. + OII) / 10**R23

        if self.use_redlines & (self.model == 'maiolino'):
            #NII = 1. / 10**O3N2
            #Ha = NII / 10**N2
            NII = 10**N2

        if self.use_redlines & (self.model == 'curti'):
            #Ha = (2.86*HB)/(10**(-0.4 * (self.calzettilaw(np.array([6562.]))-self.calzettilaw(np.array([4861.]))) *thetas[1]))
            #NII = (10**N2)*Ha
            NII = 10**N2
        
        if self.use_NeIII:
            Ne3O2 = -1.2608 - 1.0861 * x - 0.147 * x**2
            Ne3 = 10**Ne3O2 * OII

            if self.use_redlines: 
                model_y_b = np.array([OII, Ne3, HB, 1.]) * thetas[2]
                model_y_r = np.array([1., NII]) * thetas[3]
                model_y = np.hstack((model_y_b, model_y_r)) * 10**(-0.4 * Abv)
            else:
                model_y = np.array([OII, Ne3, HB, 1.]) * 10**(-0.4 * Abv) * thetas[2]
        else:
            if self.use_redlines:
                model_y_b = np.array([OII, HB, 1.]) * thetas[2]
                model_y_r = np.array([1., NII]) * thetas[3]
                model_y = np.hstack((model_y_b, model_y_r)) * 10**(-0.4 * Abv)
            else:
                model_y = np.array([OII, HB, 1.]) * 10**(-0.4 * Abv) * thetas[2]

        return model_y

    def lnprior(self, thetas):
        ''' Uniform prior '''
        flag = True
        for theta, lim in zip(thetas, self.lims):
            flag *= (theta > lim[0])
            flag *= (theta < lim[1])
        if flag:
            e_bv = thetas[1]
            mu = 0.295
            si = 0.165
            return np.log(1.0/(np.sqrt(2*np.pi)*si))*(-0.5*((e_bv-mu)**2)/(si**2)) #guassian

            # #guassian E(B-V) prior for de-reddened fluxes
            # mu = 0.0
            # si = 0.1
            # return np.log(1.0/(np.sqrt(2*np.pi)*si))*(-0.5*((e_bv-mu)**2)/(si**2)) #guassian
        else:
            return -np.inf

    def lnlike_ind(self, data, error, thetas):
        model_y = self.get_model_fluxes(thetas)
        #model_e = self.get_model_err(thetas)
        inv_sigma2 = 1.0 / (error**2 + (model_y * self.sigma_m)**2)
        chi2_term = -0.5 * (data - model_y)**2 * inv_sigma2
        parm_term = -0.5 * np.log(1 / inv_sigma2) #what is this term?
        return chi2_term + parm_term

    def lnlike(self, thetas):
        model_y = self.get_model_fluxes(thetas)
        #model_e = self.get_model_err(thetas)
        inv_sigma2 = 1.0 / (self.error**2 + (model_y * self.sigma_m)**2)
        chi2_term = -0.5 * np.sum((self.data - model_y)**2 * inv_sigma2)
        parm_term = -0.5 * np.sum(np.log(1 / inv_sigma2)) #what is this term?
        return chi2_term + parm_term

    def lnprob(self, thetas):
        ''' Combine prior and likelihood handling the out of bounds case '''
        init_term = self.lnprior(thetas)
        if np.isfinite(init_term):
            return init_term + self.lnlike(thetas)
        else:
            return -np.inf

    def fit(self):
        ''' Fit using emcee and a Gaussian ball for initial position '''
        if self.ndim == 3:
            pos = emcee.utils.sample_ball([8.5, 0.15, self.data[-1]],
                                      [0.2, 0.05, self.error[-1]],
                                      size=self.nwalkers)
        elif self.ndim == 4:
            pos = emcee.utils.sample_ball([8.5, 0.15, self.data[-2], self.data[-1]],
                                      [0.2, 0.05, self.error[-2], self.data[-1]],
                                      size=self.nwalkers)

        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, a=2.0)

        # # Do real run
        # sampler.run_mcmc(pos, self.nsteps, rstate0=np.random.get_state())
        # #tau = np.max(sampler.acor)
        # tau = 90
        # burnin_step = int(tau*3)
        # self.samples = sampler.chain[:, burnin_step:, :].reshape((-1, 3))

        print "Burning in ..."
        pos, prob, state = sampler.run_mcmc(pos, self.nsteps[0], rstate0=np.random.get_state())

        # Reset the chain to remove the burn-in samples.
        sampler.reset()

        # Starting from the final position in the burn-in chain, sample for 1000
        # steps. (rstate0 is the state of the internal random number generator)
        print "Running MCMC ..."
        pos, prob, state = sampler.run_mcmc(pos, self.nsteps[1], rstate0=state)

        self.flatchain = sampler.flatchain
        self.samples = sampler.chain[:, :, :].reshape((-1, self.ndim))

        self.results = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                        zip(*np.percentile(self.samples, [16, 50, 84], axis=0)))
        #print self.samples, self.results

    def plot_corner_results(self, outpath, show=False, save=True):
        ''' Plot results '''
        if self.ndim == 3: 
            names = ['Log (12 + O/H)', 'E(B-V)', 'OIII intrinsic']
        elif self.ndim == 4:
            names = ['Log (12 + O/H)', 'E(B-V)', 'OIII intrinsic', r'H$\alpha$ intrinsic']   

        plt.rcParams['xtick.labelsize']=20
        plt.rcParams['ytick.labelsize']=20         

        fig = corner.corner(self.samples, labels=names, range = [(7.4,9.0), (0.0, 1.0), (0,150)],
                            truths=[r[0] for r in self.results], truth_color='lightcoral',
                            label_kwargs={"fontsize": 20, "fontweight":'bold'}, show_titles=True,
                            title_kwargs={"fontsize": 20, "fontweight":'bold'},
                            quantiles=[0.16, 0.5, 0.84], bins=30)
        #plt.tick_params(axis='both', labelsize=25)
        fig.set_size_inches(13,9)

        if save:
            plt.savefig(op.join(outpath, self.obj_name+'_cornerplot.pdf'), bbox_inches='tight')
        if show:
            print(op.join(outpath,'met_cornerplots', self.obj_name+'_cornerplot.pdf'))
            plt.show()

    def plot_best_solution(self, outpath, show=False, save=True):
        flatchain = self.flatchain
        x_samples = flatchain[:,0]
        E_bv_samples = flatchain[:,1]

        mc_x    = self.results[0]
        mc_E_bv = self.results[1]

        E_sig = (mc_E_bv[1]+mc_E_bv[2])/2
        x_sig = (mc_x[1]+mc_x[2])/2

        #make the confidience interval plots seen in the maiolino paper
        cov = np.cov(E_bv_samples, x_samples)
        lambda_, v = np.linalg.eig(cov)
        lambda_ = np.sqrt(lambda_)

        # ax = plt.subplot(111)
        # for j in xrange(1, 4):
        #     ell = Ellipse(xy=(mc_E_bv[0], mc_x[0]),
        #                     width=lambda_[0]*j*2, height=lambda_[1]*j*2,
        #                     angle=np.rad2deg(np.arccos(v[0, 0])))
        #     ell.set_facecolor('none')
        #     ell.set_edgecolor('red')
        #     ax.add_artist(ell)

        # plt.scatter(mc_E_bv[0], mc_x[0])
        # plt.xlabel('E(B-V)')
        # plt.ylabel('X')
        # plt.xlim(0.0, 0.60)
        # plt.ylim(7.0, 10.0)
        # plt.text(0.5, 7.3, 'X: '+str(round(mc_x[0],2))+'\n'+'E(B-V): '+str(round(mc_E_bv[0],2)), 
        #     horizontalalignment='center', verticalalignment='center', bbox=dict(facecolor='red', alpha=0.3)) 
        # #plt.text(x, y, s, bbox=dict(facecolor='red', alpha=0.5))

        plt.figure(figsize=(20,20))
        plt.rcParams['xtick.labelsize']=25
        plt.rcParams['ytick.labelsize']=25

        sn = 3
        cp = sns.color_palette()
        #sns.set(rc={'figure.figsize':(15.0, 15.0)})
        a=sns.jointplot(x=E_bv_samples, y=x_samples, kind="kde", color=cp[0], xlim=(mc_E_bv[0]+sn*E_sig, mc_E_bv[0]-sn*E_sig), ylim=(mc_x[0]+sn*x_sig, mc_x[0]-sn*x_sig))
        print 'VALUES: ', mc_E_bv[0], mc_x[0]
        a.ax_joint.scatter(mc_E_bv[0], mc_x[0],color='white', marker='x', s=100)
        a.set_axis_labels('E(B-V)', 'Log (12 + O/H)', fontsize=30)
        #plt.tick_params(axis='both', which='major', labelsize=25)
        plt.tight_layout()
        #a.ax_joint.legend_.remove()
        #plt.scatter(mc_E_bv[0], mc_x[0],color='red', marker='x', s=100)

        if show:
            plt.show()
        if save:
            plt.savefig(op.join(outpath, self.obj_name+'_best_result.pdf'))

    # def plot_result_ratios(self, outpath, show=False, save=True):
    #     mc_x    = self.results[0]
    #     mc_E_bv = self.results[1]

    #     met = np.arange(7.5, 9.3, 0.1)
    #     met_norm = np.subtract(met, self.solar_x)
    #     model_mc_x = mc_x[0] - self.solar_x

    #     if self.use_NeIII:
    #         if self.use_redlines:
    #             OII, NeIII, Hb, OIII, Ha, NII = self.data
    #             OII_e, NeIII_e, Hb_e, OIII_e, Ha_e, NII_e = self.error
    #         else: 
    #             OII, NeIII, Hb, OIII = self.data
    #             OII_e, NeIII_e, Hb_e, OIII_e = self.error
    #     else:
    #         if self.use_redlines:
    #             OII, Hb, OIII, Ha, NII = self.data
    #             OII_e, Hb_e, OIII_e, Ha_e, NII_e = self.error
    #         else:
    #             OII, Hb, OIII = self.data
    #             OII_e, Hb_e, OIII_e = self.error

    #     x_samples = self.flatchain[:,0]
    #     E_bv_samples = self.flatchain[:,1]

    #     c_model_dict = {'R23':mrf.R23_model_c, 'O3O2':mrf.RO32_model_c, 'Ne3O2':mrf.RNeO2_model}

    #     ratio_dict = {'R23':[mrf.R23_model, 
    #                         mrf.R23_ratio(OIII, OII, Hb, E_bv_samples), 
    #                         mrf.R23_ratio(OIII, OII, Hb, mc_E_bv[0]), 
    #                         mrf.R23_ratio(OIII, OII, Hb, 0), 
    #                         mrf.R23_ratio_err(OIII, OII, Hb, OIII_e, OII_e, Hb_e, 0)], 
    #                 'O3O2':[mrf.RO32_model, 
    #                         mrf.RO32_ratio(OIII, OII, E_bv_samples),
    #                         mrf.RO32_ratio(OIII, OII, mc_E_bv[0]),
    #                         mrf.RO32_ratio(OIII, OII, 0), 
    #                         mrf.RO32_ratio_err(OIII, OII, OIII_e, OII_e, 0)]}

    #     if self.use_NeIII:

    #         ratio_dict['Ne3O2'] = [mrf.RNeO2_model, 
    #                             mrf.RNeO2_ratio(NeIII, OII, E_bv_samples),
    #                             mrf.RNeO2_ratio(NeIII, OII, mc_E_bv[0]),
    #                             mrf.RNeO2_ratio(NeIII, OII, 0), 
    #                             mrf.RNeO2_ratio_err(NeIII, OII, NeIII_e, OII_e, 0)]

    #     f, ax = plt.subplots(1,len(ratio_dict), sharex=True, figsize=(18, 6))
    #     f.suptitle('Results')

    #     for i, ratio in enumerate(ratio_dict):
    #         #plot the model curves with disperson 
    #         if self.model == 'maiolino':
    #             model_log_ratio = np.log10(ratio_dict[ratio][0](met_norm))
    #         elif self.model == 'curti':
    #             model_log_ratio = np.log10(c_model_dict[ratio](met_norm))
    #         model_log_disp  = self.sigma_m

    #         ax[i].plot(met, model_log_ratio, color='black')
    #         ax[i].plot(met, model_log_ratio+model_log_disp, ls=':', color='black')
    #         ax[i].plot(met, model_log_ratio-model_log_disp, ls=':', color='black')

    #         #plot the values derived from the model and the observed values
    #         mc_model_ratio     = np.log10(ratio_dict[ratio][0](model_mc_x))
    #         mc_obs_ratio       = np.log10(ratio_dict[ratio][3])
    #         mc_obs_ratio_deRed = np.log10(ratio_dict[ratio][2])
    #         mc_obs_ratio_err   = ratio_dict[ratio][4]
    #         #mc_log_ratio_err   = (np.log10((10**mc_obs_ratio)+mc_obs_ratio_err)-np.log10((10**mc_obs_ratio)-mc_obs_ratio_err))/2
    #         mc_log_ratio_err   = (1/np.log(10))*(mc_obs_ratio_err/(10**mc_obs_ratio))

    #         ax[i].errorbar(mc_x[0], mc_obs_ratio, xerr=(mc_x[1]+mc_x[2])/2, yerr= mc_log_ratio_err, capsize=7, capthick=2, color='green')
    #         ax[i].scatter(mc_x[0], mc_obs_ratio_deRed , marker='x', s=100, color='blue')

    #         #compute the covariance matrix of the ratios and the metallicties 
    #         #this is based on the values for X and E_bv
    #         ratio_samples = np.log10(ratio_dict[ratio][1])

    #         cov = np.cov(x_samples, ratio_samples)
    #         lambda_, v = np.linalg.eig(cov)
    #         lambda_ = np.sqrt(lambda_)

    #         ell = Ellipse(xy=(mc_x[0], mc_obs_ratio_deRed),width=lambda_[0]*2, height=lambda_[1]*2, angle=np.rad2deg(np.arccos(v[0, 0])))
    #         ell.set_facecolor('none')
    #         ell.set_edgecolor('red')
    #         ax[i].add_artist(ell)

    #         #plot the axis labels
    #         ax[i].set_xlabel('12+log(O/H)')
    #         ax[i].set_ylabel('log '+ ratio)

    #     if show:
    #         plt.show()
    #     if save:
    #         plt.savefig(op.join(outpath, self.obj_name+'_result_ratios.png'))
