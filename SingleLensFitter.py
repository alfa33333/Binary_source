import sys
import numpy as np
from scipy import linalg
from scipy import optimize
import emcee
# import george
# from george import kernels

import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from pylab import subplots_adjust

import corner
# from astropy.stats import mad_std

from scipy.interpolate import interp1d
# import nestle

class SingleLensFitter():

	"""Class definition for single microlens event fitter."""

	def __init__(self,data,initial_parameters,eigen_lightcurves=None,reference_source=None):

		"""Initialise the SingleLensFitter.

		inputs:

			data:           A dictionary with each key being a data set name string and each 
						value being a tuple of (date, flux, flux_err) each being numpy
						arrays.

			initial_parameters: A numpy array of starting guess values for u_0, t_0, t_E.

			eigen_lightcurves:      If defined, this should be a dictionary with the same keys as data,
						and each value being an n x m numpy array of n lightcurves each 
						with m data values corresponding to the dates of the input data.
						These lightcurves (detrend vectors) are fitted linearly to the data
						at the same time as the magnification and flux parameters. An
						eigen_lightucrves dictionary entry can be defined for any subset 
						of data sources, and different data sources can have different
						numbers of eigenlightcurves.

		"""

		self.data = data
		self.eigen_lightcurves = eigen_lightcurves
		self.initial_parameters = initial_parameters
		self.p = initial_parameters

		if reference_source is None:

			self.reference_source = list(self.data.keys())[0]

			print('Using',self.reference_source,'as reference.')

		else:

			if reference_source in self.data:

				self.reference_source = reference_source
				print('Using',self.reference_source,'as reference.')

			else:

				self.reference_source = list(self.data.keys())[0]

				print('Warning:',reference_source,'is not a valid data source.')
				print('Using',self.reference_source,'as reference.')


		self.parameter_labels = [r"$u_0$",r"$t_0$",r"$t_E$"]
		self.ndim = 3
		self.zp = 25 # Basically never used unless renormalization.

		self.use_finite_source = False

		self.marginalise_linear_parameters = True
		self.fit_blended = True

		self.u0_limits = (0.0,1.3)
		self.tE_limits = (0.5,200)
		self.t0_limits = None
		self.lrho_limits = (-6,0)

		self.use_gaussian_process_model = False
		self.GP_default_params = (1.0,-2.0)
		self.ln_a_limits = (-5,15)
		self.ln_tau_limits = (-5.5, 7.0)
		self.ln_a = {}
		self.ln_tau = {}

		self.use_mixture_model = False
		self.mixture_default_params = (0.0001,1.0e8,0.0)
		self.P_b_limits = (0.0,0.05)
		self.V_b_limits = (0.0,1.0e12)
		self.Y_b_limits = (-1.e5,1.e5)
		self.P_b = {}
		self.V_b = {}
		self.Y_b = {}

		self.use_finite_source = False

		self.nwalkers = 50
		self.nsteps = 50
		self.nsteps_production = 500
		self.thresh_std = 0.1
		self.thresh_mean = 0.05
		self.max_burnin_iterations = 200
		self.emcee_optimize_first = True
		self.emcee_lnp_convergence_threshold = 0.5
		self.emcee_mean_convergence_threshold = 0.1
		self.emcee_std_convergence_threshold = 0.1

		self.plotprefix = 'single_lens_fit'
		self.samples = None
		self.n_plot_samples = 30
		self.make_plots = True
		self.plotrange = None

		self.plot_colours = ['#FF0000', '#000000', '#008000', '#800080', '#FFA500', '#A52A2A', '#ff9999', '#999999', \
					'#4c804c', '#804d80', '#ffdb98', '#a56262', '#CCFF66', '#CC9900', '#9966FF', '#FF3366']

		#Nested parameters.
		self.nlive = 1000 # livepoints.
		self.bound = 'multi' #bounding options = 'single', 'multi'
		self.tol = 0.01 #stopping criterion for nested sampling  dlogz
		
		self.t0t1extlim = False
		self.t02_limits = None
		self.sampletype = 'rwalk'

		#Parallel parameters
		self.cpu = None
		self.dynestyparallel = False

		#Binary source parameters
		self.use_binary_source = False
		self.binary_both = False
		self.binsource_default_paramas =(0.5,np.pi/2.,0.8) #parameters d, theta, q
		self.binboth_default_paramas =(self.p[0]+ 0.1,self.p[1]+1.,0.8) #parameters u02, t02, q
		self.dist_limits = (0.0, 10000.0)
		self.thetabin_limits = (0.0, 2*np.pi)
		self.qlum_limits = (0.0,1.0)

		#optimization options
		self.gen_optimize_first = False

		return

		
	#------------------------------------------------
	#Extra utilities
	def fcreate(self,filename):
		with open(self.plotprefix+filename,'w') as fid:
			fid.write('\n')
		fid.close()
		return None

	def fwrite(self,parameter_labels,parameters):
		with open(self.plotprefix+'.fit_results','a+') as fid:
			for i in range(len(parameter_labels)):
				fid.write('%s %f %f %f\n'%(parameter_labels[i],
							   parameters[i][0],parameters[i][1],parameters[i][2]))

		fid.close()
		return None

	def fappend(self,label,data):
		''' label  = A message or a series of labels
			data   = A 2D numpy array to output '''
		with open(self.plotprefix+'.fit_results','a+') as fid:
			if not isinstance(data , np.ndarray):
				#This is just for messages
				fid.write('%s \n' % label )
			else:
				if len(data.shape) == 2:
					for i in range(len(label)):
						fid.write('%s '%(label[i]))
						for j in range(data.shape[1]):
							fid.write(' %.3f'%data[i,j])
						fid.write('\n')
				else:
					for i in range(len(label)):
						fid.write('%s %.3f\n'%(label[i],data[i]))
		fid.close()
		return None


	#------------------------------------------------
	# Extra Chi2 calculation help functions

	def chi2_calc(self, p = None, source=None):
		ind_chi2 = {}
		temp = self.p
		if p is None:
			p = self.p
		else:
			self.p = p

		if source is None:
			chi2 = 0.0
			# this is the case for chi2 for all combined sites

			for data_set_name in list(self.data.keys()):

				t, y, yerr = self.data[data_set_name]
				if self.use_binary_source is True:
					mag = self.mag_binsource(t)
				else:
					mag = self.magnification(t)

			
				lp = self.linear_fit_reduced(data_set_name,mag)
				chi2 += lp
				a, _ = self.linear_fit(data_set_name, mag)
				a0 = a[1]
				a1 = a[0]
				ind_chi2[data_set_name] = (mag*a0+ a1 - y)**2 / yerr**2
		else:
			# This is the case for a particular site
			t, y, yerr = self.data[source]
			if self.use_binary_source is True:
				mag = self.mag_binsource(t)
			else:
				mag = self.magnification(t)
		
			chi2 = self.linear_fit_reduced(source,mag)
			a, _ = self.linear_fit(source, mag)
			a0 = a[1]
			a1 = a[0]
			ind_chi2[source] = (mag*a0+ a1 - y)**2 / yerr**2
		self.p = temp
			
		return chi2, ind_chi2
		

	def linear_fit_reduced(self,data_key,mag):
		''' Only obtains the chi2 per site, no support for
			mixed model'''

		_, y, yerr = self.data[data_key]

		if self.use_gaussian_process_model:
			C = self.cov
		else:
			C = np.diag(yerr**2)

		C_inv = linalg.inv(C)

		if self.fit_blended:
			A = np.vstack((np.ones_like(mag),mag))
			n_params = 2
		else:
			A = (mag-1.0).reshape(1,len(mag))
			n_params = 1

		if self.eigen_lightcurves is not None:
			if data_key in self.eigen_lightcurves:
				eigs = self.eigen_lightcurves[data_key]
				for i in range(eigs.shape[0]):
					A = np.vstack((A,eigs[i,:]))
					n_params += 1

		S = np.dot(A,np.dot(C_inv,A.T)).reshape(n_params,n_params)
		b = np.dot(A, np.dot(C_inv,y).T)

		if self.marginalise_linear_parameters:

			try:
				M = linalg.inv(S)
			except linalg.LinAlgError:
				return (0,0), -np.inf

			g = np.dot(M,b)
			D = y - np.dot(A.T,g)
			chi2 = np.dot(D.T,np.dot(C_inv,D))

			return chi2

		else:

			try:
				a = linalg.solve(S,b)
			except linalg.LinAlgError:
				return (0,0), -np.inf
			D = y - np.dot(A.T,a)
			chi2 = np.dot(D.T,np.dot(C_inv,D))

			return chi2
	#------------------------------------------------
	#Binary source section
	def add_binary_source(self,system=None):
	
		self.use_binary_source =  True
		self.ndim +=3
		
		if system == 'both' :
			self.binary_both = True
			self.p = np.hstack((self.p,self.binboth_default_paramas))
			self.parameter_labels.extend([r"$u_{02}$",r"$t_{02}$",r"$q$"])
			self.dist_limits = self.u0_limits
			t0lim,_ ,_ = self.data[list(self.data.keys())[0]]
			self.thetabin_limits = (t0lim[0],t0lim[-1])
		else:
			self.p = np.hstack((self.p,self.binsource_default_paramas))
			self.parameter_labels.extend([r"$d$",r"$\theta$",r"$q$"])
		
	def mag_binsource(self,t,p=None):
	
		if p is None:
			p = self.p
			
		u0,t0,te,par1,par2,par3 = p[:6]
		
		tau = (t-t0)/te
		
		u1 = np.sqrt(tau**2 + u0**2)
		
		A1 = (u1**2 + 2.0)/(u1*np.sqrt(u1**2+4.0))
		
		if self.binary_both is True:
			tau2 =  (t-par2)/te
			
			u02 = par1
		else:
			tau2 = tau + par1*np.cos(par2) 
			u02 = u0 + par1*np.sin(par2)
		
		u2 = np.sqrt(tau2**2+u02**2)
		
		A2 = (u2**2 + 2.0)/(u2*np.sqrt(u2**2+4.0))
		
		Atot = (A1 + par3*A2)/(1.+par3)
		
		return Atot
	#------------------------------------------------
	def add_finite_source(self,lrho=None):

		from mpmath import ellipe
		
		self.use_finite_source = True
		self.finite_source_index = self.ndim

		if lrho is not None:
			self.lrho = lrho
		else:
			self.lrho = -3.0

		self.p = np.hstack((self.p,self.lrho))
		self.parameter_labels.append(r"$log_{10} rho$")
		self.ndim += 1

		lz = np.arange(-5,10,0.01)
		fsz = 10.0**lz
		ell = np.zeros_like(fsz)
		rsz = np.arcsin(1.0/fsz)
		p2 = np.pi/2.0
		for zi in range(len(lz)):
			if fsz[zi] < 1.0:
				ell[zi] = ellipe(p2,fsz[zi])
			else:
				ell[zi] = ellipe(rsz[zi],fsz[zi])
		self._ellipe_interpolator = interp1d(fsz,ell)

	def add_mixture_model(self):

		self.use_mixture_model = True
		self.mixture_index = self.ndim

		for site in list(self.data.keys()):

			self.p = np.hstack((self.p, self.mixture_default_params))
			self.ndim += 3
			self.parameter_labels.append(site+'_P_b')
			self.parameter_labels.append(site+'_V_b')
			self.parameter_labels.append(site+'_Y_b')



	def add_gaussian_process_model(self):

		self.use_gaussian_process_model = True
		self.gaussian_process_index = self.ndim

		for site in list(self.data.keys()):

			self.p = np.hstack((self.p, self.GP_default_params))
			self.ndim += 2
			self.parameter_labels.append(site+'_ln_a')
			self.parameter_labels.append(site+'_ln_tau')



	def lnprior_GP(self,GP_params):
		params = ['ln_a','ln_tau']
		for i, p in enumerate(params):
			prange = eval('self.'+p+'_limits')
			if prange:
				if GP_params[i] < prange[0] or GP_params[i] > prange[1]:
					return -np.inf
		return 0.0

	

	def lnprior_mixture(self,mixture_params):
		params = ['P_b','V_b','Y_b']
		for i, p in enumerate(params):
			prange = eval('self.'+p+'_limits')
			if prange:
				if mixture_params[i] < prange[0] or mixture_params[i] > prange[1]:
					return -np.inf
		return 0.0

	

	def lnprior_ulens(self):
		params = ['u0','t0','tE']
		for p in params:
			prange = eval('self.'+p+'_limits')
			if prange:
				param = eval('self.'+p)
				if param < prange[0] or param > prange[1]:
					return -np.inf
		return np.log(1./((1.3)*(200-0.5)*(8240-8220)))

	def lnprior_binsource(self):
		params = ['dist','thetabin','qlum']
		for p in params:
			prange = eval('self.'+p+'_limits')
			if prange:
				param = eval('self.'+p)
				if param < prange[0] or param > prange[1]:
					return -np.inf
		return 0.0

	def lnprior_lrho(self,lrho):
		if lrho < self.lrho_limits[0] or lrho > self.lrho_limits[1]:
			return -np.inf
		return 0.0


	def lnprob(self,p):

		self.p = p
		self.u0 = p[0]
		self.t0 = p[1]
		self.tE = p[2]
		
		lp = self.lnprior_ulens()
        
		if self.use_binary_source == True:
			self.dist = p[3]
			self.thetabin = p[4]
			self.qlum = p[5]
			# if  p[4] > 8232 or p[4] < 8228:
			# 	return -np.inf
			lp += self.lnprior_binsource()
		

		if self.use_finite_source:
			lp += self.lnprior_lrho(p[self.finite_source_index])

		if self.use_mixture_model:

			pi = self.mixture_index
			for data_set_name in list(self.data.keys()):

				self.P_b[data_set_name] = p[pi]
				self.V_b[data_set_name] = p[pi+1]
				self.Y_b[data_set_name] = p[pi+2]
				lp += self.lnprior_mixture(p[pi:pi+3])
				pi += 3

		if self.use_gaussian_process_model:

			pi = self.gaussian_process_index
			for data_set_name in list(self.data.keys()):

				self.ln_a[data_set_name] = p[pi]
				self.ln_tau[data_set_name] = p[pi+1]
				lp += self.lnprior_GP(p[pi:pi+2])
				pi += 2

		if np.isfinite(lp):
			lp += self.lnlikelihood()
		else:
			return -np.inf

		return lp

	def neglnprob(self,p):
		return -self.lnprob(p)


	def lnlikelihood(self):

		lnprob = 0.0

		for data_set_name in list(self.data.keys()):

			t, y, yerr = self.data[data_set_name]
			if self.use_binary_source is True:
				mag = self.mag_binsource(t)
			else:
				mag = self.magnification(t)

			if self.use_gaussian_process_model:
				continue
				# a = np.exp(self.ln_a[data_set_name])
				# tau = np.exp(self.ln_tau[data_set_name])
				# gp = george.GP(a * kernels.ExpKernel(tau))
				# gp.compute(t, yerr)
				# self.cov = gp.get_matrix(t)
				# result, lp = self.linear_fit(data_set_name,mag)
				# model = self.compute_lightcurve(data_set_name,t)
				# lnprob = gp.lnlikelihood(y-model)
			else:
				result, lp = self.linear_fit(data_set_name,mag)
				lnprob += lp

		return lnprob


	def linear_fit(self,data_key,mag):

		_, y, yerr = self.data[data_key]

		if self.use_gaussian_process_model:
			C = self.cov
		else:
			C = np.diag(yerr**2)

		C_inv = linalg.inv(C)

		if self.fit_blended:
			A = np.vstack((np.ones_like(mag),mag))
			n_params = 2
		else:
			A = (mag-1.0).reshape(1,len(mag))
			n_params = 1

		if self.eigen_lightcurves is not None:
			if data_key in self.eigen_lightcurves:
				eigs = self.eigen_lightcurves[data_key]
				for i in range(eigs.shape[0]):
					A = np.vstack((A,eigs[i,:]))
					n_params += 1

		S = np.dot(A,np.dot(C_inv,A.T)).reshape(n_params,n_params)
		b = np.dot(A, np.dot(C_inv,y).T)

		if self.marginalise_linear_parameters:

			try:
				M = linalg.inv(S)
			except linalg.LinAlgError:
				return (0,0), -np.inf

			g = np.dot(M,b)
			D = y - np.dot(A.T,g)
			chi2 = np.dot(D.T,np.dot(C_inv,D))



			if self.use_mixture_model:

				detM = linalg.det(M)
				lnprob = np.log( 2*np.pi*np.sum( (1.0 - self.P_b[data_key]) * \
						np.exp(-D**2/(2.0*yerr**2)) / \
						np.sqrt(detM) + \
						self.P_b[data_key]*np.exp(-(y-self.Y_b[data_key])**2 / \
						2*(self.V_b[data_key]+yerr**2)) / \
						np.sqrt(2*np.pi*(self.V_b[data_key]+yerr**2)) ))
			else:

				lnprob = np.log(2*np.pi) - 0.5*chi2 - 0.5*np.log(linalg.det(M))

			if self.gen_optimize_first:
				lnprob /= 0.1

			return g, lnprob

		else:

			try:
				a = linalg.solve(S,b)
			except linalg.LinAlgError:
				return (0,0), -np.inf
			D = y - np.dot(A.T,a)
			chi2 = np.dot(D.T,np.dot(C_inv,D))

			if self.use_mixture_model:

				lnprob = np.sum( np.log( (1.0 - self.P_b[data_key])*np.exp(-D**2/(2.0*yerr**2)) / \
						np.sqrt(2*np.pi*yerr**2) + \
						self.P_b[data_key]*np.exp(-(y-self.Y_b[data_key])**2 / \
						2*(self.V_b[data_key]+yerr**2)) / \
						np.sqrt(2*np.pi*(self.V_b[data_key]+yerr**2)) ))

			else:

				lnprob = -np.log(np.sum(np.sqrt(2*np.pi*yerr**2))) - 0.5*chi2

			return a, lnprob


	def magnification(self,t,p=None):
	
		if self.use_binary_source is True:
		
			return self.mag_binsource(t,p=p)
		
		else:

			if p is None:
				p = self.p

			u0, t0, tE = p[:3]

			if self.use_finite_source:
				lrho = p[self.finite_source_index]

			tau = (t-t0)/tE

			u = np.sqrt(u0**2+tau**2)

			A = (u**2 + 2.0)/(u*np.sqrt(u**2+4.0))

			if self.use_finite_source:

				rho = 10.0**lrho


				#
				#  This is from Lee et al. (2009) ApJ, 695, 200 - there seems to be a bug somewhere
				#

				# n = 10

				# for q in range(len(t)):

				# 	if u[q] <= rho:

				# 		k = np.arange(1,n)
				# 		theta = np.pi*k/n
				# 		u2 = u[q]*np.cos(theta) + np.sqrt(rho**2 - u[q]**2 * np.sin(theta)**2)
				# 		f1 = u2 * np.sqrt(u2**2 + 4.0)

				# 		k = np.arange(1,n+1)
				# 		theta = np.pi*(2.0*k - 1.0)/(2.0*n)
				# 		u2 = u[q]*np.cos(theta) + np.sqrt(rho**2 - u[q]**2 * np.sin(theta)**2)
				# 		f2 = u2 * np.sqrt(u2**2 + 4.0)


				# 		A[q] = (1.0/(2.0*n*rho**2)) * ( ((u[q]+rho)/3.0) * np.sqrt((u[q]+rho)**2+4.0) - \
				# 										((u[q]-rho)/3.0) * np.sqrt((u[q]-rho)**2+4.0) + \
				# 										(2.0/3.0) * np.sum(f1) + (4.0/3.0) * np.sum(f2) )

				# 	else:

				# 		k = np.arange(1,n/2)
				# 		theta = 2.0*k*np.arcsin(rho/u[q])/n
				# 		if u[q] <= np.arcsin(rho/u[q]):
				# 			u1 = u[q] * np.cos(theta) - np.sqrt(rho**2 - u[q]**2 * np.sin(theta)**2)
				# 			u2 = u[q] * np.cos(theta) + np.sqrt(rho**2 - u[q]**2 * np.sin(theta)**2)
				# 			f1 = u2 * np.sqrt(u2**2 + 4.0) - u1 * np.sqrt(u1**2 + 4.0)
				# 		else:
				# 			f1 = 0.0

				# 		k = np.arange(1,n/2+1)
				# 		theta = (2.0*k-1.0)*np.arcsin(rho/u[q])/n
				# 		if u[q] <= np.arcsin(rho/u[q]):
				# 			u1 = u[q] * np.cos(theta) - np.sqrt(rho**2 - u[q]**2 * np.sin(theta)**2)
				# 			u2 = u[q] * np.cos(theta) + np.sqrt(rho**2 - u[q]**2 * np.sin(theta)**2)
				# 			f2 = u2 * np.sqrt(u2**2 + 4.0) - u1 * np.sqrt(u1**2 + 4.0)
				# 		else:
				# 			f2 = 0.0

				# 		A[q] = (np.arcsin(rho/u[q])/(np.pi*n*rho**2)) * \
				# 						(   ((u[q]+rho)/3.0) * np.sqrt((u[q]+rho)**2+4.0) - \
				# 							((u[q]-rho)/3.0) * np.sqrt((u[q]-rho)**2+4.0) + \
				# 							(2.0/3.0) * np.sum(f1) + (4.0/3.0) * np.sum(f2) )



				#
				#  This is from Gould (1994) ApJ, 421, L71
				#

				z = np.abs(u)/rho

				fs_points = np.where(z < 10.0)[0]

				print('fs_points', fs_points)
				print('z[fs_points]', z[fs_points])

				B0 = 4 * z[fs_points] * self._ellipe_interpolator(z[fs_points]) / np.pi
				A[fs_points] *= B0

			return A



	def emcee_has_converged(self,sampler,n_steps=100):

		# Emcee convergence testing is not easy. The method we will adopt
		# is to test whether the parameter means and standard deviations, and ln_p, have 
		# stabilised, comparing the last n steps, with the previous n steps.

		#std_threshold = 0.01
		#mean_threshold = 0.01

		n_test = sampler.chain.shape[0]*n_steps

		lnp = sampler.lnprobability.T.ravel()
		if len(lnp) < 2*n_test:
			return False

		converged = True

		steps = sampler.chain.shape[1]

		with open(self.plotprefix+'_lnp','a') as fid:

			fid.write("After %d steps, parameter means, standard deviations, convergence metrics and ln_P:\n"%steps)

			for k in range(sampler.chain.shape[2]):

				samples = sampler.chain[:,:,k].T.ravel()
				mean2 = np.mean(samples[-2*n_test:-n_test])
				mean1 = np.mean(samples[-n_test:])
				std2 = np.std(samples[-2*n_test:-n_test:])
				std1 = np.std(samples[-n_test:])

				delta_param = np.abs(mean1 - mean2)/std1
				delta_std = np.abs(std2-std1)/std2

				fid.write("%g %g %g %g\n"%(mean1,std1,delta_param,delta_std))

				if  delta_param > self.emcee_mean_convergence_threshold:
					converged = False

				if  delta_std > self.emcee_std_convergence_threshold:
					converged = False

			lnp_delta = np.mean(lnp[-n_test:]) - np.mean(lnp[-2*n_test:-n_test])

			if lnp_delta > self.emcee_lnp_convergence_threshold:
				converged = False

			fid.write("delta lnp: %10.4f %d\n"%(lnp_delta,converged))

		return converged


	def lnprior_nest(self,ucube):

		cube = np.copy(ucube)
		#prior for hypercube for the main parameters, flat priors.
		cube[0]= cube[0]*(self.u0_limits[1] - self.u0_limits[0]) + self.u0_limits[0]
		if self.t0t1extlim:
			cube[1]= cube[1]*(self.t0_limits[1]-self.t0_limits[0]) + self.t0_limits[0]
		else:
			cube[1]= cube[1]*(self.data[list(self.data.keys())[0]][0][-1]-self.data[list(self.data.keys())[0]][0][0])+ self.data[list(self.data.keys())[0]][0][0]
		cube[2]= cube[2]*(self.tE_limits[1] - self.tE_limits[0]) + self.tE_limits[0]
		if self.use_binary_source is True:
			cube[3] = cube[3]*(self.dist_limits[1] - self.dist_limits[0]) + self.dist_limits[0]
			if self.binary_both:
				if self.t0t1extlim and (self.t02_limits is not None) :
					cube[4] = cube[4]*(self.t02_limits[1]-self.t02_limits[0]) + self.t02_limits[0]
				else:
					cube[4] = cube[4]*(self.data[list(self.data.keys())[0]][0][-1]-self.data[list(self.data.keys())[0]][0][0])+ self.data[list(self.data.keys())[0]][0][0]
			else:
				cube[4] = cube[4]*(self.thetabin_limits[1] - self.thetabin_limits[0]) + self.thetabin_limits[0]
			cube[5] = cube[5]*(self.qlum_limits[1] - self.qlum_limits[0]) + self.qlum_limits[0]

		return cube

	def lnlike_nest(self,cube):#Likelihood function to evaluate for nested. 
		lnprob = 0.0

		for data_set_name in list(self.data.keys()):

			t, y, yerr = self.data[data_set_name]
			
			mag = self.magnification(t,p = cube)

			if self.use_gaussian_process_model:
				continue
			# 	a = np.exp(self.ln_a[data_set_name])
			# 	tau = np.exp(self.ln_tau[data_set_name])
			# 	gp = george.GP(a * kernels.ExpKernel(tau))
			# 	gp.compute(t, yerr)
			# 	self.cov = gp.get_matrix(t)
			# 	result, lp = self.linear_fit(data_set_name,mag)
			# 	model = self.compute_lightcurve(data_set_name,t)
			# 	lnprob = gp.lnlikelihood(y-model)

			else:
				result, lp = self.linear_fit(data_set_name,mag)
				lnprob += lp

		return lnprob 

	def Nested(self):# Nested sampling main calling method.
		from pymultinest import solve 

		if self.p is None:
			raise Exception('Error in SingleLensFitter.fit(): No initial_parameters found.')
		#prints the boundary to be used for the prior
		if self.t0t1extlim:
			t0lim = self.t0_limits
		else:
			t0lim,_,_ = self.data[list(self.data.keys())[0]]
		print('Prior boundaries:')
		print(('u_0 = [%f , %f]'% self.u0_limits))
		print(('t_0 = [%f , %f]'% (t0lim[0],t0lim[-1])))
		print(('t_E = [%f , %f]'% self.tE_limits))

		#basic parameters to run the sampler
		#This parameters will be moved to the class.
		nlive = 600 # livepoints.
		ndim = self.ndim #number of dimensions
		tol = 0.01 #stopping criterion
		print('')
		#run the algorithm
		result =solve(LogLikelihood=self.lnlike_nest,Prior=self.lnprior_nest,
					n_dims=ndim,outputfiles_basename=self.plotprefix,
					n_live_points=nlive,
					evidence_tolerance=tol,
					verbose=True)
		u0line = result['samples'].T
		print('')
		print(('evidence:%(logZ).f +- %(logZerr).1f' % result ))
		print('')
		print('parameter values:')
		for name, col in zip(self.parameter_labels,u0line):
			print(('%15s : %.3f +- %.3f' %(name,col.mean(),col.std())))

		return

	def dnesty(self):#Nested sampling using Dynesty
		import dynesty
		from dynesty import NestedSampler
		from dynesty.utils import resample_equal
		from multiprocessing import Pool, cpu_count


		if self.p is None:
			raise Exception('Error in SingleLensFitter.fit(): No initial_parameters found.')
		#prints the boundary to be used for the prior
		if self.t0t1extlim:
			t0lim = self.t0_limits
		else:
			t0lim,_,_ = self.data[list(self.data.keys())[0]]
		print('Prior boundaries:')
		print(('u_0 = [%f , %f]'% self.u0_limits))
		print(('t_0 = [%f , %f]'% (t0lim[0],t0lim[-1])))
		print(('t_E = [%f , %f]'% self.tE_limits))
		if self.use_binary_source is True:
			if self.binary_both:
				print(('u_02 = [%f , %f]'% self.dist_limits))
				if self.t0t1extlim and (self.t02_limits is not None):
					t02lim = self.t02_limits
				else:
					t02lim,_,_ = self.data[list(self.data.keys())[0]]
				print(('t_02 = [%f , %f]'% (t02lim[0],t02lim[-1])))
				print(('q = [%f , %f]'% self.qlum_limits))
			else:
				print(('d = [%f , %f]'% self.dist_limits))
				print(('theta = [%f , %f]'% self.thetabin_limits))
				print(('q = [%f , %f]'% self.qlum_limits))
		#basic parameters to run the sampler
		print('')
		print(('dynesty version: {}'.format(dynesty.__version__)))
		#This parameters will be moved to the class.
		nlive = self.nlive # livepoints.
		bound = self.bound #bounding options
		ndim = self.ndim #number of dimensions
		tol = self.tol #stopping criterion
		sample = self.sampletype #"rwalk" #'rwalk' #sample type
		print('')
		print(('nlive: %i' % nlive))
		print(('Bounding method: %5s' % bound))
		print(('Number of dimensions: %i' % ndim))
		print(('Sampling method: {:s}'.format(sample)))
		print(('stopping tolerance dlogz: %.3f' % tol))
		print('')

		#This will use a parallel sampler if selected
		#Create sampler
		if self.dynestyparallel:
			print("Parallel section")
			cpunum = cpu_count() if self.cpu is None else (self.cpu + 1)
			with Pool(cpunum-1) as executor:
				sampler = NestedSampler(self.lnlike_nest,self.lnprior_nest,ndim,
					bound=bound,sample=sample,nlive=nlive, \
						pool=executor,
						queue_size = cpunum,
						bootstrap=0)

				sampler.run_nested(dlogz=tol, print_progress=True)#logl_max = -1e-10,dlogz = 1e-20)
		else:
			print("Single thread")
			sampler = NestedSampler(self.lnlike_nest,self.lnprior_nest,ndim,
					bound=bound,sample=sample,nlive=nlive,walks=10)
			#run sampler
			sampler.run_nested(dlogz=tol, print_progress=True) # Output the progress bar

		res = sampler.results # get results dictionary from sampler

		logZdynesty = res.logz[-1]        # value of logZ
		logZerrdynesty = res.logzerr[-1]  # estimate of the statistcal uncertainty on logZ

		# draw posterior samples
		weights = np.exp(res['logwt'] - res['logz'][-1])
		samples_dynesty = resample_equal(res.samples, weights)
		logl_equal = resample_equal(res.logl,weights)
		self.samples = samples_dynesty
		ML = logl_equal[np.where(logl_equal == np.max(logl_equal))[0]]

		#Saving
		bufferarray = np.zeros((res.samples.shape[0],
						  res.samples.shape[1]+2))
		bufferarray[:,0:res.samples.shape[1]] = res.samples
		bufferarray[:,res.samples.shape[1]] = res.logwt
		bufferarray[:,res.samples.shape[1]+1] = res.logl
		np.save(self.plotprefix+'-vol',res.logvol)

		#saving complete sample data
		np.save(self.plotprefix+'-Samples',bufferarray)

		#saving equally  weighted samples

		np.savetxt(self.plotprefix+'-output_equal.txt',samples_dynesty)
		np.savetxt(self.plotprefix+'-evidence.txt',np.array([logZdynesty, logZerrdynesty]))

		#print results

		output = (logZdynesty, logZerrdynesty)
		print()
		print(('evidence: %.3f +- %.3f' % (output[0],output[1]) ))
		print()
		print((res.summary()))
		print()

		if self.use_binary_source is True:
			params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(samples_dynesty[:,:7], \
								[16, 50, 84], axis=0))]
		else:
			params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(samples_dynesty[:,:4], \
								[16, 50, 84], axis=0))]

		self.p = np.asarray(params)[:,0]
		print(self.p)
		self.u0 = self.p[0]
		self.t0 = self.p[1]
		self.tE = self.p[2]
		if self.use_binary_source is True:
			self.dist = self.p[3]
			self.thetabin = self.p[4]
			self.qlum = self.p[5]
        	#rewriting the initial parameters  for plotting
		if self.use_binary_source is True:
			self.initial_parameters = [self.u0,self.t0,self.tE,self.dist,self.thetabin,self.qlum]
		else:
			self.initial_parameters = [self.u0,self.t0,self.tE]

		if self.use_finite_source:
			self.rho = self.p[self.finite_source_index]


		if self.make_plots:

			self.plot_combined_lightcurves()
			self.plot_lightcurves()
			self.plot_chain_corner()


		#printing results
		print('Parameter values:')
		self.fcreate('.fit_results')
		self.fappend('Parameter values:',None)
		self.fwrite(self.parameter_labels,params)
		self.fappend('Maximum likelihood parameters:',None)
		self.fappend(np.array(['']),ML)
		if self.use_finite_source:
			pi = self.finite_source_index
			pitemp = np.asarray([params[pi][0],params[pi][1],params[pi][2]])
			self.fappend(np.array(['rho: ']),pitemp)
		for i in range(len(self.parameter_labels)):
			print(('%10s : %.3f' % (self.parameter_labels[i],self.p[i])))
			if self.use_finite_source:
				print('rho', params[self.finite_source_index])
		self.fappend(np.asarray(['LogZ:','Logzerr:']),np.asarray([res.logz,res.logzerr]))
		#print('parameter values:')
		#for name, col in zip(parameters, samples_dynesty.transpose()):
		#	print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))
		#print(('%15s : %.3f +- %.3f' % (self.parameter_labels[0],np.mean(samples_dynesty[:,0]),np.std(samples_dynesty[:,0])) ))
		#print(('%15s : %.3f +- %.3f' % (self.parameter_labels[1],np.mean(samples_dynesty[:,1]),np.std(samples_dynesty[:,1])) ))
		#print(('%15s : %.3f +- %.3f' % (self.parameter_labels[2],np.mean(samples_dynesty[:,2]),np.std(samples_dynesty[:,2])) ))



		return

	# def nestling(self):
	# 	if (self.add_finite_source or self.add_mixture_model or \
	# 	self.add_mixture_model) is True:
	# 		raise Exception('Error in SingleLensFitter.nestling(): feature is not yet implemented.')
	# 		return None

	# 	if self.p is None:
	# 		raise Exception('Error in SingleLensFitter.fit(): No initial_parameters found.')
	# 		return None
	# 	#prints the boundary to be used for the prior
	# 	t0lim ,_ ,_ = self.data[list(self.data.keys())[0]]
	# 	print('Prior boundaries:')
	# 	print(('u_0 = [%f , %f]'% self.u0_limits))
	# 	print(('t_0 = [%f , %f]'% (t0lim[0],t0lim[-1])))
	# 	print(('t_E = [%f , %f]'% self.tE_limits))
	# 	if self.use_binary_source is True:
	# 		if self.binary_both:
	# 			print(('u_02 = [%f , %f]'% self.dist_limits))
	# 			print(('t_02 = [%f , %f]'% self.thetabin_limits))
	# 			print(('q = [%f , %f]'% self.qlum_limits))
	# 		else:
	# 			print(('d = [%f , %f]'% self.dist_limits))
	# 			print(('theta = [%f , %f]'% self.thetabin_limits))
	# 			print(('q = [%f , %f]'% self.qlum_limits))
	# 	#basic parameters to run the sampler
	# 	print('')
	# 	print(('Nestle version: {}'.format(nestle.__version__)))
	# 	#This parameters will be moved to the class.
	# 	nlive = self.nlive # livepoints.
	# 	bound = self.bound #bounding options
	# 	ndim = self.ndim #number of dimensions
	# 	tol = self.tol #stopping criterion

	# 	print('')
	# 	print(('nlive: %i' % nlive))
	# 	print(('Bounding method: %5s' % bound))
	# 	print(('Number of dimensions: %i' % ndim))
	# 	print(('stopping tolerance dlogz: %.3f' % tol))
	# 	print('')

	# 	# Run nested sampling.
	# 	#rstate = np.random.RandomState(1)
	# 	result = nestle.sample(self.lnlike_nest, self.lnprior_nest,
	# 			ndim, npoints=nlive,method= bound,
	# 			callback = nestle.print_progress)
	# 	print('')
	# 	print('Results summary:')
	# 	print((result.summary()))

	# 	# re-scale weights to have a maximum of one
	# 	nweights = result.weights/np.max(result.weights)

	# 	# get the probability of keeping a sample from the weights
	# 	keepidx = np.where(np.random.rand(len(nweights)) < nweights)[0]

	# 	# get the posterior samples
	# 	samples_nestle = result.samples[keepidx,:]
	# 	equallike_nestle = result.logl[keepidx]
	# 	self.samples = samples_nestle
	# 	ML = samples_nestle[np.where(equallike_nestle == \
	# 						   np.max(equallike_nestle))[0]]
	# 	#Saving
	# 	bufferarray = np.zeros((result.samples.shape[0],
	# 					  result.samples.shape[1]+2))
	# 	bufferarray[:,0:result.samples.shape[1]] = result.samples
	# 	bufferarray[:,result.samples.shape[1]] = result.weights
	# 	bufferarray[:,result.samples.shape[1]+1] = result.logl
	# 	np.save(self.plotprefix+'-vol',result.logvol)

	# 		#saving complete sample data
	# 	np.save(self.plotprefix+'-Samples',bufferarray)

	# 		#saving equally  weighted samples

	# 	bufferarray = np.zeros((samples_nestle.shape[0],
	# 					 samples_nestle.shape[1]+1))
	# 	bufferarray[:,0:samples_nestle.shape[1]] = samples_nestle
	# 	bufferarray[:,samples_nestle.shape[1]] = equallike_nestle
	# 	np.savetxt(self.plotprefix+'-output_equal.txt',bufferarray)

	# 	if self.use_binary_source is True:
	# 		params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(samples_nestle[:,:7], \
	# 							[16, 50, 84], axis=0))]
	# 	else:
	# 		params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(samples_nestle[:,:4], \
	# 							[16, 50, 84], axis=0))]

	# 	self.p = np.asarray(params)[:,0]
	# 	print(self.p)
	# 	self.u0 = self.p[0]
	# 	self.t0 = self.p[1]
	# 	self.tE = self.p[2]
	# 	if self.use_binary_source is True:
	# 		self.dist = self.p[3]
	# 		self.thetabin = self.p[4]
	# 		self.qlum = self.p[5]
    #     	#rewriting the initial parameters  for plotting
	# 	if self.use_binary_source is True:
	# 		self.initial_parameters = [self.u0,self.t0,self.tE,self.dist,self.thetabin,self.qlum]
	# 	else:
	# 		self.initial_parameters = [self.u0,self.t0,self.tE]

	# 	if self.use_finite_source:
	# 		self.rho = self.p[self.finite_source_index]


	# 	if self.make_plots:

	# 		self.plot_combined_lightcurves()
	# 		self.plot_lightcurves()
	# 		self.plot_chain_corner()


	# 	#printing results
	# 	print('Parameter values:')
	# 	self.fcreate('.fit_results')
	# 	self.fappend('Parameter values:',None)
	# 	self.fwrite(self.parameter_labels,params)
	# 	self.fappend('Maximum likelihood parameters:',None)
	# 	self.fappend(np.array(['']),ML)
	# 	if self.use_finite_source:
	# 		pi = self.finite_source_index
	# 		pitemp = np.asarray([params[pi][0],params[pi][1],params[pi][2]])
	# 		self.fappend(np.array(['rho: ']),pitemp)
	# 	for i in range(len(self.parameter_labels)):
	# 		print(('%10s : %.3f' % (self.parameter_labels[i],self.p[i])))
	# 		if self.use_finite_source:
	# 			print('rho', params[self.finite_source_index])
	# 	self.fappend(np.asarray(['LogZ:','Logzerr:']),np.asarray([result.logz,result.logzerr]))
	# 	return

	def fitparallel(self, cpu_cores=1, max_opt_iterations=1):
		from multiprocessing import Pool
		
		if self.p is None:
			raise Exception('Error in SingleLensFitter.fit(): No initial_parameters found.')

		print('Initial parameters:', self.p)
		print('ln Prob = ',self.lnprob(self.p))

		ndim = self.ndim

		self.state = [self.p + 1e-8 * np.random.randn(ndim) \
						for i in range(self.nwalkers)]

		if self.emcee_optimize_first:

			print('Optimising...')

		# 	optimize.minimize(self.neglnprob,self.p,method='Nelder-Mead')

		# 	print('Optimized parameters:', self.p)
		# 	print('ln Prob = ',self.lnprob(self.p))
			print(ndim, len(self.p), self.nwalkers)

		
			with Pool(cpu_cores) as pool:

				sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.lnprob,pool=pool)

				print("Running optimize...")

				iteration = 0
				converged = False
				steps = 0

				self.count = 0

				print('ndim, walkers, nsteps, max_iterations:', ndim, self.nwalkers, self.nsteps, max_opt_iterations)

				while not converged and iteration < max_opt_iterations:

					self.state, lnp ,_=sampler.run_mcmc(self.state,self.nsteps)

					iteration += 1
					print('iteration', iteration, 'completed')

					kmax = np.argmax(sampler.flatlnprobability)
					self.p = sampler.flatchain[kmax,:]

					np.save(self.plotprefix+'-state-optimising',np.asarray(self.state))

					if self.make_plots:

						self.plot_chain(sampler,suffix='-optimising.png')

						ind = 3

						if self.use_finite_source:
							ind += 1

						npar = 0
						if self.use_mixture_model:
							npar += 3
						if self.use_gaussian_process_model:
							npar += 2

						if npar > 0:

							for data_set_name in list(self.data.keys()):
							
								self.plot_chain(sampler,index=list(range(ind,npar+ind)),  \
										suffix='-optimising-'+data_set_name+'.png', \
										labels=self.parameter_labels[ind:ind+npar])
								ind += npar

						self.plot_combined_lightcurves(t_range = self.plotrange)

					converged = self.emcee_has_converged(sampler,n_steps=self.nsteps)
					

				np.save(self.plotprefix+'-chain-optimising', sampler.flatchain)
				np.save(self.plotprefix+'-lnp-optimising',sampler.flatlnprobability)
			self.pmin2o = np.asarray(sampler.flatchain[np.argmax(sampler.flatlnprobability)])
			print('Best optimized value: {}'.format(self.pmin2o))
		###########

		self.gen_optimize_first = False

		print(ndim, len(self.p), self.nwalkers)

		# self.state = [self.p + 1e-8 * np.random.randn(ndim) \
		# 				for i in range(self.nwalkers)]
		with Pool(cpu_cores) as pool:

			sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.lnprob,pool=pool, a=1.25)

			print("Running burn-in...")

			iteration = 0
			converged = False
			steps = 0

			self.count = 0

			print('ndim, walkers, nsteps, max_iterations:', ndim, self.nwalkers, self.nsteps, self.max_burnin_iterations)

			while not converged and iteration < self.max_burnin_iterations:

				self.state, lnp ,_=sampler.run_mcmc(self.state,self.nsteps)

				iteration += 1
				print('iteration', iteration, 'completed')

				kmax = np.argmax(sampler.flatlnprobability)
				self.p = sampler.flatchain[kmax,:]

				np.save(self.plotprefix+'-state-burnin',np.asarray(self.state))

				if self.make_plots:

					self.plot_chain(sampler,suffix='-burnin.png')

					ind = 3

					if self.use_finite_source:
						ind += 1

					npar = 0
					if self.use_mixture_model:
						npar += 3
					if self.use_gaussian_process_model:
						npar += 2

					if npar > 0:

						for data_set_name in list(self.data.keys()):
						
							self.plot_chain(sampler,index=list(range(ind,npar+ind)),  \
									suffix='-burnin-'+data_set_name+'.png', \
									labels=self.parameter_labels[ind:ind+npar])
							ind += npar

					self.plot_combined_lightcurves(t_range = self.plotrange)

				converged = self.emcee_has_converged(sampler,n_steps=self.nsteps)
				

			np.save(self.plotprefix+'-chain-burnin', sampler.flatchain)
			np.save(self.plotprefix+'-lnp-burnin',sampler.flatlnprobability)
			print("Running production...")

			sampler.reset()

			self.state, lnp, _ = sampler.run_mcmc(self.state,self.nsteps_production)

		self.samples = sampler.flatchain

		if self.use_binary_source is True:
			params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(self.samples[:,:7], \
								[16, 50, 84], axis=0))]
		else:
			params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(self.samples[:,:4], \
								[16, 50, 84], axis=0))]

		self.p = np.asarray(params)[:,0]

		self.u0 = self.p[0]
		self.t0 = self.p[1]
		self.tE = self.p[2]
		if self.use_binary_source is True:
			self.dist = self.p[3]
			self.thetabin = self.p[4]
			self.qlum = self.p[5]

		if self.use_finite_source:
			self.rho = self.p[self.finite_source_index]

		if self.make_plots:

			self.plot_chain(sampler,suffix='-final.png')
			ind = 3

			if self.use_finite_source:
				ind += 1

			if npar > 0:

				for data_set_name in list(self.data.keys()):
				
					self.plot_chain(sampler,index=list(range(ind,npar+ind)),  \
							suffix='-final-'+data_set_name+'.png', \
							labels=self.parameter_labels[ind:ind+npar])
					ind += npar
			self.plot_combined_lightcurves(t_range = self.plotrange, prefix = self.plotprefix+'final_combined')
			self.plot_lightcurves()
			self.plot_chain_corner()

		print('Results:')
		print('u0', params[0])
		print('t0', params[1])
		print('tE', params[2])
		
		if self.use_binary_source:
			if self.binary_both is True:
				print('u02', params[3])
				print('t02', params[4])
				print('q', params[5])
			else:
				print('d', params[3])
				print('theta', params[4])
				print('q', params[5])

		if self.use_finite_source:
			print('rho', params[self.finite_source_index])

		with open(self.plotprefix+'.fit_results','w') as fid:
			fid.write('u0 %f %f %f\n'%(params[0][0],params[0][1],params[0][2]))
			fid.write('t0 %f %f %f\n'%(params[1][0],params[1][1],params[1][2]))
			fid.write('tE %f %f %f\n'%(params[2][0],params[2][1],params[2][2]))
			if self.use_binary_source:
				fid.write('par1 %f %f %f\n'%(params[3][0],params[3][1],params[3][2]))
				fid.write('par2 %f %f %f\n'%(params[4][0],params[4][1],params[4][2]))
				fid.write('q %f %f %f\n'%(params[5][0],params[5][1],params[5][2]))
				
			if self.use_finite_source:
				pi = self.finite_source_index
				fid.write('rho %f %f %f\n'%(params[pi][0],params[pi][1],params[pi][2]))


		np.save(self.plotprefix+'-chain-production', sampler.flatchain)
		np.save(self.plotprefix+'-lnp-production',sampler.flatlnprobability)
		np.save(self.plotprefix+'-state-production',np.asarray(self.state))
		np.save(self.plotprefix+'-min_chi2-production',np.asarray(sampler.flatchain[np.argmax(sampler.flatlnprobability)]))
		self.pmin2 = np.asarray(sampler.flatchain[np.argmax(sampler.flatlnprobability)])

		return

	def fit(self):
		from multiprocessing import Pool
		
		if self.p is None:
			raise Exception('Error in SingleLensFitter.fit(): No initial_parameters found.')

		print('Initial parameters:', self.p)
		print('ln Prob = ',self.lnprob(self.p))

		ndim = self.ndim

		if self.emcee_optimize_first:

			print('Optimising...')

			optimize.minimize(self.neglnprob,self.p,method='Nelder-Mead')

			print('Optimized parameters:', self.p)
			print('ln Prob = ',self.lnprob(self.p))

		print(ndim, len(self.p), self.nwalkers)

		self.state = [self.p + 1e-8 * np.random.randn(ndim) \
						for i in range(self.nwalkers)]

		sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.lnprob)

		print("Running burn-in...")
		

		iteration = 0
		converged = False
		steps = 0

		self.count = 0

		print('ndim, walkers, nsteps, max_iterations:', ndim, self.nwalkers, self.nsteps, self.max_burnin_iterations)

		while not converged and iteration < self.max_burnin_iterations:

			self.state, lnp ,_=sampler.run_mcmc(self.state,self.nsteps)

			iteration += 1
			print('iteration', iteration, 'completed')

			kmax = np.argmax(sampler.flatlnprobability)
			self.p = sampler.flatchain[kmax,:]

			np.save(self.plotprefix+'-state-burnin',np.asarray(self.state))

			if self.make_plots:

				self.plot_chain(sampler,suffix='-burnin.png')

				ind = 3

				if self.use_finite_source:
					ind += 1

				npar = 0
				if self.use_mixture_model:
					npar += 3
				if self.use_gaussian_process_model:
					npar += 2

				if npar > 0:

					for data_set_name in list(self.data.keys()):
					
						self.plot_chain(sampler,index=list(range(ind,npar+ind)),  \
								suffix='-burnin-'+data_set_name+'.png', \
								labels=self.parameter_labels[ind:ind+npar])
						ind += npar

				self.plot_combined_lightcurves(t_range = self.plotrange)

			converged = self.emcee_has_converged(sampler,n_steps=self.nsteps)
			

		np.save(self.plotprefix+'-chain-burnin', sampler.flatchain)
		np.save(self.plotprefix+'-lnp-burnin',sampler.flatlnprobability)
		print("Running production...")

		sampler.reset()

		self.state, lnp, _ = sampler.run_mcmc(self.state,self.nsteps_production)

		self.samples = sampler.flatchain

		if self.use_binary_source is True:
			params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(self.samples[:,:7], \
								[16, 50, 84], axis=0))]
		else:
			params = [(v[1], v[2]-v[1], v[1]-v[0]) for v in zip(*np.percentile(self.samples[:,:4], \
								[16, 50, 84], axis=0))]

		self.p = np.asarray(params)[:,0]

		self.u0 = self.p[0]
		self.t0 = self.p[1]
		self.tE = self.p[2]
		if self.use_binary_source is True:
			self.dist = self.p[3]
			self.thetabin = self.p[4]
			self.qlum = self.p[5]

		if self.use_finite_source:
			self.rho = self.p[self.finite_source_index]

		if self.make_plots:

			self.plot_chain(sampler,suffix='-final.png')
			ind = 3

			if self.use_finite_source:
				ind += 1

			if npar > 0:

				for data_set_name in list(self.data.keys()):
				
					self.plot_chain(sampler,index=list(range(ind,npar+ind)),  \
							suffix='-final-'+data_set_name+'.png', \
							labels=self.parameter_labels[ind:ind+npar])
					ind += npar
			self.plot_combined_lightcurves(prefix = 'final_combined')
			self.plot_lightcurves()
			self.plot_chain_corner()

		print('Results:')
		print('u0', params[0])
		print('t0', params[1])
		print('tE', params[2])
		
		if self.use_binary_source:
			if self.binary_both is True:
				print('u02', params[3])
				print('t02', params[4])
				print('q', params[5])
			else:
				print('d', params[3])
				print('theta', params[4])
				print('q', params[5])

		if self.use_finite_source:
			print('rho', params[self.finite_source_index])

		with open(self.plotprefix+'.fit_results','w') as fid:
			fid.write('u0 %f %f %f\n'%(params[0][0],params[0][1],params[0][2]))
			fid.write('t0 %f %f %f\n'%(params[1][0],params[1][1],params[1][2]))
			fid.write('tE %f %f %f\n'%(params[2][0],params[2][1],params[2][2]))
			if self.use_binary_source:
				fid.write('par1 %f %f %f\n'%(params[3][0],params[3][1],params[3][2]))
				fid.write('par2 %f %f %f\n'%(params[4][0],params[4][1],params[4][2]))
				fid.write('q %f %f %f\n'%(params[5][0],params[5][1],params[5][2]))
				
			if self.use_finite_source:
				pi = self.finite_source_index
				fid.write('rho %f %f %f\n'%(params[pi][0],params[pi][1],params[pi][2]))


		np.save(self.plotprefix+'-chain-production', sampler.flatchain)
		np.save(self.plotprefix+'-lnp-production',sampler.flatlnprobability)
		np.save(self.plotprefix+'-state-production',np.asarray(self.state))
		np.save(self.plotprefix+'-min_chi2-production',np.asarray(sampler.flatchain[np.argmax(sampler.flatlnprobability)]))
		self.pmin2 = np.asarray(sampler.flatchain[np.argmax(sampler.flatlnprobability)])

		return


	def plot_chain(self,s,index=None,plot_lnprob=True,suffix='',labels=None):

		if index is None:
			index = list(range(self.ndim))

		if labels is None:
			labels = self.parameter_labels

		ndim = len(index)

		plt.figure(figsize=(8,11))
		
		subplots_adjust(hspace=0.0001)

		for i in range(ndim):

			if i == 0:
				plt.subplot(ndim+plot_lnprob,1,i+1)
				ax1 = plt.gca()
			else:
				plt.subplot(ndim+plot_lnprob,1,i+1,sharex=ax1)

			plt.plot(s.chain[:,:,index[i]].T, '-', color='k', alpha=0.3)

			if labels:
				plt.ylabel(labels[i])

			ax = plt.gca()

			if i < ndim-1+plot_lnprob:
				plt.setp(ax.get_xticklabels(), visible=False)
				ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
				ax.locator_params(axis='y',nbins=4)

		if plot_lnprob:
			plt.subplot(ndim+plot_lnprob,1,ndim+plot_lnprob,sharex=ax1)
			plt.plot(s.lnprobability.T, '-', color='r', alpha=0.3)
			plt.ylabel(r"$ln P$")
			ax = plt.gca()
			ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
			ax.locator_params(axis='y',nbins=4)

		plt.savefig(self.plotprefix+suffix)
		plt.close()



	def compute_lightcurve(self,data_key, x, params=None):

		t, _, _ = self.data[data_key]
		coeffs, _ = self.linear_fit(data_key,self.magnification(t,params))

		if self.fit_blended:
			fx = coeffs[0]+coeffs[1]*self.magnification(x,params)
		else:
			fx = coeffs[0]*(self.magnification(x,params)-1)

		return fx


	def plot_combined_lightcurves(self,t_range=None,y_range=None, prefix = None):

		plt.figure(figsize=(8,11))
	
		colour = iter(plt.cm.jet(np.linspace(0,1,len(self.data))))
	
		if t_range is None:
			t_min = self.p[1]-2*self.p[2]
			t_max = self.p[1]+2*self.p[2]
		else:
			t_min, t_max = t_range

		n_data = len(self.data)

		gs = gridspec.GridSpec(2,1,height_ratios=(3,1))

		# Main lightcurve plot

		ax0 = plt.subplot(gs[0])

		a0 = {}
		a1 = {}

		for site in list(self.data.keys()):

			t, y, yerr = self.data[site]
			mag = self.magnification(t)
			a, lnprob = self.linear_fit(site, mag)
			a0[site] = a[1]
			a1[site] = a[0]

		for k, site in enumerate(self.data.keys()):

			scaled_dflux = a0[self.reference_source]*((self.data[site][1] - a1[site])/a0[site]) + a1[self.reference_source]
			scaled_dflux_err = a0[self.reference_source]*((self.data[site][1] + self.data[site][2] - a1[site])/a0[site]) + \
								a1[self.reference_source] - scaled_dflux

			data_merge = 25 - 2.5*np.log10(scaled_dflux)
			sigs_merge = np.abs(25 - 2.5*np.log10(scaled_dflux+scaled_dflux_err) - data_merge)

			ax0.errorbar(self.data[site][0], data_merge, sigs_merge, fmt='.', ms=2, mec=self.plot_colours[k], \
				c=self.plot_colours[k], label=site)

		# Plot the model

		t_plot = np.linspace(t_min,t_max,10001)

		A = self.magnification(t_plot)

		ax0.plot(t_plot,25-2.5*np.log10(a0[self.reference_source]*A+a1[self.reference_source]),'b-')
		ax0.invert_yaxis()
		plt.ylabel(r'$I_{'+self.reference_source+r'}$')

		ax0.grid()

		plt.legend()
		plt.xlabel('HJD-2450000')

		if y_range is not None:
			plt.ylim(y_range)

		if t_range is not None:
			plt.xlim(t_range)

		xlim = ax0.get_xlim()

		# Residuals plot

		ax1 = plt.subplot(gs[1],sharex=ax0)

		for k, site in enumerate(self.data.keys()):
	
			A = self.magnification(self.data[site][0])

			scaled_dflux = a0[self.reference_source]*((self.data[site][1] - a1[site])/a0[site]) + a1[self.reference_source]
			scaled_dflux_err = a0[self.reference_source]*((self.data[site][1] + self.data[site][2] - a1[site])/a0[site]) + \
								a1[self.reference_source] - scaled_dflux

			scaled_model = 25 -2.5*np.log10(a0[self.reference_source]*A + a1[self.reference_source])
			data_merge = 25 - 2.5*np.log10(scaled_dflux) 
			sigs_merge = np.abs(25 - 2.5*np.log10(scaled_dflux+scaled_dflux_err) - data_merge)

			data_merge -= scaled_model

			ax1.errorbar(self.data[site][0], data_merge, sigs_merge, fmt='.', ms=2, mec=self.plot_colours[k], \
				c=self.plot_colours[k], label=site)

		ax1.set_xlim(xlim)
		ax1.grid()

		if y_range is not None:
			ymean = np.mean(y_range)
			y_range = (y_range[0]-ymean,y_range[1]-ymean)
			plt.ylim(y_range)

		plt.xlabel('HJD-2450000')
		ax1.invert_yaxis()
		plt.ylabel(r'$\Delta I_{'+self.reference_source+r'}$')

		plt.tight_layout()

		if prefix is None:
			plt.savefig(self.plotprefix+'-combined-lightcurve')
		else:
			plt.savefig(prefix+'-combined-lightcurve')
		plt.close()




	def plot_lightcurves(self):

		plt.figure(figsize=(8,11))
	
		colour = iter(plt.cm.jet(np.linspace(0,1,len(self.data))))

		xmin = self.initial_parameters[1]-2*self.initial_parameters[2]
		xmax = self.initial_parameters[1]+2*self.initial_parameters[2]

		n_data = len(self.data)
		for i, data_set_name in enumerate(self.data.keys()):

			t, y, yerr = self.data[data_set_name]
			#c=next(colour)
			c = 'r'

			if i == 0:
				plt.subplot(n_data,1,i+1)
				ax1 = plt.gca()
			else:
				plt.subplot(n_data,1,i+1,sharex=ax1)

			y_cond = y
			if self.eigen_lightcurves is not None:
				if data_set_name in self.eigen_lightcurves:
					coeffs, _ = self.linear_fit(data_set_name,self.magnification(t))
					ci = 1
					if self.fit_blended:
						ci = 2
					eigs = self.eigen_lightcurves[data_set_name]
					for j in range(eigs.shape[0]):
						y_cond -= coeffs[ci+j]*eigs[j,:]

			plt.errorbar(t, y_cond, yerr=yerr, fmt=".", color=c, capsize=0)
			ax=plt.gca()
			ax.set_xlim(xmin,xmax)
			plt.xlabel(r"$\Delta t (d)$")
			plt.ylabel(data_set_name+r"  $\Delta F$")
	
			x = np.linspace(xmin,xmax, 3000)

			if not(self.use_gaussian_process_model):
				plt.plot(x, self.compute_lightcurve(data_set_name,x),color="k")
				ylim = ax.get_ylim()
		
				# Plot posterior samples.
			for s in self.samples[np.random.randint(len(self.samples), size=self.n_plot_samples)]:

				if self.use_gaussian_process_model:
					continue

				# 	# Set up the GP for this sample.
				# 	a, tau = np.exp(s[3+2*i:3+2*i+2])
				# 	gp = george.GP(a * kernels.ExpKernel(tau))
				# 	gp.compute(t, yerr)
				# 	self.cov = gp.get_matrix(t)
				# 	modelt = self.compute_lightcurve(data_set_name,t,params=s)
				# 	modelx = self.compute_lightcurve(data_set_name,x,params=s)

				# 	# Compute the prediction conditioned on the observations
				# 	# and plot it.
				# 	m = gp.sample_conditional(y - modelt,x) + modelx
				# 	plt.plot(x, m, color="#4682b4", alpha=0.3)

				else:

					plt.plot(x, self.compute_lightcurve(data_set_name,x,params=s), \
							color="#4682b4",alpha=0.3)

			if not(self.use_gaussian_process_model):
				ax.set_ylim(ylim)

		plt.savefig(self.plotprefix+'-lc.png')

		plt.close()

 
	def plot_chain_corner(self):

		if self.use_finite_source:

			figure = corner.corner(self.samples[:,:4],
						labels=[r"$u_0$",r"$t_0$",r"$t_E$",r"$\rho$"],
						quantiles=[0.16, 0.5, 0.84],
						truths=(self.u0,self.t0,self.tE,self.rho),
						show_titles=True, title_args={"fontsize": 12})

		elif self.use_binary_source:
			if self.binary_both:
				figure = corner.corner(self.samples[:,:6],
						labels=[r"$u_0$",r"$t_0$",r"$t_E$",r"$u_{02}$",r"$t_{02}$",r"$q$"],
						quantiles=[0.16, 0.5, 0.84],
						truths=(self.u0,self.t0,self.tE,self.dist,self.thetabin,self.qlum),
						show_titles=True, title_args={"fontsize": 12})
			else:
				figure = corner.corner(self.samples[:,:6],
						labels=[r"$u_0$",r"$t_0$",r"$t_E$",r"$d$",r"$\theta$",r"$q$"],
						quantiles=[0.16, 0.5, 0.84],
						truths=(self.u0,self.t0,self.tE,self.dist,self.thetabin,self.qlum),
						show_titles=True, title_args={"fontsize": 12})
		else:

			figure = corner.corner(self.samples[:,:3],
						labels=[r"$u_0$",r"$t_0$",r"$t_E$"],
						quantiles=[0.16, 0.5, 0.84],
						truths=(self.u0,self.t0,self.tE),
						show_titles=True, title_args={"fontsize": 12})

		figure.savefig(self.plotprefix+'-pdist.png')

