#!/home/ch3246/.conda/envs/dalma-python3-CPU/bin/python3

#########################################################

# This routine fits the inertial mode power spectrum with a Lorentzian using Markov Chain Monte Carlo


#########################################################

from Routines import *
import numpy as np

mkdir_p("FIGURES/MCMC")

#-------------------------------------------------------
# user Inputs
#-------------------------------------------------------


MSHIFT 		= 0					# Which channel to fit ell = m + MSHIFT
Mgrid 		= np.arange(3,16)	# The m grid to fit
PMODEgrid 	= [0,1,2]			# Which PMODE spectrum data should we combine togther

#-------------------------------------------------------
# Load in the data
#-------------------------------------------------------
# Load in the full ell,m power spectrum and average over p
POW = 0
for PMODE in PMODEgrid:
	with np.load('/scratch/ch3246/Private/InertialWaves/PowerSpectra_VORT_[12.84years]_1Spect_HMI5deg_FULL_PMODE%i.npz' % PMODE) as rossby_dict:
		nuGrid = rossby_dict['freqGrid']*1e9
		ell_grid  = rossby_dict['Lgrid']
		mm_grid   = rossby_dict['Mgrid']
		POW       += rossby_dict['POWS'].squeeze()[:,mm_grid>=0]
		ell_grid,mm_grid = np.meshgrid(ell_grid,mm_grid[mm_grid>=0],indexing='ij')		
POW = POW/len(PMODEgrid)

# Extract the channel ell = m + MSHIFT
POWt = []
for ii in range(POW.shape[-1]):
	POWt.append(np.diag(POW[...,ii],-MSHIFT))
POW = np.array(POWt).T
ell_grid = np.diag(ell_grid)[MSHIFT:];mm_grid = np.diag(mm_grid,-MSHIFT)

# Clean up
del POWt


#-------------------------------------------------------
# Perform the fits
res_tot = [];res_tot_mcmc = [];ell_fit = [];amp_factors=[];DAT = []
for mm in Mgrid:
	print(dashLine,'Computing mm=%i' % mm)
	plt.close('all')

	# Guess the dispersion and then extract data in window of +/- 225 nHz either side of guess
	dispersion = -2*453.1/(mm+1)+3*mm

	indl = np.argmin(abs(nuGrid - (dispersion-225)));indu = np.argmin(abs(nuGrid - (dispersion+225)))
	freq = nuGrid[indl:indu]
	data = POW[mm,indl:indu]
	amp_factor = np.amax(data)
	data = data/amp_factor

	# Plot the data
	plt.figure(1)
	plt.plot(freq,data,label = 'data')
	plt.ylim(top=1.5)

	#-------------------------------------------------------
	# For each m define the guess and prior model parameters
	# A 	= 	Mode amplitude
	# B 	= 	Background amplitude
	# G 	= 	Line width
	# x0 	=	Central frequency
	#-------------------------------------------------------
	Aguess = 1; Gguess = 20;x0guess = dispersion;Bguess = 0.1;
	Aprior = [0.01,5];Gprior = [np.diff(nuGrid)[0]*2,200];Bprior = [1e-4,5];


	theta_guess = [];priors = []

	theta_guess.append(np.log(Aguess))
	theta_guess.append(np.log(Gguess))
	theta_guess.append(x0guess)

	priors.append(np.log(Aprior))
	priors.append(np.log(Gprior))
	priors.append([nuGrid[indl]+50,nuGrid[indu]-50])


	theta_guess.append(np.log(Bguess))
	priors.append(np.log(Bprior))

	

	#-------------------------------------------------------
	# Define the model to fit
	#-------------------------------------------------------
	def Model_Lorenz(theta,x):
		if (len(theta)-1)%3 != 0:
			raise Exception('Number of entries for A does match number for B or x0')
		nLor = (len(theta)-1)//3

		sumLor = 0.
		for ii in range(nLor):
			Ap     = np.exp(theta[ii*3]);
			Gammap = np.exp(theta[ii*3 + 1])
			x0p    = theta[ii*3 + 2]
			sumLor += Ap*(Gammap/2)**2 / ( (x-x0p)**2 + (Gammap/2)**2 )

		return sumLor + np.exp(theta[-1])


	# plot the guess profile
	plt.plot(freq,Model_Lorenz(theta_guess,freq),label='initial guess')

	def lnlike(theta,x,y):
		return -np.sum(np.log(Model_Lorenz(theta,x)) + y/Model_Lorenz(theta,x))


	# perform a fit using scipy optizmie and plot as well
	import scipy.optimize as op
	nll = lambda *args: -lnlike(*args)
	result = op.minimize(nll, theta_guess, args=(freq, data))
	res = result["x"]
	plt.plot(subSampleVector(freq,0.1),Model_Lorenz(res,subSampleVector(freq,0.1)),label='Minimized -LnLike result')

	#-------------------------------------------------------
	# Define the log priors
	#-------------------------------------------------------
	def lnprior(theta,priors):
		if (len(theta)-1)%3 != 0:
			raise Exception('Number of entries for A does match number for B or x0')

		nLor = (len(theta)-1)//3

		A      = np.exp(theta[0::3])[:-1] # got to ignore the B
		Gamma  = np.exp(theta[1::3])
		x0     = theta[2::3]
		B      = np.exp(theta[-1])

		nChecks = 0

		for ii in range(nLor):
			Aprior  = np.exp(priors[ii*3]);
			Gprior  = np.exp(priors[ii*3+1]);
			x0prior = priors[ii*3+2];
			if Aprior[0] < A[ii] < Aprior[1] and Gprior[0] < Gamma[ii] < Gprior[1] and x0prior[0] < x0[ii] < x0prior[1]:
				nChecks +=1
	

		Bprior  = np.exp(priors[-1])
		if nChecks == nLor and Bprior[0] < B < Bprior[1]:
			return 0
		return -np.inf

	#-------------------------------------------------------
	# define the log probability
	#-------------------------------------------------------

	def lnprob(theta, x, y,priors):
	    lp = lnprior(theta,priors)
	    if not np.isfinite(lp):
	        return -np.inf
	    return lp + lnlike(theta, x, y)

	#-------------------------------------------------------
	# Define initial walker parameters for emcee and run
	#-------------------------------------------------------
	print("Running MCMC")
	ndim, nwalkers = len(priors), 100
	pos = [theta_guess + 1e-2*np.random.randn(ndim)*np.diff(priors,axis=1).squeeze() for i in range(nwalkers)]

	import emcee
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(freq, data,priors))

	sampler.run_mcmc(pos, 1000,rstate0=np.random.get_state())

	# Plot the walker chain
	fig,ax = plt.subplots(4,1,sharex=True)
	for ii in [0,1,2,3]:
		ax[ii].plot(sampler.chain[:,:,ii].T,'k',linewidth=0.1)


	plt.savefig('FIGURES/MCMC/mm%i_walkers.png' % (mm))

	# Cut out the first 500 (burn in)
	samples = sampler.chain[:, 500:, :].reshape((-1, ndim))

	

	# Convert log params back 
	import copy
	samples_ln = copy.copy(samples)
	for ii in range(samples.shape[-1]):
		if (ii+1)%3 != 0:
			samples[:,ii] = np.exp(samples[:,ii])


	# Compute percentile ranges
	res_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
	                             zip(*np.percentile(samples, [16, 50, 84],
	                                                axis=0)))
	res_mcmc = np.array(list(res_mcmc))
	
	# Get data ready for plotting
	res_mcmc_plot = [];res_plot = []
	for ii in range(len(res_mcmc)):
		if (ii+1)%3 != 0:
			res_mcmc_plot.append(np.log(res_mcmc[ii,0]))
			res_plot.append(np.log(res[ii]))
		else:
			res_mcmc_plot.append(res_mcmc[ii,0])
			res_plot.append(res[ii])
	res_mcmc_plot = np.array(res_mcmc_plot)
	res_plot = np.array(res_plot)

	plt.figure(1)
	plt.plot(subSampleVector(freq,0.1),Model_Lorenz(res_mcmc_plot,subSampleVector(freq,0.1)),label='MCMC result')
	plt.legend()


	plt.savefig('FIGURES/MCMC/mm%i_data_fits.png' % (mm))
	

	plt.figure()
	plt.plot(freq,data/Model_Lorenz(res_mcmc_plot,freq),label='mcmc')
	plt.plot(freq,data/Model_Lorenz(res,freq),label='LnLike')
	plt.legend()

	for ii in range((len(priors)-1)//3):
		plt.axvline(x=res_mcmc[2::3,0][ii],color='r')
	plt.axhline(y=np.mean(data/Model_Lorenz(res_mcmc_plot,freq)),color='g')

	plt.savefig('FIGURES/MCMC/mm%i_Minimization.png' % (mm))
	
	res_tot.append(res)
	res_tot_mcmc.append(res_mcmc)
	ell_fit.append(mm)
	amp_factors.append(amp_factor)
	DAT.append([freq,data,Model_Lorenz(res_mcmc_plot,freq),Model_Lorenz(res,freq)])

	print("Fits for m = %i Complete" % mm)

	



np.savez_compressed('data/RDA5deg_fitting_data_ell=m+%i.npz' % MSHIFT ,\
									fit_results_mcmc = res_tot_mcmc,\
									fit_results = res_tot,\
									ell = ell_fit,\
									percentile = [50,16,84],\
									amp_factor = amp_factors,\
									DAT = DAT)