#!/home/ch3246/.conda/envs/dalma-python3-CPU/bin/python3

#########################################################

# This routine reads in a data cube and computes the inertial mode power spectrum


#########################################################

from Routines import *
import numpy as np

#------------------------------------
# User inputs

# Directory containing the flow data cubes. Also location where spectrum is saved
outdir		= "/scratch/ch3246/Private/InertialWaves"

windowSize 	= 12456 						# How many times steps to be used. if smaller than data cube time length, will cut into smaller cubes and compute spectrum of each
mshift 		= 0								# For plotting, which channel do you want to see ell = m+ mshift
lMax 		= 40							# Maximum ell to compute to
OBS 		= ['UX','UY','VORT','DIV'][2]	# Compute spectrum of which observable
#------------------------------------


for PMODE in [0,1,2]:
	print(dashLine,'Computing for P%i mode' % PMODE) 

	#------------------------------------------------------------
	#-----Load in Data
	#------------------------------------------------------------
	with np.load(outdir + '/Ring_Maps_HMI_5deg_PMODE%i.npz' % PMODE) as DICT:
		Tgrid = DICT['Tgrid']
		Lat_HG = DICT['Lat_HG'];Lon_HG = DICT['Lon_HG'];Lon_SH = DICT['Lon_SH']
		if Lon_HG[-1] != 360.:
			Lon_HG = np.append(Lon_HG,360.)

		if OBS.upper() == 'UY':
			observable = DICT['Uy_eqrot']
		elif OBS.upper() == 'UX':
			observable = DICT['Ux_eqrot']
		elif OBS.upper() == 'VORT':
			observable = DICT['Vorticity']
		elif OBS.upper() == 'DIV':
			observable = DICT['Divergence']


	# Compute the number of sub cubes
	Nwindows = observable.shape[-1]//windowSize 


	# For each sub cube, compute the power spectrum
	POWS = [];Tgrid_c_grid = []
	for wind_ind in range(Nwindows):
		#---------------------------------------------------
		# Obtain time indices of sub cubes and extract data
		#---------------------------------------------------
		print('Dividing into sub cubes if windowSize < timeGrid')

		ind_time  = [wind_ind*windowSize,wind_ind*windowSize + windowSize]
		DATA_xt   = observable[...,ind_time[0]:ind_time[1]]
		if DATA_xt.shape[-1] != windowSize:
			break
		Tgrid_new = Tgrid[ind_time[0]:ind_time[1]]
		Tgrid_c_grid.append(np.mean(Tgrid_new))


		#---------------------------------------------------
		# decompose into l m
		#---------------------------------------------------
		print('Performing spherical transform')
		DATA_kt = []
		Ntimes = len(Tgrid_new)
		
		# Define function to compute spherical harmoic coefficients in parallel
		def computeSHT(DATA,null):
			proj = projectOnSphericalHarmonics(DATA,np.arange(lMax+1),np.arange(-lMax,lMax+1),\
				                                theta=Lat_HG*np.pi/180.+np.pi/2,phi=Lon_HG*np.pi/180.,normalized=True,\
				                                axisTheta=-2,axisPhi=-1,pgBar=False,nbCores=1)
			return proj

		# Perform calculation in parallel
		DATA_kt = reduce(computeSHT,(DATA_xt,0),DATA_xt.shape[-1],96,None,True)

		#---------------------------------------------------
		# fourier transform in time
		#---------------------------------------------------
		print('Performing Fourier Transform in time')
		Tgrid_tmp = (Tgrid_new - Tgrid_new[0])*3600 *24*365.25
		omegaGrid = np.fft.fftshift(np.fft.fftfreq(Ntimes,np.diff(Tgrid_tmp)[0]),axes=-1)
		DATA_kw = np.fft.fftshift(np.fft.ifft(DATA_kt,axis=-1),axes=-1) * len(Tgrid_tmp) * np.diff(Tgrid_tmp)[0]
		
		#---------------------------------------------------
		# Power spectrum
		#---------------------------------------------------
		print('Compute Power Spectrum')
		POW = np.nan_to_num(abs(DATA_kw*np.conj(DATA_kw)))

		# extract the channel for plotting
		POW2 = POW[:,np.arange(-lMax,lMax+1)>=0]
		POW_tmp = []
		for ii in range(POW2.shape[-1]):
			POW_tmp .append(np.diag(POW2[...,ii],-mshift))
		POW2 = np.array(POW_tmp).T




		#---------------------------------------------------
		# plot function of Power vs (nu and m) for ell==mm (each depth)
		#---------------------------------------------------
		print('Plotting for ell = m+%i' % mshift)
		# Compute normalization for plotting purposes
		indl = np.argmin(abs(omegaGrid*1e9 +300)); indu = np.argmin(abs(omegaGrid*1e9-100))
		norms = np.mean(POW2[:,indl:indu],axis=1)

		plt.figure()
		vmax = 0.1*np.amax((POW2/norms[:,None]).T[np.argmin(abs(omegaGrid*1e9+500)):np.argmin(abs(omegaGrid*1e9-100)),:])
		plt.pcolormesh(np.arange(POW2.shape[0])-0.5,omegaGrid*1e9,(POW2/norms[:,None]).T,cmap='binary',vmax=vmax)
		plt.plot(np.arange(lMax+1),-2*(453.1e-9)/(np.arange(lMax+1)+1)*1e9,'r',linewidth=2)
		plt.plot(np.arange(lMax+1),-3*2*(453.1e-9)/(np.arange(lMax+1)+1)*1e9,'g',linewidth=2)
		plt.ylim([-500,500])
		plt.xlim([0,lMax ])
		plt.ylabel('Frequency [nHz]',fontsize=20)
		plt.xlabel('Azimuthal order $m$',fontsize=20)
		plt.tight_layout()

		POWS.append(POW)

	POWS = np.array(POWS)

	print("Saving")
	np.savez_compressed(outdir + '/PowerSpectra_%s_[%1.2fyears]_%iSpect%s_FULL_PMODE%i.npz' % (OBS.upper(),Tgrid_tmp[-1] / (3600*24*365.25),len(Tgrid_c_grid),'_HMI5deg',PMODE),\
											freqGrid = omegaGrid,Lgrid = np.arange(lMax+1),\
											Mgrid = np.arange(-lMax,lMax+1),\
											POWS = POWS)
	print("Done")
