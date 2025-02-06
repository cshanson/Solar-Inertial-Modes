#!/home/ch3246/.conda/envs/dalma-python3-CPU/bin/python3

#########################################################

# This routine takes all the flow map tiles for a single n and puts them into a data cube or lat, lon and time
# Inclues removal of temporal variations (B angle etc)

#########################################################

from Routines import *
import numpy as np
import scipy.interpolate as interp


#---------------------------------------------------
# User inputs
#---------------------------------------------------

# Specify the PMODE, can be done interactively or as input on command line
# PMODE = int(sys.argv[1])
PMODE = 0

# Specify the MASK size, to remove data with angular distance MASK (in degree) from disk center
MASK = 65

# Define the path where you want the data cubes saved
outdir = "/scratch/ch3246/Private/InertialWaves/"

# Specify array of Carrington rotations to load
CarrRot = np.arange(2097,2270)

# specify a reference time
# time of first HMI RDA tile for CR 2097
T0 = dt_to_dec(datetime.datetime(2010,5,20,0,0,0)) 

#---------------------------------------------------
# Load in Data
#---------------------------------------------------
print(dashLine,'Loading in Data')


# initialize massive lists for storing lat, lon, time and flow data
CRs = [];CMLons = [];LONHGs = [];LATHGs = [];SHLons = []
DATA_UX = []; DATA_UY = []

PG = progressBar(len(CarrRot),'serial')
for CR in CarrRot:
	with np.load('/scratch/ch3246/Private/InertialWaves/5deg_data/processed/CR%04d_processed.npz' % CR) as npzFile:
		CMLons.append(npzFile['CMLon'])	# Longitude of central meridian
		LONHGs.append(npzFile['LonHG'])	# Carrington Longitude
		LATHGs.append(npzFile['LatHG'])	# Carrington Latitude
		SHLons.append(npzFile['LonCM'])	# Stoney Hurst Longitude

		tmp = npzFile['flows']			
		DATA_UX.append(tmp[:,PMODE,0])	# UX flows
		DATA_UY.append(tmp[:,PMODE,2])	# UY flows
	CRs.append(np.ones(len(CMLons[-1]))*CR) # Carrington Rotation number

	PG.update()
del PG

# currently is a list of lists (one for each CR). Flatten into single list

CRs 	= np.concatenate(CRs,axis=0)
CMLons 	= np.concatenate(CMLons,axis=0)
LONHGs 	= np.concatenate(LONHGs,axis=0)
LATHGs 	= np.concatenate(LATHGs,axis=0)
SHLons 	= np.concatenate(SHLons,axis=0)
DATA_UX = np.concatenate(DATA_UX,axis=0)
DATA_UY = np.concatenate(DATA_UY,axis=0)

# Raise exception if lists are wrong size, stops further erros
if len(DATA_UX) != len(DATA_UY):
	raise Exception('Different File numbers for UX and UY')


#---------------------------------------------------
# Build coord grids
#---------------------------------------------------
print(dashLine,'Build coordinate grids and initialize data cubes')

# Build grids off all unique values in loaded lists
CR_grid 	= np.unique(CRs); 
CMLon_grid 	= np.unique(CMLons)[::-1]; 
LONHG_grid 	= np.unique(LONHGs); 
LATHG_grid 	= np.unique(LATHGs); 
SHLon_grid 	= np.linspace(-177.5,180,144)

CR_meshgrid,CMLon_meshgrid = np.meshgrid(CR_grid,CMLon_grid,indexing='ij')
CR_meshgrid = CR_meshgrid.reshape(-1);CMLon_meshgrid = CMLon_meshgrid.reshape(-1)


# Now we initialize the data cube (lat, lon, time)
# Time is defined by the Carrington rotation and the central meridian
# Cubes are tracked for 1/72 of a Carrinton rotation. Every unique pair of CR and CM is a unique time for a patch
Ux_Table = np.zeros((len(LATHG_grid),len(LONHG_grid),len(CR_grid)*len(CMLon_grid))) * np.nan
Uy_Table = np.zeros((len(LATHG_grid),len(LONHG_grid),len(CR_grid)*len(CMLon_grid))) * np.nan

print("Done")

#---------------------------------------------------
# Populate the data cube in the Carrington coordinate system

# Note we could populate straight into stoney hurst, since thats where the preprocessing occurs
# But we keep a copy of the unprocessed carrington grid from reference.
#---------------------------------------------------
print(dashLine,'Populating data cubes (Carrington + Stoney)')
PG = progressBar(len(CRs)/100,'serial')

# Initialize lists to store indices in order to swap between stoney and carrington
ind_TIME = [];ind_LONHG = [];ind_LATHG = [];ind_LONSH = []

for i in range(len(CRs)):
	# Find the indices in the grids that correspond to the coords for this item in list
	ind1 = np.argmin(abs(LATHG_grid - LATHGs[i]))
	ind2 = np.argmin(abs(LONHG_grid - LONHGs[i]))
	ind3 = np.argmin(abs(SHLon_grid - SHLons[i]))
	ind4 = np.where((CR_meshgrid == CRs[i]) * (CMLon_meshgrid == CMLons[i]))[0][0]

	ind_LATHG.append(ind1); ind_LONHG.append(ind2)
	ind_LONSH.append(ind3); ind_TIME.append(ind4);

	# populate the data cube
	Ux_Table[ind1,ind2,ind4] = DATA_UX[i]
	Uy_Table[ind1,ind2,ind4] = DATA_UY[i]

	if not i%100:
		PG.update()
del PG

# initialize a data cube for in the stoney hurst Coords
Ux_Table_sh = np.zeros((len(LATHG_grid),len(SHLon_grid),len(CR_grid)*len(CMLon_grid))) * np.nan
Uy_Table_sh = np.zeros((len(LATHG_grid),len(SHLon_grid),len(CR_grid)*len(CMLon_grid))) * np.nan

# Populate the stoney hurst coords by translating the Carrington to stoney hurst
for i in range(len(ind_TIME)):
	Ux_Table_sh[ind_LATHG[i],ind_LONSH[i],ind_TIME[i]] = Ux_Table[ind_LATHG[i],ind_LONHG[i],ind_TIME[i]]
	Uy_Table_sh[ind_LATHG[i],ind_LONSH[i],ind_TIME[i]] = Uy_Table[ind_LATHG[i],ind_LONHG[i],ind_TIME[i]]



#---------------------------------------------------
# Remove the one year frequency
# See Liang et al. 2019

# We fit a sinusoidal function (in to) to remove annual and constant signals
#---------------------------------------------------
print(dashLine,'Removing One year frequency')

# Define a time grid in fractional years, stating from T0
Tgrid_year = np.arange(Ux_Table.shape[-1])/72*27.2753 / 365.25 + T0


Tref = dt_to_dec(datetime.datetime(2010,6,7))
Teq  = dt_to_dec(datetime.datetime(2009,12,21))


Ntimes = Uy_Table_sh.shape[-1]
Tgrid = np.linspace(0,CR_grid[-1]-CR_grid[0],Ntimes) * 27.27 / 365.25 + T0


# Functional form of Liang et al. 2019 temporal variation
def systematic_fit_ZCL(A,B,C,D,E,F,G):
	omega = 2*np.pi
	Tgrid_tmp = (Tgrid - Tref)
	Tgrid_tmp_eq = Tgrid-Teq

	return A + B*(Tgrid_tmp) + (C + D*(Tgrid_tmp))*np.sin(omega*Tgrid_tmp) + E*np.cos(omega*Tgrid_tmp) + F*np.sin(omega*365.25*Tgrid_tmp)+ G*np.cos(omega*365.25*Tgrid_tmp)

# Perform removal and interp on each lat and lon (in Stoney Hurst) as function of time
PG = progressBar(Ux_Table_sh.shape[0],'serial')
for j in range(Uy_Table_sh.shape[0]):
	for k in range(Uy_Table_sh.shape[1]):
		for FlowComp in [0,1]: 
			B          = [Ux_Table_sh,Uy_Table_sh][FlowComp][j,k,:]
			goodinds   = ~np.isnan(B)
			if np.sum(goodinds) > 3:
				# we will use least squares to solve Ax = B for coefficients of temporal function
				A = np.array([systematic_fit_ZCL(1,0,0,0,0,0,0),\
					          systematic_fit_ZCL(0,1,0,0,0,0,0),\
					          systematic_fit_ZCL(0,0,1,0,0,0,0),\
					          systematic_fit_ZCL(0,0,0,1,0,0,0),\
					          systematic_fit_ZCL(0,0,0,0,1,0,0),\
					          systematic_fit_ZCL(0,0,0,0,0,1,0),\
					          systematic_fit_ZCL(0,0,0,0,0,0,1)]).T
				Bls = B[goodinds]
				Als = A[goodinds]
				coeffs = np.linalg.lstsq(Als,Bls,rcond=None)[0]

				# Now subtract the fitted terms from the signal
				B = B - systematic_fit_ZCL(coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4],coeffs[5],coeffs[6])			

				# interpolate the time array for missing data
				B = interp.griddata((np.arange(len(B))[goodinds],), B[goodinds], (np.arange(len(B))), method='linear',fill_value=np.nan)

			[Ux_Table_sh,Uy_Table_sh][FlowComp][j,k,:] = B
	PG.update()
del PG


#---------------------------------------------------
# Interpolate across the map for each time step in case missing data persists
# Could be redudant given previous step, but maintained due to consistancy with Proxauf
#---------------------------------------------------

print(dashLine,'Interpolating across map for each time step')

PG = progressBar(Ntimes,'serial')
for i in range(Ntimes):
	for FlowComp in [0,1]:

		lats,lons = np.meshgrid(LATHG_grid,SHLon_grid,indexing='ij')
		lons = lons.ravel();lats = lats.ravel()
		data = [Ux_Table_sh,Uy_Table_sh][FlowComp][:,:,i].ravel()


		lons_g = lons[np.isnan(data) == False]
		lats_g = lats[np.isnan(data) == False]
		data_g = data[np.isnan(data) == False]

		if len(lons_g) > 3:
			[Ux_Table_sh,Uy_Table_sh][FlowComp][:,:,i] = interp.griddata((lats_g, lons_g), data_g, (LATHG_grid[:,None], SHLon_grid[None,:]), method='linear',fill_value=np.nan)
	PG.update()
del PG


#---------------------------------------------------
# cut to MASK degree from the disk center (fill nan)
#---------------------------------------------------

print(dashLine,'Applying Mask, removing data > %i degree from disk center' % MASK)
lats,lons = np.meshgrid(LATHG_grid,SHLon_grid,indexing='ij')
dists = np.sqrt(lats**2 + lons**2)
mask  = np.where(dists <= MASK, 1, np.nan)

Ux_Table_sh = Ux_Table_sh * mask[:,:,None]
Uy_Table_sh = Uy_Table_sh * mask[:,:,None]


#---------------------------------------------------
# Make a copy of the data cube prior to shifting is to Carrington
#---------------------------------------------------

import copy
Ux_sh_save = copy.copy(Ux_Table_sh)
Uy_sh_save = copy.copy(Uy_Table_sh)


#---------------------------------------------------
# shift to equitorial rotation rate
# Data is tracked at Carrington rotation rate
# We shift it onto grid consistant with equatorial rotation rate 453.1nHz
#---------------------------------------------------
print(dashLine,'Shifting from Carrington to the Equatorial Rotation Rate')


# Calucalte the shift for each time step
Tgrid_tmp   = (Tgrid - Tgrid[0]) * 365.25*24*3600 
shift       = -0.0178 * Tgrid_tmp *1e-6  *180/np.pi # deltaLon = -0.0178 uRads/s  #delta_omega_murad = -0.0178

# Caluclate the max shift possible, in order to pad and avoid wrapping at this stage
MaxShift    = np.amax(np.abs(shift))
if MaxShift > 100:
	NpixelsPad  = np.ceil(MaxShift/5*1.5).astype(int)
	Pads        = np.zeros((Uy_Table_sh.shape[0],NpixelsPad))*np.nan
else:
	NpixelsPad = 0

# Extend the Stoney grid to account for wrapping, then initialize the new data cube
SHLon_gridE = np.linspace(SHLon_grid[0] - (NpixelsPad)*5,SHLon_grid[-1] + NpixelsPad*5,len(SHLon_grid)+2*NpixelsPad)
Ux_Table_v3 = np.zeros((len(LATHG_grid),len(SHLon_gridE),len(Tgrid))) * np.nan
Uy_Table_v3 = np.zeros((len(LATHG_grid),len(SHLon_gridE),len(Tgrid))) * np.nan


# For each time and flow component, shift from stoney hurst the coordinates (carrington - Eq) as if we were tracking at Eq rotation
PG = progressBar(Ntimes,'serial')
for i in range(Ntimes):
	for FlowComp in [0,1]:
		lats,lons = np.meshgrid(LATHG_grid,SHLon_gridE,indexing='ij')
		lons = lons.ravel();lats = lats.ravel()
		if MaxShift > 100:
			data = np.concatenate((Pads,[Ux_Table_sh,Uy_Table_sh][FlowComp][:,:,i],Pads),axis=1).ravel()
		else:
			data = [Ux_Table_sh,Uy_Table_sh][FlowComp][:,:,i].ravel()


		lons_g = lons[np.isnan(data) == False]
		lats_g = lats[np.isnan(data) == False]
		data_g = data[np.isnan(data) == False]

		if len(lons_g) > 3:
			[Ux_Table_v3,Uy_Table_v3][FlowComp][:,:,i] = interp.griddata((lats_g, lons_g), data_g, (LATHG_grid[:,None], (SHLon_gridE+shift[i])[None,:]), method='linear',fill_value=np.nan)

	PG.update()
del PG



# Now since we padded the data, we need to bring back incase the data went from e.g. 345 deg to 15 deg. i.e. wrap the data
# wrap around the sphere
if MaxShift > 100:
	nan_ux                                       = (~np.isnan(Ux_Table_v3)).astype(int)
	Ux_Table_v3                                  = np.nan_to_num(Ux_Table_v3)
	Ux_Table_v3[:,-2*NpixelsPad:-NpixelsPad,:] 	+= Ux_Table_v3[:,:NpixelsPad,:]
	Ux_Table_v3[:,NpixelsPad:2*NpixelsPad,:]   	+= Ux_Table_v3[:,-NpixelsPad:,:]
	nan_ux[:,-2*NpixelsPad:-NpixelsPad,:]      	+= nan_ux[:,:NpixelsPad,:]
	nan_ux[:,NpixelsPad:2*NpixelsPad,:]        	+= nan_ux[:,-NpixelsPad:,:]

	Ux_Table_v3 = Ux_Table_v3 * np.where(nan_ux !=0 ,1,np.nan)

	nan_uy                                       = (~np.isnan(Uy_Table_v3)).astype(int)
	Uy_Table_v3                                  = np.nan_to_num(Uy_Table_v3)
	Uy_Table_v3[:,-2*NpixelsPad:-NpixelsPad,:] 	+= Uy_Table_v3[:,:NpixelsPad,:]
	Uy_Table_v3[:,NpixelsPad:2*NpixelsPad,:]   	+= Uy_Table_v3[:,-NpixelsPad:,:]
	nan_uy[:,-2*NpixelsPad:-NpixelsPad,:]      	+= nan_uy[:,:NpixelsPad,:]
	nan_uy[:,NpixelsPad:2*NpixelsPad,:]        	+= nan_uy[:,-NpixelsPad:,:]

	Uy_Table_v3 = Uy_Table_v3 * np.where(nan_uy !=0 ,1,np.nan)

	Ux_Table_v3 = Ux_Table_v3[:,NpixelsPad:-NpixelsPad,:]
	Uy_Table_v3 = Uy_Table_v3[:,NpixelsPad:-NpixelsPad,:]
	del nan_ux,nan_uy


# Clean up
Ux_Table_sh = Ux_Table_v3;Uy_Table_sh = Uy_Table_v3;
del Ux_Table_v3,Uy_Table_v3




#----------------------------------------------------
# put back on carrington grid
# Data is on stoney grid, with (carr - eq) shift. We bring back to carrington grid
# Data will be in the Equatorial frame hence forth
#----------------------------------------------------
print(dashLine,'Transform from StoneyHurst back to Carrington')


Ux_Table_v2 = np.zeros((len(LATHG_grid),len(LONHG_grid),len(Tgrid))) * np.nan
Uy_Table_v2 = np.zeros((len(LATHG_grid),len(LONHG_grid),len(Tgrid))) * np.nan


def SH_to_Carr(DATA,Tind):
	MAT   = np.array(ind_TIME)== Tind
	if np.sum(MAT) > 0:
		rolln = np.amax(np.array(ind_LONHG)[MAT] - np.array(ind_LONSH)[MAT])
		DATA = np.roll(DATA,rolln,axis=-1)
	return np.append(DATA,DATA[:,:,0][:,:,None],axis=2)
Ux_Table_v2,Uy_Table_v2 = reduce(SH_to_Carr, (np.array([Ux_Table_sh,Uy_Table_sh]),np.arange(Ntimes)), Ntimes, 16, type=None, progressBar=True)


#---------------------------------------------------
# Subtract Lon Mean
#---------------------------------------------------

print(dashLine,'Subtracting the mean in Lon')
Ux_Table_v3 = np.nan_to_num(Ux_Table_v2 - np.nanmean(Ux_Table_v2,axis=1)[:,None,:])
Uy_Table_v3 = np.nan_to_num(Uy_Table_v2 - np.nanmean(Uy_Table_v2,axis=1)[:,None,:])

#---------------------------------------------------
# Vorticity Maps (FDM) np.gradient
#---------------------------------------------------
print(dashLine,'Computing Radial Vorticity and Horizontal Divergence')
LONHG_grid_rad = LONHG_grid*np.pi/180; 
LATHG_grid_rad = LATHG_grid*np.pi/180 +np.pi/2
dLon = np.diff(LONHG_grid_rad)[0]
dLat = np.diff(LATHG_grid_rad)[0]


LATHG_grid_rad = (LATHG_grid)*np.pi/180
Divergence = (np.gradient(Uy_Table_v3*np.cos(LATHG_grid_rad)[:,None,None],dLat,axis=0)\
              + np.gradient(Ux_Table_v3,dLon,axis=1))/(RSUN*np.cos(LATHG_grid_rad)[:,None,None])

Vorticity = (-np.gradient(Ux_Table_v3*np.cos(LATHG_grid_rad)[:,None,None],dLat,axis=0)\
              + np.gradient(Uy_Table_v3,dLon,axis=1))/(RSUN*np.cos(LATHG_grid_rad)[:,None,None])



#----------------------------------------------
# save Flow data cubes for spectrum calculation and post processing
print(dashLine,'Saving Data')
np.savez_compressed(outdir + '/Ring_Maps_HMI_5deg_PMODE%i.npz' % PMODE,Ux = Ux_sh_save,\
																		Uy = Uy_sh_save,\
																		Ux_eqrot = np.nan_to_num(Ux_Table_v3 - np.nanmean(Ux_Table_v3,axis=1)[:,None,:]),\
																		Uy_eqrot = np.nan_to_num(Uy_Table_v3 - np.nanmean(Uy_Table_v3,axis=1)[:,None,:]),\
																		Lat_HG = LATHG_grid,Lon_HG = LONHG_grid,\
																		Lon_SH = SHLon_grid,\
																		Tgrid = Tgrid,\
																		Divergence = Divergence,\
																		Vorticity = Vorticity)
