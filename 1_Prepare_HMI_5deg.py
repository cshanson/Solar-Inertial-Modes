#!/home/ch3246/.conda/envs/dalma-python3-CPU/bin/python3

#########################################################

# This routine reads in the NetDRMS data files constructed during extraction in 0_ExtractFromNetDRMS.py, 
# averages the ell dependancy and keywords and stores ready for constructing the data cubes

#########################################################

from Routines import *
import subprocess
import os
import h5py
import numpy as np

OVERWRITE = True



# Specify the directory where the data from netdrms is stored.
# Must be in h5 format as per script 0_ExtractFromNetDRMS.py
outdir = '/scratch/ch3246/Private/InertialWaves/5deg_data/'


# Obtain a list of all the files in the directory
subP = subprocess.Popen('ls %s/CR_*' % outdir,\
                              shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
out,err  = subP.communicate()
out = out.split()
out = [x.decode('utf-8') for x in out]



# Go through For loop for each file
for ii in range(len(out)):
    print('Processing %s' % out[ii])

    # If the file exists, and we don't want to overwrite, then continue
    if os.path.isfile(outdir + '/processed/CR%s_processed.npz' % out[ii].split('/')[-1].split('_')[-1][:-3]) and not OVERWRITE:
        print('CR %s already processed, moving on....' % out[ii].split('/')[-1].split('_')[-1][:-3])
        continue

    # Initialize lists which will contain the keywords of the data
    CMLon = []; LonHG = []; LatHG = []; LonCM = []
    with h5py.File(out[ii]) as h5f:
        # read keywords in
        keywords = [x for x in h5f.keys()]

        # initialize list to contain the ell averaged flows
        flows = []
        PG = progressBar(len(keywords),'serial')
        for jj in range(len(keywords)):

            # break the key word string down into a list of each key
            kw = keywords[jj].replace('[','').split(']')


            # load the data
            dat = np.array(h5f[keywords[jj]])

            # If data is missing in file (very rare) continue
            try:
                ngrid = dat[:,0]
            except:
                continue

            # Store the keywords into respective keyword lists
            CarrRot = int(kw[0])        # Carrington rotation
            CMLon.append(float(kw[1]))  # Central Meridian Longitude
            LonHG.append(float(kw[2]))  # Carrington Longitude 
            LatHG.append(float(kw[3]))  # Latitude
            LonCM.append(float(kw[4]))  # Stoneyhurst Longitude
            
            # extract the ell averaged mode fits for each n
            flows_nn = []
            for nn in range(3):

                ells = dat[ngrid == nn,1]   # ell
                ux   = dat[ngrid == nn,5]   # ux mode fits
                dux  = dat[ngrid == nn,6]   # error in ux
                uy   = dat[ngrid == nn,7]   # uy mode fits
                duy  = dat[ngrid == nn,8]   # error in uy

                # average over the ell range for ell > 500
                flows_nn.append([np.mean(ux[ells>500]),np.mean(dux[ells>500]),\
                                    np.mean(uy[ells>500]),np.mean(duy[ells>500])])
            flows.append(flows_nn)
            PG.update()
        del PG
    
    # Turn lists into numpy arrays
    CMLon = np.array(CMLon)
    LonHG = np.array(LonHG)
    LatHG = np.array(LatHG)
    LonCM = np.array(LonCM)

    flows = np.array(flows)

    # Saved the data in the processed subfolder of outdir
    mkdir_p(outdir + '/processed/')
    np.savez_compressed(outdir + '/processed/CR%i_processed.npz' % CarrRot,\
                        CMLon = CMLon,\
                        LonHG = LonHG,\
                        LatHG = LatHG,\
                        LonCM = LonCM,
                        flows = flows)










