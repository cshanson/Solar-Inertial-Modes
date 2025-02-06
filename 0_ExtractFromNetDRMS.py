#!/opt/local/anaconda/anaconda-2.2.0/bin/python2.7

#########################################################

# Routine to combine all 5deg data into single files for each Carrington Rotation
# Since there ~15 million files, this saves on file number burdens
# Run this routine on a system that has NetDRMS installed and the data (e.g. Stanford or MPS)

#########################################################

import numpy as np
import h5py

# Specify the directory where you want the combined files to be stored
destDIR = '/scratch/seismo/hanson/ToDalma/5degTiles/'

# In a for loop over the carrington rotations you want (cr)
for cr in range(2123,2270):
        print('Export CR %04d' % cr)

        # Call the NetDRMS routine show info to get the primary keys and the file location
        subP = subprocess.Popen('show_info hmi.rdvfitsc_fd05[%i][][][] -kP' % cr,\
                                      shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err  = subP.communicate()

        # Store the keys and location into variables
        netDRMSkey = out.split()[1::3]
        netDRMSdir = out.split()[2::3]

        # Open a h5py file to store the data
        with h5py.File(destDIR + '/CR_%04d.h5' % cr,'w') as h5f:
                for ii in range(len(netDRMSdir)):
                        # If file does not exist, continue the loop
                        try:
                                tmp = np.genfromtxt(netDRMSdir[ii][6:] + '/fit.out')
                        except:
                                continue;

                        # Create a data set with the primary keywords and store the data
                        h5f.create_dataset(netDRMSkey[ii].split('fd05')[1],data = tmp)

