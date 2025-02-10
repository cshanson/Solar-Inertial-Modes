# Hanson et al. 2025: Solar Inertial Mode Codes
 These codes were developed iteratively through the Hanson 2020-2025 papers. The ones available on this repository are specific to the ring-diagram 5 deg tile analysis. Codes for RDA 15 degree and LCT are available upon request.

 ## Installation
I have provided the yml file for the conda environment I use for these codes. It is likely more than needed, so feel free to parse it down if you need to keep your environments tidy. In the Conda folder you will find the **conda_env_hanson.yml** file. Edit the first and last lines of the yml file to change the environment name and location. Create this conda environment through

`>conda env create -f conda_env_hanson.yml`

Activate the environment and you should be ready to go.

## Getting the data
The ring fit files are available online at [JSOC](http://jsoc.stanford.edu/ajax/lookdata.html?ds=hmi.rdVfitsc_fd05). However, there are too many to download online. There are two options available

1) If you have access to a NetDRMS data base (e.g. Stanford or MPS) you can access the files there. Please use this [script](0_ExtractFromNetDRMS.py) 
2) Otherwise, I have already run the extraction code and provide a tar ball of the files on Dropbox.

## Computational requirements
These procedures are best used on high performace computing systems. For instance, I use the NYUAD Jubail cluster, often requesting 50Gb in memory and around 27 cores. The codes are parallelized, and can be run on a few cores. However, the memory requirements may prohibit use on personal machines (e.g. laptops).

## Running the codes
I have synthesized the power spectrum codes into 5 scripts, that need to be run in order. I have provided a [jupyter notebook](Script_Notebook.ipynb) as well. Here is a basic description of each script
|Script Name | Description |
|---|---|
| [0_ExtractFromNetDRMS.py](0_ExtractFromNetDRMS.py) | Only run if you have access to a NetDRMS database with the `hmi.rdVfitsc_fd05` data series available. Script will read **all** the available mode fit files, and associated keywords, and store them in a h5py dictionary. A h5py file is generated for each Carrington rotation.|
|[1_Prepare_HMI_5deg.py](1_Prepare_HMI_5deg.py) | By providing the path to the h5py files, this script will parse and process the relevant flow data as well as the information relevant to the flow patch location (i.e. lat, lon, Rotation). Here is where you can choose which p modes to average (for u<sub>x</sub> and u<sub>y</sub>) and over what $\ell$ range.|
| [2_BuildFlowMaps_HMI_RDA_5deg.py](2_BuildFlowMaps_HMI_RDA_5deg.py)| With the flow information extracted from the ring fit data, this script creates the data cubes, or flow maps. For each observation time step (~8 hrs) we populate a grid in lat and lon with all the available flow tiles. This is repeated for each time step, generating a data cube (lat, lon, time). This is repeated for each p mode the user wishes to analyse. Each data cube is processed as per the procedure in [Hanson and Hanasoge 2024](https://ui.adsabs.harvard.edu/abs/2024PhFl...36h6626H/abstract). This includes the removal of annual and constant flow signals as specified by [Liang et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...626A...3L/abstract). Data cubes of the horizontal flows and their radial vorticity and horizontal divergence are stored for the next step.|
| [3_PowerSpect_HMI_RDA_5deg.py](3_PowerSpect_HMI_RDA_5deg.py)| This script reads in the data cubes generated in the previous step, computes the spherical harmonic $Y_\ell^m$ coefficients and stores the power spectrum of these coefficients. A plot of the $\ell=m$ power spectrum is also provided to ensure the routines are running correctly.|
| [4_SpectrumFitting_MCMC.py](4_SpectrumFitting_MCMC.py)| This final step should only be run if fitting the modes. **Note:** this is less of a black box than the previous scripts. Several parts of this script should be tailored to the user's requirements. The fitting routines are utilizing MCMC, so be careful in the priors, posterior, burn in etc.. Images of the fits and the MCMC performace will be generated for each $m$ in the grid chosen by the user. |

## Referencing this work
Users of the entire or part of the code should cite the following papers: [Hanson et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.109H/abstract), [Hanson et al. 2022](https://ui.adsabs.harvard.edu/abs/2022NatAs...6..708H/abstract), [Hanson and Hanasoge 2024](https://ui.adsabs.harvard.edu/abs/2024PhFl...36h6626H/abstract) and Hanson et al. 2025.
