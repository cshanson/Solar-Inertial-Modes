# Solar Inertial Modes
 These codes were developed iteratively through the Hanson 2020-2025 papers. The ones available on this repository are specific to the ring-diagram 5 deg tile analysis. Codes for RDA 15 degree and LCT are available upon request.

 ## Installation
I have provided the yml file for the conda environment I use for these codes. It is likely more than needed, so feel free to parse it done if you need to keep your environments tidy. In the Conda folder you will find the **conda_env_hanson.yml** file. Edit the first and last lines of the yml file to change the environment name and location. Create this conda environment through

`>conda env create -f conda_env_hanson.yml`

Activate the environment and you should be ready to go.

## Getting the data
The ring fit files are available online at [JSOC](http://jsoc.stanford.edu/ajax/lookdata.html?ds=hmi.rdVfitsc_fd05). However, there are too many to download online. There are two options available

1) If you have access to a NetDRMS data base (e.g. Stanford or MPS) you can access the files there. Please use this [script](0_ExtractFromNetDRMS.py) 
