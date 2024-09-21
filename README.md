# BioArgo_PAR_Reconstruction

This is the code to reconstruct below surface PAR at different depth ${z}$, denoted as $PAR({z})$, using depth ${z}$ and the correpsonding $E_{d}$ values at wavelengths ${\lambda_i}$ (where ${\lambda_i}$=380, 443, 490, and 555 nm). It employs a generalized additive model (GAM) with coefficients as functions of depth, i.e., $PAR(z)$= $\sum{f_i(z)*E_{d}(\lambda_i,z)}$. To each PAR estimate is associated the correponding theoretical uncertainty. The code is only valid from below surface 0 to 100 m. 

## Files needed to run the program

- _pygam_ folder: includes the necessary scripts of the pyGAM package (Servén D. and Brummitt C., 2018. pyGAM: Generalized Additive Models in Python. Zenodo. DOI: 10.5281/zenodo.1208723), which is used to develop GAM for PAR(z). Details can be viewed at https://github.com/dswah/pyGAM/blob/master/doc/source/index.rst and https://pygam.readthedocs.io/en/latest/#
- _gam_3u_3c_ folder: includes 100 GAMs for PAR and one GAM for uncertainty, in pickle file format.
- _reconstruct_par.py_: the main code to read in $E_{d}(z)$ and z to generate estimated PAR and the associated uncertainty, using the example data.
- _example folder_: includes the example data file for demonstration.

## Installation
The installation of dependencies can be done easily using Anaconda/Miniconda. The dependencies are defined in environment.yml. A new environment (python installation) called gam can be created as follows:
```
conda env create -f environment.yml
conda activate gam
```
## Usage

To run the code, simple do:

```
python reconstruct_par.py
```
The estimated PAR will be displayed on screen and a plot of measured/modeled PAR vs depth will be generated in the _example_ folder. Note that it may take a bit longer time to run the code for the first time after setting up the environement. After that it becomes fast.

**Make sure to activate the gam environment before running the code**. Once processing is done, you can choose to deactivate the enviroment by:
```
conda deactivate
```


