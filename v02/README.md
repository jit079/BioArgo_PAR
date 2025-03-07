# BioArgo_PAR_Reconstruction version 2

This is the code to reconstruct below surface PAR at different depth ${z}$, denoted as $PAR({z})$, using depth ${z}$ and the corresponding $E_{d}$ values at wavelengths ${\lambda_i}$ (where ${\lambda_i}$=380, 443, 490, and 555 nm). It employs a generalized additive model (GAM) with coefficients as functions of depth, i.e., $PAR(z)$= $\sum{f_i(z)*E_{d}(\lambda_i,z)}$. Theoretical uncertainty is associated to each PAR estimate. The code is only valid from below surface 0 to $${\color{red}200}$$ m. In the code, $E_{d}(\lambda_i,z)$ are expressed in ${mW/cm^2/{\mu}m}$ and $PAR({z})$ in ${{\mu}E/m^2/s}$.

## Files needed to run the program

- _reconstruct_par.py_: the main code to read in $E_{d}(z)$ and z to generate estimated PAR and the associated uncertainty, using the example data.
- _example_ folder: includes the example data file for demonstration.
- _aux_ folder: includes _coeff.txt_ and _LUT_par_uncertainty.txt_.

## Usage

To run the code, simple do:

```
python reconstruct_par.py

```
The estimated PAR will be displayed on screen and a plot of measured/modeled PAR vs depth as well as uncertainty will be generated in the _example_ folder. 


