# BioArgo_PAR_Reconstruction version 2, Matlab code

This is the Matlab version of code to reconstruct below surface PAR at different depth ${z}$, denoted as $PAR({z})$, using depth ${z}$ and the corresponding $E_{d}$ values at wavelengths ${\lambda_i}$ (where ${\lambda_i}$=380, 443, 490, and 555 nm). 

## Files needed to run the program

- _b_spine_basis.m_: the function to generate spline basis based on the input depth z.
- _calc_uncertainty.m_: the function to calculate uncertainty based on depth z and PAR values.
- _reconstruct_par.m_: the function to read in $E_{d}(z)$ and z to generate estimated PAR and the associated uncertainty. It calls _b_spine_basis.m_ and _calc_uncertainty.m_
- _example.m_: the code to demonstrate the use of _reconstruct_par.m_

## Usage

To run the code, simply execute _example.m_ in Matlab.

The estimated PAR will be displayed on screen and a plot of measured/modeled PAR vs depth as well as uncertainty will be generated in the _example_ folder. 



