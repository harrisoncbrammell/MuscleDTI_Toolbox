# Help for tpca_denoising

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md#1-Purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md#5-Acknowledgements)

## 1. Purpose 
This function performs 4-dimensional image denoising by exploiting data redundancy in the Principal Components Analysis (PCA) domain using threshold PCA denoising, as proposed by Henriques et al., Imaging Neuroscience, 2023 

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md)

## 2. Usage
The user inputs the image data, a mask defining the region to be denoised (optional), a window size for the PCA analysis (optional), and a matrix containing spatail noise variance estimates. The denosed data and the number of non-noise principal components are returned.

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md)

## 3. Syntax
[Signal, n_comps] = TPCA_denoising(data, mask, kernel, sampling, centering, noise_var) 

The input arguments are:
* <i>data:</i> [x, y, z, M] 4D DTI volume
* <i>mask:</i> (optional) region-of-interest (boolean)
* <i>kernel:</i> (optional) window size for PCA analyses
* <i>noise_var:</i> [x,y,z] data matrix containing the spatial noise variance estimates of the data

The output arguments are:
* <i>Signal:</i> [x, y, z, M] denoised data matrix
* <i>n_comps:</i> [x, y, z] number of non-noise pricipal components identified

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md)

## 4. Example code

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md)

## 5. Acknowledgements
Adapted from: MPPCAdenoising by Jelle Veraart (jelle.veraart@nyumc.org), Copyright (c) 2016 New York University and University of Antwerp and PCAdenoising by Rafael Neto Henriques (https://github.com/RafaelNH/PCAdenoising)
    
Permission is hereby granted, free of charge, to any non-commercial entity ('Recipient') obtaining a copy of this software and associated documentation files (the 'Software'), to the Software solely for non-commercial research, including the rights to use, copy and modify the Software, subject to the following conditions: 
   
1. The above copyright notice and this permission notice shall be included by Recipient in all copies or substantial portions of the Software. 
     
2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESOF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

3. In no event shall NYU be liable for direct, indirect, special, incidental or consequential damages in connection with the Software. Recipient will defend, indemnify and hold NYU harmless from any claims or liability resulting from the use of the Software by recipient. 
 
4. Neither anything contained herein nor the delivery of the Software to recipient shall be deemed to grant the Recipient any right or licenses under any patents or patent application owned by NYU. 
 
5. The Software may only be used for non-commercial research and may not be used for clinical care. 

6. Any publication by Recipient of research involving the Software shall cite the references listed below.

Final version written by: Roberto Pineda Guzman (roberto.guzman@carle.com)

REFERENCES
Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping using random matrix theory Magn. Res. Med., 2016, early view, doi: 10.1002/mrm.26059
Henriques, R.N.; Ianus, A.; Novello, L.,; Jovicich, J.; Jespersen, S.N.; Shemesh, N. Efficient PCA denoising of spatially correlated redundant MRI data. Imaging Neuroscience, 2023, doi:
10.1162/imag_a_0049
Manjon, J.V.; Coupe, P.; Concha, L.; Buades, A.; Collins, D.L.;
Robles, M. Diffusion weighted image denoising using overcomplete local PCA. PLoS One, 2013, doi: 10.1371/journal.pone.0073021

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-TPCA_denoising.md)

