# MuscleDTI_Toolbox
## A MATLAB Toolbox for Skeletal Muscle Diffusion Tensor MRI Fiber Tractography 

The MuscleDTI_Toolbox consists of a series of custom-written MATLAB functions for performing diffusion-tensor MRI fiber tractography in skeletal muscle. This README file contains
  1) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#1-acknowledgements)
  2) [License information](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#2-license-information)
  3) [A list of MATLAB requirements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#3-matlab-requirements)
  4) [A list of the conventions assumed regarding data acquisition](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#4-data-acquisition-conventions-assumed)
  5) [An overview of a typical workflow using the toolbox](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#5-overview-of-a-typical-workflow)
  6) [Installation instructions](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#6-installation-instructions)
  7) [Becoming a MuscleDTI_Toolbox contributor](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#7-Becoming-a-MuscleDTI_Toolbox-contributor)
  8) [Links to other resources in the toolbox and online](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#8-other-resources)

Help files for all of the functions described below are available [here](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/README.md).

## 1. Acknowledgements
The functions in this toolbox reflect the collective contributions of many individuals over many years, including: Adam Anderson, Amanda Buck, Crystal Coolbaugh, Bruce Damon, Zhaohua Ding, Mark George, John Gore, Hannah Kilpatrick, Anneriet Heemskerk, Melissa Hooijmans, Drew Lansdown, Bennett Landman, Carly Lockard, Justin Montenegro, Roberto Pineda Guzman, Tuhin Sinha, Kevin Wilson, and Xingyu Zhou. Details regarding authorship and other individual contributions are noted in each function.

This work was supported by NIH grants NIH/NIAMS R01 AR050101 and NIH/NIAMS R01 AR073831. By using this software, users agree to acknowledge NIH/NIAMS R01 AR073831 in presentations and publications, cite the paper in which this toolbox is described (see [Section 8B](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#8-other-resources)), and adhere to [NIH policies](https://grants.nih.gov/policy-and-compliance/policy-topics/public-access) regarding open access to their publications. 

## 2. License Information
This work is covered under a [GNU General Public License](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/LICENSE.md), v. 3 or later.

## 3. MATLAB Requirements
The functions have been tested using MATLAB v. 2019b.  The toolbox consists primarily of custom-written functions, but also calls some built-in MATLAB functions.

## 4. Data Acquisition Conventions Assumed 
In-plane Geometry
* In-plane field of view is square
* In-plane reconstructed imaging matrix size is square

Slice Geometry
* Slices were acquired in the axial anatomical plane
* Slices numbers ascend in the +Z direction 
* Slice geometries of structural and DT images match

Whole-image Properties
* Signal patterns of structural and DT images match
* The slices cover the entire muscle of interest
* Acquisition volumes of structural and DT images match

## 5. Overview of a Typical Workflow
Muscle DTI tractography includes pre-processing and fiber-tract processing steps, as described below.

### A. Pre-processing
Before performing fiber tractography, several pre-processing steps must be performed:
* <i>File input</i>: Depending on the image format, this could be accomplished using the built-in MATLAB function <i>dicomread</i>, the built-in MATLAB function <i>niftiread</i>, or custom-written functions for proprietary image formats such as PAR/REC or XML/REC (Philips Medical Systems). 
* <i>Concatenation of multiple image acqusitions into a single dataset</i>: Many users acquire DTI images in multiple stacks.  In thes cases, the image stacks will need to be concatenated in a manner consistent with the slice number conventions described above. 
* <i>Mask formation</i>: Several steps in the toolbox require the formation of a binary image mask of the muscle of interest. An image mask is a matrix having the same dimensions of the image(s) or interest. In the matrix, 1's are used to represent the location of structures of interest and 0's are used to indicate other regions.
* <i>Image registration</i>: Switching of the diffusion-encoding gradients induces distortions in the images. In the toolbox, these distortions are corrected using image registration; the example in this toolbox uses the Demons registration technique, called using the <i>imregdemons</i> function in MATLAB.
* <i>De-noising</i>: Some level of noise in the data is inevitable. This noise adds variability to the estimation of the diffusion tensor; at a sufficiently low signal-to-noise ratio, it can also add bias. Two options are presented here. One uses the custom-written function anisotropic smoothing method, [<i>aniso4d_smooth</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Preprocessing-Functions/aniso4D_smoothing.m). A Threshold-based Principal Components Analysis, [<i>TPCA_denoising</i>](Preprocessing-Functions/TPCA_denoising.m) is also available.
* <i>Estimation of the diffusion tensor throughout the muscle of interest</i>: The example given here uses a weighted least squares method to estimate the diffusion tensor that best matches the observed signals, given the diffusion-encoding matrix and diffusion encoding (b-) value. This is performed in the function [<i>signal2tensor2</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Preprocessing-Functions/signal2tensor2.m), help for which is found [here](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md).

Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Sample-Scripts) for a MATLAB script and [this link](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Preprocessing-Functions) for the custom-written MATLAB functions that perform these tasks. Other required functions are part of MATLAB's proprietary toolboxes and cannot be distributed here. 

Switching of the phase-encoding gradients can induce distortions in the images. These can be corrected.  If the user wishes to perform this correction, they must have obtained a second set of non-diffusion weighted images with the phase-encoding direction reversed. The eddy-current correction scheme illustrated here uses freeware called [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki). A MATLAB function that calls the correction scheme in FSL is located [here](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Preprocessing-Functions/eddy_correct.m), and help for this function is located [here](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md).

### B. Fiber-tract Processing
* <i>Define muscle boundaries</i>: Real muscle fibers are assumed to be contained entirely within a single muscle of interest. The <i>fiber_track</i> function therefore requires the user to input a binary image mask demarcating the muscle boundaries; this mask is used to prevent fiber tracts from exiting the muscle of interest. The function [<i>define_muscle</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/define_muscle.m) is used to define this mask. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md) for detailed help on this function. Other programs, such as [ITK-SNAP](http://www.itksnap.org/pmwiki/pmwiki.php), can also be used to define a mask (which is subsequently imported into MATLAB.
* <i>Define the seed points</i>: In DTI fiber-tracking, the tracts are propagated from a set of points, commonly called "seed points." The MuscleDTI_Toolbox presents several options for defining these points:
    * <i> Aponeurosis seeding</i>: The user can seed the points along the anatomical structure into which the muscle fibers insert (a flattened tendinous structure called an aponeurosis). We refer to this as APO seeding. The function [<i>define_roi</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/define_roi.m) is used to digitize the {row, column, slice} coordinates of the aponeurosis; these points are used to define the seed surface. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md) for detailed help on this function.
    * <i> Voxel-based seeding</i>: The user can seed the points at the edges of the muscle mask (EDGE) or in uniformly distributed points throughout the muscle mask (VXL).  These options are specified within the fiber tracking function, described next.
* <i>Generate the fiber tracts</i>: Fiber tracts are propagated from the seed points by integrating the first eigenvector of the diffusion tensor. The function <i>fiber_track</i> is used for this purpose. The user can select from several seeding methods, two propagation algorithms, and several methods for determining when to stop propagating a tract. The major output of this function is a matrix containing the {row, column, slice} coordinates of each point along each fiber tract.  The function [<i>fiber_track</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_track.m) calls the function [<i>retrieve_tensor</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/retrieve_tensor.m), which finds the diffusion tensor in each voxel of the image. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md) for detailed help on this function.
* <i>Smooth the fiber tracts</i>: Fiber tract points are subject to errors in position because of the presence of noise and artifacts in the images. To mitigate these effects, the function [<i>fiber_smoother</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_smoother.m) performs a polynomial fit to each fiber tract. This also allows the interpolation of the fiber tract positions at a resolution higher than the original tracts.  This step is not required, but is strongly recommended prior to calling the [<i>fiber_quantifier</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_quantifier.m) function. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md) for detailed help on this function.
* <i>Quantify the tracts' structural properties</i>: After the fiber tracts have been polynomial-fitted, their structural properties are quantified using the function [<i>fiber_quantifier</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_quantifier.m).  The properties quantified include the pennation angle, curvature, and length. These properties are calculated in a pointwise manner along the fiber tracts. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md) for detailed help on this function.
* <i> Uniformly sample the whole muscle volume</i>. We have recently implemented the farthest streamline sampling method, introduced by [Li et al. (2024)](https://pubmed.ncbi.nlm.nih.gov/38963745/), which allows the user to maximize the uniformity with which the fiber tracts sample the muscle volume by finding the set of N tracts that are maximally separated from each other. Follow this link to the [<i> far_stream_sampling</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/far_stream_sampling.m) function and this link to a detailed help file. 
* <i>Eliminate erroneous results</i>: The quantitative results are examined and obviously wrong results are eliminated from the dataset. The function [<i>fiber_goodness</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_goodness.m) eliminates tracts that ascend and descend (an error due to overfitting); that have architectural properties that exceed certain limits; and/or that vary greatly from their neighbors. A final dataset is calulated, including the mean properties for each tract and for the entire muscle. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md) for detailed help on this function.
* <i>Visualize the results</i>: At several of these stages, the results can be visualized using the function [<i>fiber_visualizer</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_visualizer.m). The user can select the mask, seed surface, and/or fiber tracts for display.  The user can also select which image slices to display for anatomical reference. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md) for detailed help on this function.

Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Sample-Scripts/process_tracts.m) for a MATLAB script and [this link](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Tractography-Functions) for the custom-written MATLAB functions that perform these tasks. Other required functions are part of MATLAB's proprietary toolboxes and cannot be distributed here.

## 6. Installation Instructions
We assume that you have installed MATLAB v. 2019b or later, including the Image Processing toolbox.

1) From this page, click the green Code button.
2) Select Download Zip.
3) On your system, navigate to your Downloads folder.  
4) Unzip the file folder, saving it in a directory of your own choosing.
5) In MATLAB, add a search path to the folders with the toolbox's custom functions. There are separate folders for pre-processing and fiber tractography functions. You may wish to rename the folders or even combine all functions into a single folder.

## 7. Becoming a MuscleDTI_Toolbox Contributor
We encourage users  to contribute to and improve the MuscleDTI_Toolbox and become part of our community.  Together, we can make the software into a useful and robust tool for muscle research. Whether you are a coding wizard or someone who has never opened MATLAB before this project, we can use your help. Here are a few ways how you can get involved:

### A. Do you like to code?
* Fork the repository through GitHub. 
* Make changes (e.g. bug fixes, feature enhancements) in MATLAB. 
* Send a pull request to share your change with the community - don't forget to include a detailed commit text!
* Please use [existing conventions and styles](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/General%20Documents/CODING_CONVENTIONS.md) for variable names, versioning, etc. to keep the code readable. 

By contributing your code, you agree to license your contribution under the terms of the GNU General Public License.

### B. Do you like to write?
* Write improvements to the project's documentation and wiki. 
* Create tutorial videos or tips to make the project work. 
* Develop a folder of examples to show how the MuscleDTI_Toolbox was used in your research. 

### C. Do you like to help people?
* Help answer questions for people on open Issues. 
* Offer to help mentor another member of the community on how to use the MuscleDTI_Toolbox or contribute to the project. 

### D. Bugs, Features, Etc.
Experiencing a bug?  Have an idea? Please submit a bug report or feature request [here](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Issues). You can also see a list of our planned feature additions/upgrades [here](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Issues/Planned-Upgrades.md).

Lastly, please be aware of and commit to following the [Code of Conduct](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/General%20Documents/CODE%20OF%20CONDUCT.md).

## 8. Other Resources
### A. Within the toolbox:
Links to scripts, functions, and their help files:
* [Here's a link to the sample scripts](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Sample-Scripts)
* [Here's a link to all of the preprocessing functions](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Preprocessing-Functions)
* [Here's a link to all of the tractography functions](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Tractography-Functions)
* [Here's a link to all of the help files](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Help)
* [Here's a link to code to run the simulations described in the Damon et al (2021) paper](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Simulations)


### B. External to the toolbox:
Please cite the following publication when using this toolbox:
* Damon BM, et al. (2021). [A MATLAB toolbox for muscle diffusion-tensor MRI tractography.](https://doi.org/10.1016/j.jbiomech.2021.110540) <i>Journal of Biomechanics</i>.

Several recent reviews on muscle DTMRI include:
* Damon BM, et al. (2017). [Skeletal muscle DT-MRI fiber tracking: Rationale, data acquisition and analysis methods, applications, and future directions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5136336/). <i>NMR in Biomedicine</i>.

* Oudeman J, et al. (2016). [Techniques and applications of skeletal muscle diffusion tensor imaging: A review](https://www.ncbi.nlm.nih.gov/pubmed/26221741). <i> Journal of Magnetic REsonance Imaging</i>.

Our most recent technical updates to the Toolbox are described in:
* Damon BM, et al. (2025). [A comparison of skeletal muscle diffusion tensor imaging tractography seeding methods](https://pubmed.ncbi.nlm.nih.gov/40660805/). <i> NMR in Biomedicine</i>. Introduces VXL and EDGE seeding methods; compares the outcomes of these methods and APO seeding.
* Pineda Guzman RA, et al. (2025). [Improved DTI-based skeletal muscle architecture estimation via diffusion weighted image denoising and first eigenvector map smoothing](https://pubmed.ncbi.nlm.nih.gov/40660805/). <i> NMR in Biomedicine</i>. Compares TPCA-based denoising with anisotropic smoothing, with or without additional first eigenvector field smoothing.
* Hooijmans MT, et al. (2025). [A registration strategy to characterize DTI-observed changes in skeletal muscle architecture due to passive shortening](https://pubmed.ncbi.nlm.nih.gov/40063556/). <i>PLoS One</i>. Demonstrates ability to use image registration to transform fiber-tracts from a reference state to another structural state.
