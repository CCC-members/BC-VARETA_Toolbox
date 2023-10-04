# BC-VARETA toolbox

Tool for MEEG data processing based on Brain Connectivity Variable Resolution Tomographic Analysis (BC-VARETA) Model. 

## Summary
Electrophysiological Source Imaging (ESI) is the technical term encompassing the type of reconstruction or inverse methods for the functional images of the brain or body defined as a current source profile causing their remote observations as a peripheral electromagnetic field profile adquired by sensors. ESI resolution in time and frequency is ideal and promises a mean to elucidate in these domains the intrincate mechanisms that govern brain function during resting state or task. inverse solutions for the corresponding inverse problem of electromagnetism
ESI    achieving a  a group of imaging methods with the aim of uncovering the mechanisms underpinnning brain function with appropiate temporal and spectral resolution.  and time : Brain Connectivity Variable Resolution Electromagnetic Tomographic Analysis (BC-VARETA). BC-VARETA is meant to be the  distribute our recent advances  developed methods on the third generation of nonlinear methods for MEEG Time Series analysis. Into the state of the art of MEEG analysis, the methodology underlying our tool (BC-VARETA) brings out several assets. First: Constitutes a truly Bayesian Identification approach of Linear Dynamical Systems in the Frequency Domain, grounded in more consistent models (third generation) for the joint nonlinear estimation of MEEG Sources Activity and Connectivity. Second: Achieves Super-Resolution, through the iterative solution of a Sparse Hermitian Sources Graphical Model that underlies the Connectivity Target Function. Third: Tackles efficiently in High Dimensional and Complex set up the estimation of connectivity, those constituting technical issues that challenge current MEEG source analysis methods. Fourth: Incorporates priors at the connectivity level by penalizing the groups of variables, corresponding to the Gray Matter anatomical segmentation, and including a probability mask of the anatomically plausible connections, given by synaptic transmission in the short-range (spatially invariant empirical Kernel of the connections strength decay with distance) and long-range (White Matter tracks connectivity strength from Diffusion Tensor Imaging). Along with the implementation of our method, we include in this toolbox a benchmark for the validation of MEEG source analysis methods, that would serve for the evaluation of sophisticated methodologies (third generation). It incorporates two elements. First: A realistic simulation framework, for the generation of MEEG synthetic data, given an underlying source connectivity structure. Second: Sensitive quality measures that allow for a reliable evaluation of the source activity and connectivity reconstruction performance, based on the Spatial Dispersion and Earth Movers’ Distance, in both source and connectivity space.

References:

    Paz-Linares, D., Gonzalez-Moreira, E., Martinez-Montes, E. and Valdes-Sosa, P.A., 2018. Note on the Estimation of Embedded Hermitian Gaussian Graphical Models for MEEG Source Activity and Connectivity Analysis in the Frequency Domain. Part I: Single Frequency Component and Subject. arXiv preprint arXiv:1810.01174. https://arxiv.org/abs/1810.01174
    Paz-Linares Deirel, Gonzalez-Moreira Eduardo, Areces-Gonzalez Ariosky, et al. (2023),Minimizing the distortions in electrophysiological source imaging of cortical oscillatory activity via Spectral Structured Sparse Bayesian Learning. URL: https://www.frontiersin.org/articles/10.3389/fnins.2023.978527 DOI:10.3389/fnins.2023.978527
   
    Paz-Linares, D., Gonzalez-Moreira, E., Martinez-Montes, E., Valdes-Hernandez, P.A., Bosch-Bayard, J., Bringas-Vega, M.L. and Valdes-Sosa, P.A., 2018. Caulking the Leakage Effect in MEEG Source Connectivity Analysis. arXiv preprint arXiv:1810.00786. https://arxiv.org/abs/1810.00786
    Paz-Linares D., Gonzalez-Moreira E., Areces-Gonzalez A., Wang Y., Li M., Martinez-Montes E., et al.. (2022). Identification of oscillatory brain networks with hidden gaussian graphical spectral models of EEG/MEG. ArXiv [Preprint]. arXiv:1810.01174. 10.48550/arXiv.1810.01174

CCC-members/BC-VARETA_Toolbox direct sourse:

    https://codeload.github.com/CCC-members/BC-VARETA_Toolbox/zip/master

Example of data structure (time series, leadfield, surface, and electrodes) is hosted in Onedrive:

    https://lstneuro-my.sharepoint.com/:u:/g/personal/cc-lab_neuroinformatics-collaboratory_org/EQVy7Y3oL9lDqS4_aNwglCsBMngspSuQ6yVudDj1xUOhgA?download=1

Main Function for MEEG real data analysis
    - Main.m      (**call this function for run**)

Execute BC-VARETA Toolbox
    - Main          % For graphical interface
    - Main nogui    % For batch processing
    - Main update   % Update to latest version     
  
Inputs for bash:
    - configure files:
        bcv_properties/general_params.json
        bcv_properties/sensor_params.json
        bcv_properties/activation_params.json
        bcv_properties/connectivity_params.json
 
   Outputs:
    - results: subfolder containing the bc-vareta outputs
  
Complementary Functions
    - xspectrum: computes the spectra of the simulated scalp activity 
    - bcvareta: executes BC-VARETA method
    - bcvareta_initial_values: computes 'bcvareta' initialization
    - screening_ssbl: extracts the posibly active generators as part of 'bcvareta_initial_values', using the Elastic Net Structured Sparse Bayesian Learning
    - trascendent_term: nonlinear function for regularization parameters estimation within the function 'screening_ssbl'     
    - screening: applies a smoothing to the outputs of 'screening_ssbl'


Authors:
   - Deirel Paz Linares
   - Eduardo Gonzalez Moreira
   - Ariosky Areces Gonzalez
   - Pedro A. Valdes Sosa

Date: September 15, 2018
