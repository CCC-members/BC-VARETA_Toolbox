# BC-VARETA toolbox

Tool for MEEG data processing based on Brain Connectivity Variable Resolution Tomographic Analysis (BC-VARETA) Model. 
See description of BC-VARETA and example in simulations at the link https://github.com/dpazlinares/BC-VARETA.

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
  
Inputs for bash:
    - configure files:
        app/properties.json
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
