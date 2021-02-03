# COVID_NYC_Epi_Model

This repository contains the supporting data and code for the PNAS publication "Quantifying asymptomatic infection and transmission of COVID-19 in New York City using observed cases,serology, and testing capacity" by R. Subramanian, Q. He, and M. Pascual. 

#Organization
The "Code" folder contains code used to implement the model.
All raw data is stored in the Down_Data folder.
Data generated during the analysis is stored in the Generated_Data folder.
Figures generated during the analysis are stored in the Figures folder.


# Core Code Pipeline
"NYC_COVID_Data_Exploration.Rmd" -> "NYC_Covid_Model_SEIAR_and_SEPIAR.Rmd" ->  "Figures_3_and_4_Plots.Rmd" -> "Surface_Plot_Panel_Script_add_colorbar_2D.m"

# Data Processing Code
For code used to process the COVID-19 testing data, see the file "NYC_COVID_Data_Exploration.Rmd".

For code used to process the respiratory syndrome surveillance data and confirmed influenza cases, see the file "NYC_Flu_and_Syndrome_Surveillance_Analysis.Rmd".

# Model Implementation

The principal file ("NYC_Covid_Model_SEIAR_and_SEPIAR.Rmd") contains the code used to implement and fit the SEIAR and SEPIAR models, as well as to generate the panels for Figures 2 and 5.

# Plotting 

The file "Figures_3_and_4_Plots.Rmd" contains code used to generate the panels for Figures 3 and 4. 

The full Figure 3 is then plotted in Matlab using the code in "Surface_Plot_Panel_Script_add_colorbar_2D.m".

# Other Supporting Files
## Implementation of SEPIR model
The file "NYC_Covid_Model_SEPIR.Rmd" contains the code used to implement and fit the SEPIR model.

## Validation Analysis Files.
The file "Representative_Simulation_Analysis_N_12.Rmd" contains the code used to generate sample simulated trajectories for use in the validation analysis.

The file "Model_N_12_Simulated_Data_Big_b_a_param_Combine_Grid_Outputs.Rmd" has the code used to combine the results of the initial MIF run for the grid search using the SEIAR model fitting to the simulated trajectory from the  big b_a parameter combination. 

The file "Model_N_12_Simulated_Data_Small_b_a_param_Combine_Grid_Outputs.Rmd" has the code used to combine the results of the initial MIF run for the grid search using the SEIAR model fitting to the simulated trajectory from the  small b_a parameter combination. 

The file "Model_N_12_Combine_Profile_Outputs_Fit_to_Sim_Data_Big_b_a_param_combination.Rmd" has the code used to combine the results of the MIF run for the b_a profile using the SEP
IAR model fitting to the simulated trajectory from the  big b_a parameter combination. 

The file "Model_N_12_Combine_Profile_Outputs_Fit_to_Sim_Data_Small_b_a_param_combination.Rmd" has the code used to combine the results of the MIF run for the b_a profile using the SEP
IAR model fitting to the simulated trajectory from the  small b_a parameter combination. 

See the file "Plot_code_for_p_S_profile_serology_LL_sim_traj_vs_obs_fit.Rmd" for code used to plot the results of the validation analysis.

For the validation analysis code used for MIF fitting to big/small b_a parameter simulated trajectories for the SEIAR grid search and SEPIAR b_a profile search along with code used to generate the profile combinations for b_a profile, please see the last section of the main implementation file (NYC_Covid_Model_SEIAR_and_SEPIAR.Rmd).

## IFR Calculations
The file "IFR_Calculations" contains the code used to calculate the infection fatality rate.
