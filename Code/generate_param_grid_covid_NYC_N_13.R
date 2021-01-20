# ---- generate_profile_combinataions ----


# Header ------------------------------------------------------------------
## Name: generate_param_grid_covid_NYC_N_13.R
## Author: Rahul Subramanian
## Description: Creates 25,000 combination parameter grid via LHS for 
## Model N_13
rm(list = ls())

ptm <- proc.time()

#Load Libraries
source("load_libraries_essential.R")
source("rahul_theme.R")
library(pomp)
library(stringr)


args = commandArgs(trailingOnly=TRUE)

#model_name = "N_13"

model_name = as.character(args[1])
print(model_name)

require(tgp)

x <- lhs(25000, 
         rbind(
           c(5,5), # M_0
           c(13,13), # V_0
           c(14,14), # K_0
           c(2.0,8.0), # R_0
           c(0,1), # b_q
           c(0,0), # b_a
           c(0,1), # b_p
           c(0,1), # p_S
           c(.05,.40), # p_H_cond_S
           c(1.09,1.09), # phi_E
           c(1.09,1.09), # phi_U
           c(1/5,1/5), # phi_S
           c(1/8,1/8), # h_V
           c(1/1,1/5), # gamma
           c(8.0e6,8.0e6), # N_0
           c(0,2e4), #E_0
           c(0,2e4), #z_0
           c(0,0), #C_0
           c(1.7e+01,1.7e+01), #social_distancing_start_time
           c(2.2e+01, 2.2e+01), # quarantine_start_time
           c(9.0e-01,9.0e-01), # PCR_sens
           c(0,0.50), # sigma_M
           c(1.215073e-02 ,1.215073e-02 ), # beta_w_3
           c(9.810086e-01,9.810086e-01), # beta_w_2
           c(-3.723481e+01,-3.723481e+01), # beta_w_1
           c(2.294094e+02,2.294094e+02), # beta_w_0
           c(1.183300e+03,1.183300e+03), # g_0
           c(1.162005e-01,1.162005e-01), # g_F
           c(1.091121e+02,1.091121e+02), # sigma_epsilon
           c(1.62e-01,1.62e-01) # G_w_y_scaling
           )
)

names <- c("M_0","V_0", "K_0","R_0","b_q","b_a","b_p","p_S","p_H_cond_S","phi_E",
           "phi_U",
           "phi_S","h_V","gamma","N_0","E_0","z_0","C_0",
           "social_distancing_start_time","quarantine_start_time",
           "PCR_sens","sigma_M", "beta_w_3", "beta_w_2",
           "beta_w_1", "beta_w_0", "g_0", "g_F", "sigma_epsilon",
           "G_w_y_scaling")

x <- as.data.frame(x)
colnames(x) <- names

write.csv(x, file = paste0("../Generated_Data/Profile_Combination_Lists/",
                            model_name,"_Model/",
                            model_name,
                            "_param_grid.csv"),
          append = FALSE, row.names = FALSE)
proc.time() - ptm

