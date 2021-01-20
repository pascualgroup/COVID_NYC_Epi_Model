# ---- generate_profile_combinataions ----


# Header ------------------------------------------------------------------
## Name: generate_profile_combinations_covid_NYC_N_12.R
## Author: Rahul Subramanian
## Description: Creates 30*40-combination list for given by profile_var as 1st command line argument
rm(list = ls())

ptm <- proc.time()

#Load Libraries
source("load_libraries_essential.R")
source("rahul_theme.R")
library(pomp)
library(stringr)


args = commandArgs(trailingOnly=TRUE)

#model_name = "N_12"
#profile_var = "b_a"

profile_var = as.character(args[1])
print(profile_var)

model_name = as.character(args[2])
print(model_name)

#Load box
top_20_LL_box = read.csv(
  file = paste0("../Generated_Data/Profile_Combination_Lists/",
  model_name,
  "_Model/orignal_20_LL_param_box_from_1st_MIF_run.csv"))

#Modify G_w_y_scaling box boundaries
par_box_boundaries = top_20_LL_box %>%
  dplyr::select(-msg, -iter_num, -param_index, -loglik, -nfail, -trace_num,
                -loglist.se) 

if(profile_var == "G_w_y_scaling"){
  par_box_boundaries$G_w_y_scaling = c(0,0.33)
}else{
  if(profile_var == 'b_a'){
    par_box_boundaries$b_a = c(0,1)
    par_box_boundaries$b_p = c(0,1)
  }else{
    
  }
}



par_box_boundaries_clean = dplyr::select(par_box_boundaries, -one_of(profile_var) )
theta.t.lo = as.numeric(as.vector(par_box_boundaries_clean[1,]))
theta.t.hi = as.numeric(as.vector(par_box_boundaries_clean[2,]))
names(theta.t.lo) = colnames(par_box_boundaries_clean)
names(theta.t.hi) = colnames(par_box_boundaries_clean)

prof_var_boundaries = dplyr::select(par_box_boundaries, one_of(profile_var))
profileDesign(
  prof_var=seq(from=prof_var_boundaries[1,],to=prof_var_boundaries[2,],length=30),
  lower=theta.t.lo,upper=theta.t.hi,nprof=40
) -> pd
pd_col = colnames(pd)
colnames(pd) = c(profile_var, pd_col[2:length(pd_col)])

write.csv(pd, file = paste0("../Generated_Data/Profile_Combination_Lists/",
                            model_name,"_Model/", profile_var,"_",
                            model_name,
                            "_profile_combination_list.csv"),
          append = FALSE, row.names = FALSE)
proc.time() - ptm
