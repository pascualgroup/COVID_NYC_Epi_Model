# ---- sim_b_a_profile_peak_N_12 ----


#rm(list = ls())
ptm <- proc.time()

#Load Libraries
source("load_libraries_essential.R")
source("rahul_theme.R")
library(pomp)

args = commandArgs(trailingOnly = TRUE)
#param_index = as.numeric(args[1]) + as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#model_name = as.character(args[1])
#print(model_name)
profile_var = "b_a"
model_name = "N_12"
#param_index = 1
#i = 1
#Load Observed NYC case data
Observed_data = read.csv(paste0(
  "../Generated_Data/observed_data_",
  model_name, ".csv"))
head(Observed_data)

### Define start date
true_start_date = as.Date("2020-03-01")
t0 = 0
start_of_year = as.Date("2020-01-01")
first_saturday_in_year = as.Date("2020-01-04")

## Compartment/Queue Cohort Numbers
M = 5
V = 13
K = 14


#Declare Csnippets and data
source("Csnippet_nyc_coronavirus_model_N_12.R")


## Load NYC covariate data
covariate_df = read.csv(file =
                          paste0("../Generated_Data/covariate_data_",
                                 model_name, ".csv"))



### Create covariate table
covar=covariate_table(
  time=covariate_df$times,
  L_advanced_2_days=covariate_df$L_advanced_2_days,
  F_w_y = covariate_df$F_w_y,
  L_orig = covariate_df$L_orig,
  w = covariate_df$Week,
  y = covariate_df$Year,
  times="time"
)

param_index = 1


head(profile_peak_data_b_a)
##load(param_grid)
load(file = paste0(
  "../Generated_Data/Profiles/", model_name,
  "_Model/",
  profile_var,
  "_Profile/top_2_LL_of_",
  profile_var,
  "_profile.RData"))
profile_peak_data = profile_peak_data_b_a


midway_max_jobs = 1
group_size = nrow(profile_peak_data) / midway_max_jobs
start_index = (param_index - 1) * group_size + 1
end_index = param_index * group_size
Num_sim_runs_per_start = 1
top_2_LL_end_data_subset_act = profile_peak_data[start_index:end_index,]
top_2_LL_end_data_subset = top_2_LL_end_data_subset_act[rep(
  seq_len(nrow(top_2_LL_end_data_subset_act)),
  each = Num_sim_runs_per_start),]

## Load Antibdoy data
nyc_antibdoy_df = read.csv("../Generated_Data/antibody_data_from_nyc_study_with_RS_calc_CI.csv")
head(nyc_antibdoy_df)






# Top 2 LL




top_2_LL_end_subset_with_antibody_LL =
  data.frame(matrix(nrow = 0,
                    ncol = ncol(top_2_LL_end_data_subset) + 5))
colnames(top_2_LL_end_subset_with_antibody_LL) = 
  c(colnames(top_2_LL_end_data_subset), "Antibody_Mean_LL", "Antibody_LL_SE","Median_Herd_Immunity",
    "sim_subset_index", "combo_num")

all_combo_data = data.frame(matrix(nrow = 0, ncol = 6))
colnames(all_combo_data) = c("time", "sim_data_median ",  "sim_data_low_Q",
                             "sim_data_high_Q","combo_num", "sim_subset_index")
all_combo_S_data = data.frame(matrix(nrow = 0, ncol = 6))
colnames(all_combo_S_data) = c("time", "sim_data_S_over_N_median ",  "sim_data_S_over_N_low_Q",
                               "sim_data_S_over_N_high_Q","combo_num", "sim_subset_index")
all_combo_beta_t_data = data.frame(matrix(nrow = 0, ncol = 4))
colnames(all_combo_beta_t_data) = c("time", "sim_data_beta_t_median ",
                                    "combo_num", "sim_subset_index")

all_combo_C_Q1_data = data.frame(matrix(nrow = 0, ncol = 6))
colnames(all_combo_C_Q1_data) = c("time", "sim_data_C_Q1_median ",  "sim_data_C_Q1_low_Q",
                                  "sim_data_C_Q1_high_Q","combo_num", "sim_subset_index")

all_combo_R_data = data.frame(matrix(nrow = 0, ncol = 6))
colnames(all_combo_R_data) = c("time", "sim_data_R_over_N_median ",  "sim_data_R_over_N_low_Q",
                               "sim_data_R_over_N_high_Q","combo_num", "sim_subset_index")

  
  
  for(combo_index in seq(1:nrow(top_2_LL_end_data_subset))){
      #print(combo_index)
    
    combo_params = top_2_LL_end_data_subset[combo_index,]
    combo_params = dplyr::select(combo_params, -one_of(
      "msg", "iter_num", "param_index","loglik", "nfail", "trace_num", "loglist.se"))
    sim_data_sample_param = simulate(nsim = 100,
                                     seed = 12345,
                                     times = Observed_data$times,
                                     t0 = t0,
                                     rprocess = pomp::euler(rproc,delta.t = 1),
                                     params = combo_params,
                                     paramnames = paramnames,
                                     statenames = statenames,
                                     obsnames = obsnames,
                                     accumvars = acumvarnames,
                                     rinit = init,
                                     rmeas = rmeas,
                                     covar = covar,
                                     partrans = par_trans,
                                     format = "data.frame")
    #head(sim_data)
    sim_data_sample_param_median_Y = aggregate(Y ~ time, sim_data_sample_param, median)
    sim_data_sample_param_quant = aggregate(Y ~ time, sim_data_sample_param, quantile, probs = c(0.025, 0.975))
    sim_data_sample_param_quant$Y = as.data.frame(sim_data_sample_param_quant$Y)
    colnames(sim_data_sample_param_quant$Y) = c("Q2.5", "Q97.5")
    
    combo_num = rep(combo_index, nrow(sim_data_sample_param_median_Y))
    sim_subset_index = rep(param_index, nrow(sim_data_sample_param_median_Y))
    single_combo_data = data.frame(time =  sim_data_sample_param_median_Y$time,
                                   sim_data_median = sim_data_sample_param_median_Y$Y,
                                   sim_data_low_Q = sim_data_sample_param_quant$Y$Q2.5,
                                   sim_data_high_Q = sim_data_sample_param_quant$Y$Q97.5,
                                   combo_num = combo_num,
                                   sim_subset_index = sim_subset_index)
    all_combo_data = rbind(all_combo_data, single_combo_data)
    
    sim_data_sample_param$S_over_N = sim_data_sample_param$S/sim_data_sample_param$N
    
    sim_data_S_over_N_median = aggregate(S_over_N ~ time, sim_data_sample_param, median)
    sim_data_sample_param_S_over_N_quant = aggregate(S_over_N ~ time, sim_data_sample_param, quantile, probs = c(0.025, 0.975))
    sim_data_sample_param_S_over_N_quant$S_over_N = as.data.frame(sim_data_sample_param_S_over_N_quant$S_over_N)
    colnames(sim_data_sample_param_S_over_N_quant$S_over_N) = c("Q2.5", "Q97.5")
    
    
    sim_data_sample_param_S_over_N_quant = aggregate(S_over_N ~ time, sim_data_sample_param, quantile, probs = c(0.025, 0.975))
    sim_data_sample_param_S_over_N_quant$S_over_N = as.data.frame(sim_data_sample_param_S_over_N_quant$S_over_N)
    colnames(sim_data_sample_param_S_over_N_quant$S_over_N) = c("Q2.5", "Q97.5")
    
    
    
    
    single_combo_S_data = data.frame(time =  sim_data_sample_param_median_Y$time,
                                     sim_data_S_over_N_median = sim_data_S_over_N_median$S_over_N,
                                     sim_data_S_over_N_low_Q = sim_data_sample_param_S_over_N_quant$S_over_N$Q2.5,
                                     sim_data_S_over_N_high_Q = sim_data_sample_param_S_over_N_quant$S_over_N$Q97.5,
                                     combo_num = combo_num,
                                     sim_subset_index = sim_subset_index)
    all_combo_S_data = rbind(all_combo_S_data, single_combo_S_data)
    
    sim_data_beta_t_median = aggregate(beta_t ~ time, sim_data_sample_param, median)
    single_combo_beta_t_data = data.frame(time =  sim_data_sample_param_median_Y$time,
                                     sim_data_beta_t_median = sim_data_beta_t_median$beta_t,
                                     combo_num = combo_num,
                                     sim_subset_index = sim_subset_index)
    all_combo_beta_t_data = rbind(all_combo_beta_t_data, single_combo_beta_t_data)

    sim_data_C_Q1_median = aggregate(C_Q1 ~ time, sim_data_sample_param, median)
    sim_data_sample_param_C_Q1_quant = aggregate(C_Q1 ~ time, sim_data_sample_param, quantile, probs = c(0.025, 0.975))
    sim_data_sample_param_C_Q1_quant$C_Q1 = as.data.frame(sim_data_sample_param_C_Q1_quant$C_Q1)
    colnames(sim_data_sample_param_C_Q1_quant$C_Q1) = c("Q2.5", "Q97.5")
    
    single_combo_C_Q1_data = data.frame(time =  sim_data_sample_param_median_Y$time,
                                        sim_data_C_Q1_median = sim_data_C_Q1_median$C_Q1,
                                        sim_data_C_Q1_low_Q = sim_data_sample_param_C_Q1_quant$C_Q1$Q2.5,
                                        sim_data_C_Q1_high_Q = sim_data_sample_param_C_Q1_quant$C_Q1$Q97.5,
                                        combo_num = combo_num,
                                        sim_subset_index = sim_subset_index)
    all_combo_C_Q1_data = rbind(all_combo_C_Q1_data, single_combo_C_Q1_data)
    
    rel_columns = sim_data_sample_param %>%
      dplyr::select(R_A, R_F, R_H, time, .id, N)
    
    sim_data_sample_param_modified = rel_columns %>%
      mutate(R_sum = R_A + R_F + R_H)
    
    sim_data_sample_param_modified$R_over_N = sim_data_sample_param_modified$R_sum/sim_data_sample_param_modified$N 

    sim_data_R_over_N_median = aggregate(R_over_N ~ time, sim_data_sample_param_modified, median)
    sim_data_sample_param_R_over_N_quant = aggregate(R_over_N ~ time, sim_data_sample_param_modified,
                                                     quantile, probs = c(0.025, 0.975))
    sim_data_sample_param_R_over_N_quant$R_over_N = as.data.frame(sim_data_sample_param_R_over_N_quant$R_over_N)
    colnames(sim_data_sample_param_R_over_N_quant$R_over_N) = c("Q2.5", "Q97.5")
    single_combo_R_data = data.frame(
      time =  sim_data_sample_param_median_Y$time,
      sim_data_R_over_N_median = sim_data_R_over_N_median$R_over_N,
      sim_data_R_over_N_low_Q = sim_data_sample_param_R_over_N_quant$R_over_N$Q2.5,
      sim_data_R_over_N_high_Q = sim_data_sample_param_R_over_N_quant$R_over_N$Q97.5,
      combo_num = combo_num,
      sim_subset_index = sim_subset_index)
    all_combo_R_data = rbind(all_combo_R_data, single_combo_R_data)
    

    nyc_antibody_df = nyc_antibdoy_df %>%
      mutate(time = times)
    
    
    sim_data_sample_param_for_antibody_comp = sim_data_sample_param_modified %>%
      dplyr::select(time, R_over_N, sim_id = .id)
    
    sim_data_sample_param_with_antibody_df = inner_join(
      sim_data_sample_param_for_antibody_comp,
      nyc_antibody_df,
      by = c("time"))
    

    ### Exclude first antibody observation on March 1st-The simulation just started
    #on that date.
    sim_data_sample_param_with_antibody_df = sim_data_sample_param_with_antibody_df %>%
      filter(time > 0)
    
    ## Calculate LL
    sim_data_sample_param_with_antibody_df = sim_data_sample_param_with_antibody_df %>%
      mutate(Antibody_LL = dbinom(x = Num_Positive, p = R_over_N, size = Num_Sampled,
                                  log = TRUE))
    
    antibody_LL_per_sim_run = sim_data_sample_param_with_antibody_df %>%
      group_by(sim_id) %>%
      summarize(LL_per_run = sum(Antibody_LL)) %>%
      as.data.frame()
    
    total_antibody_LL_for_combination = logmeanexp(antibody_LL_per_sim_run$LL_per_run,
                                                   se = TRUE)
    single_param_with_antibody_LL = top_2_LL_end_data_subset[combo_index,]
    single_param_with_antibody_LL$Antibody_Mean_LL = total_antibody_LL_for_combination[[1]]
    single_param_with_antibody_LL$Antibody_LL_SE = total_antibody_LL_for_combination[[2]]
    single_param_with_antibody_LL$Median_Herd_Immunity =
      sim_data_R_over_N_median$R_over_N[nrow(sim_data_R_over_N_median)]
    single_param_with_antibody_LL$combo_num = combo_index
    single_param_with_antibody_LL$sim_subset_index = param_index
    
    top_2_LL_end_subset_with_antibody_LL = rbind(top_2_LL_end_subset_with_antibody_LL,
                                                 single_param_with_antibody_LL)

    
  }

  



save(all_combo_data,
     file = paste0(
       "../Generated_Data/Profiles/",
       model_name, "_Model/", profile_var, "_Profile/", profile_var,
       "_profile_top_2_LL_all_params_sim_cases_data.RData"))


save(all_combo_S_data,
     file = paste0(
       "../Generated_Data/Profiles/",
       model_name, "_Model/", profile_var, "_Profile/", profile_var,
       "_profile_top_2_LL_all_params_sim_S_over_N_data.RData"
       ))

save(all_combo_beta_t_data,
     file = paste0(
       "../Generated_Data/Profiles/",
       model_name, "_Model/", profile_var, "_Profile/", profile_var,
       "_profile_top_2_LL_all_params_sim_beta_t_data.RData"
     ))

save(all_combo_R_data,
     file = paste0(
       "../Generated_Data/Profiles/",
       model_name, "_Model/", profile_var, "_Profile/", profile_var,
       "_profile_top_2_LL_all_params_sim_R_over_N_data.RData"
     ))

save(all_combo_C_Q1_data,
     file = paste0(
       "../Generated_Data/Profiles/",
       model_name, "_Model/", profile_var, "_Profile/", profile_var,
       "_profile_top_2_LL_all_params_sim_C_Q_1_data.RData"))


save(top_2_LL_end_subset_with_antibody_LL,
     file = paste0(
       "../Generated_Data/Profiles/",
       model_name, "_Model/", profile_var, "_Profile/", profile_var,
       "_profile_top_2_LL_all_params_with_antibody_LL.RData"))



