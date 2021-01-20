


# ---- mif_code ----

# Header ------------------------------------------------------------------
## Name: MIF_run_Model_N_12.R
## Author: Rahul Subramanian
## Description: Runs parameter combinations on midway for profile from original param grid
## for Model N_12

rm(list = ls())
ptm <- proc.time()

#Load Libraries
source("load_libraries_essential.R")
source("rahul_theme.R")
library(pomp)

args = commandArgs(trailingOnly = TRUE)
#param_index = as.numeric(args[1]) + as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

profile_var = as.character(args[1])
print(profile_var)

model_name = as.character(args[2])
print(model_name)

#model_name = "N_12"
#profile_var = "b_a"
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



require(foreach)
require(doParallel)
require(deSolve)

#Core management
no_cores <- detectCores()
cat("no_cores = ", no_cores, "\n")
assinged_cores = 1
cat("assinged_cores = ", assinged_cores, "\n")

cl <- makeCluster(assinged_cores)
registerDoParallel(cl)


param_index = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print("param_index")
print(param_index)


##load(param_grid)
pd = read.csv(
  file = paste0(
    "../Generated_Data/Profile_Combination_Lists/",
    model_name,
    "_Model/",
    profile_var,
    "_",
    model_name,
    "_profile_combination_list.csv"
  ),
  header = TRUE
)
head(pd)

midway_max_jobs = 500
group_size = nrow(pd) / midway_max_jobs
start_index = (param_index - 1) * group_size + 1
end_index = param_index * group_size
Num_mif_runs_per_start = 1
param_data_subset_act = pd[start_index:end_index,]
param_data_subset = param_data_subset_act[rep(
  seq_len(nrow(param_data_subset_act)),
  each = Num_mif_runs_per_start),]


rw_sd_list_default = rw.sd(
  V_0 = 0,
  K_0 = 0,
  phi_E = 0,
  phi_S = 0,
  h_V = 0,
  p_S = 0.02,
  p_H_cond_S = 0.02,
  gamma = 0.02,
  social_distancing_start_time = 0,
  quarantine_start_time = 0,
  z_0 = ivp(0.02),
  E_0 = ivp(0.02),
  N_0 = ivp(0),
  C_0 = ivp(0),
  PCR_sens = 0,
  b_q = 0.02,
  b_a = 0.02,
  b_p = 0.02,
  R_0 = 0.02,
  sigma_M = 0.02,
  beta_w_3 = 0,
  beta_w_2 = 0,
  beta_w_1 = 0,
  beta_w_0 = 0,
  g_0 = 0,
  g_F = 0,
  sigma_epsilon = 0,
  G_w_y_scaling = 0.02)


get_rwsd = function(profile_var){
  if(profile_var == "G_w_y_scaling"){
    rw.sd = rw.sd(
      V_0 = 0,
      K_0 = 0,
      phi_E = 0,
      phi_S = 0,
      h_V = 0,
      p_S = 0.02,
      p_H_cond_S = 0.02,
      gamma = 0.02,
      social_distancing_start_time = 0,
      quarantine_start_time = 0,
      z_0 = ivp(0.02),
      E_0 = ivp(0.02),
      N_0 = ivp(0),
      C_0 = ivp(0),
      PCR_sens = 0,
      b_q = 0.02,
      b_a = 0.02,
      b_p = 0,
      R_0 = 0.02,
      sigma_M = 0.02,
      beta_w_3 = 0,
      beta_w_2 = 0,
      beta_w_1 = 0,
      beta_w_0 = 0,
      g_0 = 0,
      g_F = 0,
      sigma_epsilon = 0,
      G_w_y_scaling = 0,
      M_0 = 0,
      phi_U = 0)
  }else{
    if(profile_var  == "R_0"){
      rw.sd = rw.sd(
        V_0 = 0,
        K_0 = 0,
        phi_E = 0,
        phi_S = 0,
        h_V = 0,
        p_S = 0.02,
        p_H_cond_S = 0.02,
        gamma = 0.02,
        social_distancing_start_time = 0,
        quarantine_start_time = 0,
        z_0 = ivp(0.02),
        E_0 = ivp(0.02),
        N_0 = ivp(0),
        C_0 = ivp(0),
        PCR_sens = 0,
        b_q = 0.02,
        b_a = 0.02,
        R_0 = 0,
        sigma_M = 0.02,
        beta_w_3 = 0,
        beta_w_2 = 0,
        beta_w_1 = 0,
        beta_w_0 = 0,
        g_0 = 0,
        g_F = 0,
        sigma_epsilon = 0,
        G_w_y_scaling = 0.02,
        M_0 = 0,
        phi_U = 0,)
    }else{
      if(profile_var == "b_a"){
        rw.sd = rw.sd(
          M_0 = 0,
          V_0 = 0,
          K_0 = 0,
          phi_E = 0,
          phi_U = 0,
          phi_S = 0,
          h_V = 0,
          p_S = 0.02,
          b_p = 0.02,
          p_H_cond_S = 0.02,
          gamma = 0.02,
          social_distancing_start_time = 0,
          quarantine_start_time = 0,
          z_0 = ivp(0.02),
          E_0 = ivp(0.02),
          N_0 = ivp(0),
          C_0 = ivp(0),
          PCR_sens = 0,
          b_q = 0.02,
          b_a = 0,
          R_0 = 0.02,
          sigma_M = 0.02,
          beta_w_3 = 0,
          beta_w_2 = 0,
          beta_w_1 = 0,
          beta_w_0 = 0,
          g_0 = 0,
          g_F = 0,
          sigma_epsilon = 0,
          G_w_y_scaling = 0)
      }else{
          if(profile_var == "p_S"){
            rw.sd = rw.sd(
              V_0 = 0,
              K_0 = 0,
              phi_E = 0,
              phi_S = 0,
              h_V = 0,
              p_S = 0,
              p_H_cond_S = 0.02,
              b_p = 0.02,
              gamma = 0.02,
              social_distancing_start_time = 0,
              quarantine_start_time = 0,
              z_0 = ivp(0.02),
              E_0 = ivp(0.02),
              N_0 = ivp(0),
              C_0 = ivp(0),
              PCR_sens = 0,
              b_q = 0.02,
              b_a = 0.02,
              R_0 = 0.02,
              sigma_M = 0.02,
              beta_w_3 = 0,
              beta_w_2 = 0,
              beta_w_1 = 0,
              beta_w_0 = 0,
              g_0 = 0,
              g_F = 0,
              sigma_epsilon = 0,
              G_w_y_scaling = 0.02)
          }else{
            if(profile_var == "p_H_cond_S"){
              rw.sd = rw.sd(
                V_0 = 0,
                K_0 = 0,
                phi_E = 0,
                b_p = 0.02,
                phi_S = 0,
                h_V = 0,
                p_S = 0.02,
                p_H_cond_S = 0,
                gamma = 0.02,
                social_distancing_start_time = 0,
                quarantine_start_time = 0,
                z_0 = ivp(0.02),
                E_0 = ivp(0.02),
                N_0 = ivp(0),
                C_0 = ivp(0),
                PCR_sens = 0,
                b_q = 0.02,
                b_a = 0.02,
                R_0 = 0.02,
                sigma_M = 0.02,
                beta_w_3 = 0,
                beta_w_2 = 0,
                beta_w_1 = 0,
                beta_w_0 = 0,
                g_0 = 0,
                g_F = 0,
                sigma_epsilon = 0,
                G_w_y_scaling = 0.02)
            }else{
              if(profile_var == "E_0"){
                rw.sd = rw.sd(
                  V_0 = 0,
                  K_0 = 0,
                  phi_E = 0,
                  phi_S = 0,
                  h_V = 0,
                  p_S = 0.02,
                  p_H_cond_S = 0.02,
                  gamma = 0.02,
                  social_distancing_start_time = 0,
                  quarantine_start_time = 0,
                  z_0 = ivp(0.02),
                  E_0 = ivp(0),
                  N_0 = ivp(0),
                  C_0 = ivp(0),
                  PCR_sens = 0,
                  b_q = 0.02,
                  b_a = 0.02,
                  b_p = 0.02,
                  R_0 = 0.02,
                  sigma_M = 0.02,
                  beta_w_3 = 0,
                  beta_w_2 = 0,
                  beta_w_1 = 0,
                  beta_w_0 = 0,
                  g_0 = 0,
                  g_F = 0,
                  sigma_epsilon = 0,
                  G_w_y_scaling = 0.02)
              }else{
                  if(profile_var == "z_0"){
                    rw.sd = rw.sd(
                      V_0 = 0,
                      K_0 = 0,
                      phi_E = 0,
                      phi_S = 0,
                      h_V = 0,
                      p_S = 0.02,
                      b_p = 0.02,
                      p_H_cond_S = 0.02,
                      gamma = 0.02,
                      social_distancing_start_time = 0,
                      quarantine_start_time = 0,
                      z_0 = ivp(0),
                      E_0 = ivp(0.02),
                      N_0 = ivp(0),
                      C_0 = ivp(0),
                      PCR_sens = 0,
                      b_q = 0.02,
                      b_a = 0.02,
                      R_0 = 0.02,
                      sigma_M = 0.02,
                      beta_w_3 = 0,
                      beta_w_2 = 0,
                      beta_w_1 = 0,
                      beta_w_0 = 0,
                      g_0 = 0,
                      g_F = 0,
                      sigma_epsilon = 0,
                      G_w_y_scaling = 0.02)
                  }else{
                      if(profile_var == "gamma"){
                        rw.sd = rw.sd(
                          V_0 = 0,
                          K_0 = 0,
                          phi_E = 0,
                          phi_S = 0,
                          h_V = 0,
                          p_S = 0.02,
                          p_H_cond_S = 0.02,
                          b_p = 0.02,
                          gamma = 0,
                          social_distancing_start_time = 0,
                          quarantine_start_time = 0,
                          z_0 = ivp(0.02),
                          E_0 = ivp(0.02),
                          N_0 = ivp(0),
                          C_0 = ivp(0),
                          PCR_sens = 0,
                          b_q = 0.02,
                          b_a = 0.02,
                          R_0 = 0.02,
                          sigma_M = 0.02,
                          beta_w_3 = 0,
                          beta_w_2 = 0,
                          beta_w_1 = 0,
                          beta_w_0 = 0,
                          g_0 = 0,
                          g_F = 0,
                          sigma_epsilon = 0,
                          G_w_y_scaling = 0.02)
                      }else{
                        if(profile_var == "b_q"){
                          rw.sd = rw.sd(
                            V_0 = 0,
                            K_0 = 0,
                            phi_E = 0,
                            phi_S = 0,
                            h_V = 0,
                            p_S = 0.02,
                            p_H_cond_S = 0.02,
                            gamma = 0.02,
                            b_p = 0.02,
                            social_distancing_start_time = 0,
                            quarantine_start_time = 0,
                            z_0 = ivp(0.02),
                            E_0 = ivp(0.02),
                            N_0 = ivp(0),
                            C_0 = ivp(0),
                            PCR_sens = 0,
                            b_q = 0,
                            b_a = 0.02,
                            R_0 = 0.02,
                            sigma_M = 0.02,
                            beta_w_3 = 0,
                            beta_w_2 = 0,
                            beta_w_1 = 0,
                            beta_w_0 = 0,
                            g_0 = 0,
                            g_F = 0,
                            sigma_epsilon = 0,
                            G_w_y_scaling = 0.02)
                        }else{
                          stop("Profile var not specified in rwsd wrapper function")
                        }
                        
                      }
                  }
              }
            }
          }
      }
    }
  }
}

rw.sd = get_rwsd(profile_var = profile_var)

detail_log = FALSE

if (detail_log == TRUE) {
  detailed_log_file_name = paste0(
    "../Generated_Data/Profiles/",
    model_name,
    "_Model/",
    profile_var,
    "_Profile/Detailed_Log/log_file_subset_",
    param_index,
    ".txt"
  )
  write(file = detailed_log_file_name,
        paste0("Log generated on ", Sys.time(), " \n"),
        append = FALSE)
}


mif_single_subset_data <-
  foreach(
    i = 1:nrow(param_data_subset),
    .combine = rbind,
    .packages = c('pomp', 'dplyr'),
    .export = c(
      "rproc",
      "rmeas",
      "dmeas",
      "init",
      "paramnames",
      "statenames",
      "obsnames",
      "param_data_subset",
      "par_trans",
      "acumvarnames",
      "covar"
    )
  )  %dopar%
  {
    tryCatch({
      print(param_data_subset[i,])
      print("iter_num")
      print(i)
      print("param_index")
      print(param_index)
      params =  param_data_subset[i,]
      start = param_data_subset[i,]
      start$msg = "start"
      start$iter_num = i
      start$param_index = param_index
      seed <- round(runif(1, min = 1, max = 2 ^ 30))
      #seed = 565013131
      mif_single_param_output <- mif2(
        data = Observed_data,
        times = Observed_data$times,
        t0 = t0,
        seed = seed,
        rproc = pomp::euler(rproc, delta.t = 1),
        params = params,
        paramnames = paramnames,
        statenames = statenames,
        obsnames = obsnames,
        dmeas = dmeas,
        accumvars = acumvarnames,
        rinit = init,
        tol = 0,
        rmeas = rmeas,
        partrans = par_trans,
        covar = covar,
        start =  params,
        Np = 10000,
        Nmif = 50,
        cooling.fraction.50 = 0.5,
        rw.sd = rw.sd
      )
      
      
      first_trace_df = mif_single_param_output@traces %>%
        as.data.frame()
      
      first_trace_df$trace_num = seq(1:nrow(first_trace_df))
      # trace_df_ll = trace_df %>%
      #   dplyr::select(loglik, nfail)
      # trace_df_no_ll = trace_df %>%
      #   dplyr::select(-loglik, -nfail)
      # trace_df = trace_df_no_ll %>%
      #   mutate(nfail = trace_df_ll$nfail,
      #          loglik = trace_df_ll$loglik)
      first_trace_df$loglik
      first_trace_df$loglist.se = NA
      first_trace_df$iter_num = i
      first_trace_df$param_index = param_index
      first_trace_df$msg = "first_trace"
      
      mif_second_round = mif_single_param_output %>%
        mif2(Nmif = 50)
      
      second_trace_df = mif_second_round@traces %>%
        as.data.frame()
      
      second_trace_df$trace_num = seq(1:nrow(second_trace_df))
      
      second_trace_df$loglik
      second_trace_df$loglist.se = NA
      second_trace_df$iter_num = i
      second_trace_df$param_index = param_index
      second_trace_df$msg = "second_trace"
      
      ll <- tryCatch(
        replicate(n = 10, logLik(
          pfilter(
            data = Observed_data,
            times = Observed_data$times,
            t0 = t0,
            rprocess = pomp::euler(rproc, delta.t = 1),
            paramnames = paramnames,
            statenames = statenames,
            obsnames = obsnames,
            dmeas = dmeas,
            accumvars = acumvarnames,
            rinit = init,
            rmeas = rmeas,
            partrans = par_trans,
            covar = covar,
            format = "data.frame",
            Np = 50000,
            params = coef(mif_second_round)
          )
        )),
        error = function(e)
          e
      )
      
      fin  = mif_second_round %>% coef() %>% rbind() %>% as.data.frame()
      
      
      if (is(ll, "error")) {
      } else{
        ll_with_se = logmeanexp(ll, se = TRUE)
        
        if (detail_log == TRUE) {
          log_str = paste0(log_str,
                           "pfilter_warnings: \n ",
                           warnings(),
                           " \n Done with warnings \n")
        }
        
      }
      if (is.na(ll_with_se[[1]])) {
      } else{
        fin$loglik  = ll_with_se[[1]]
        fin$loglist.se = ll_with_se[[2]]
      }
      
      
      
      
      fin$iter_num = i
      fin$param_index = param_index
      
      fin$msg = "mif1"
      
      start_and_trace = bind_rows(start, first_trace_df)
      start_and_trace = bind_rows(start_and_trace, second_trace_df)
      bind_rows(start_and_trace, fin)
    },
    error = function (e) {
      warning("Inside error function")
      print("iter_num")
      print(i)
      print("param_index")
      print(param_index)
      start = param_data_subset[i,]
      start$msg = "start"
      start$iter_num = i
      start$param_index = param_index
      start$loglik = NA
      start$nfail = NA
      start$trace_num = NA
      start$loglist.se = NA
      
      fin = start
      fin$msg = conditionMessage(e)
      
      full_join(start, fin, by = names(start))
    })
  } -> res


output_name = paste(
  "../Generated_Data/Profiles/",
  model_name,
  "_Model/",
  profile_var,
  "_Profile/Subset_Outputs/",
  profile_var,
  "_",
  model_name,
  "_Profile_subset_",
  param_index,
  ".RData",
  sep = ""
)


if (detail_log == TRUE) {
  write(file = detailed_log_file_name, log_output, append = TRUE)
}

save(res, file = output_name)
res

proc.time() - ptm
