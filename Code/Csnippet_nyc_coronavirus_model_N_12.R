# ---- statenames ----
statenames = c("S",
               sprintf("E_%d",1:M),
               "I_A", "I_P",
               "I_S_1", "I_S_2", "I_H",
               "R_A", "R_F", "R_H",
               "A_T","beta_t",
               sprintf("Q_1_%d",1:V),
               sprintf("P_Q1_%d",1:V),
               sprintf("total_samples_for_PCR_Testing_Q1_%d",1:V),
               sprintf("total_samples_for_PCR_Testing_lag_1_Q1_%d",1:V),
               sprintf("total_samples_for_PCR_Testing_lag_2_Q1_%d",1:V),
               "First_re_test_Q1",
               "Second_re_test_Q1",
               "Backlog_Queue_1",
               "total_samples_for_PCR_Testing_backlog_Q1",
               "total_samples_for_PCR_Testing_backlog_lag_1_Q1",
               "total_samples_for_PCR_Testing_backlog_lag_2_Q1",
               sprintf("Q_NC_%d",1:V),
               sprintf("total_samples_for_PCR_Testing_QNC_%d",1:V),
               sprintf("total_samples_for_PCR_Testing_lag_1_QNC_%d",1:V),
               sprintf("total_samples_for_PCR_Testing_lag_2_QNC_%d",1:V),
               "neg_samples_Q1",
               "neg_samples_Q3",
               "neg_samples_Q5",
               "total_neg_samples_all_queues",
               "Backlog_Queue_NC",
               "total_samples_for_PCR_Testing_backlog_QNC",
               "total_samples_for_PCR_Testing_backlog_lag_1_QNC",
               "total_samples_for_PCR_Testing_backlog_lag_2_QNC",
               "Q_2",
               "total_samples_for_PCR_Testing_Q2",
               sprintf("Q_3_%d",1:K),
               sprintf("P_Q3_%d",1:K),
               sprintf("total_samples_for_PCR_Testing_Q3_%d",1:K),
               sprintf("total_samples_for_PCR_Testing_lag_1_Q3_%d",1:K),
               sprintf("total_samples_for_PCR_Testing_lag_2_Q3_%d",1:K),
               "First_re_test_Q3",
               "Second_re_test_Q3",
               "Backlog_Queue_3",
               "total_samples_for_PCR_Testing_backlog_Q3",
               "total_samples_for_PCR_Testing_backlog_lag_1_Q3",
               "total_samples_for_PCR_Testing_backlog_lag_2_Q3",
               "Q_4",
               "total_samples_for_PCR_Testing_Q4",
               "infected_sample_size",
               "total_sample_size",
               "prob_infected_Q5",
               "total_samples_for_PCR_Testing_Q5",
               "N",
               "L_int",
               "L_1",
               "L_2",
               "L_3",
               "L_4",
               "Prop_Positive_Tests_Track",
               "C_Q1", "C_Q2",
               "C_Q3", "C_Q4",
               "C_QNC",
               sprintf("Y_Q1_%d",1:V), sprintf("Y_QNC_%d",1:V),
               sprintf("Y_Q3_%d",1:K), 
               "Y_Q1_backlog",
               "Y_QNC_backlog",
               "Y_Q3_backlog",
               "Y_sum",
               "G_w_y",
               "Y_Q5",
               "Neg_State_Value_Detected",
               "NAN_State_Value_Detected",
               "Error_Printing_Complete")
acumvarnames = c("C_Q1", "C_QNC", "C_Q2",
                 "C_Q3", "C_Q4",
                 sprintf("Y_Q1_%d",1:V), sprintf("Y_QNC_%d",1:V),
                 sprintf("Y_Q3_%d",1:K), 
                 "Y_Q1_backlog",
                 "Y_QNC_backlog",
                 "Y_Q3_backlog",
                 "Y_Q5",
                 "Y_sum", "total_neg_samples_all_queues")
obsnames = c("Y", "obs_prop_positive")

# ---- paramnames ----
paramnames = c("M_0","V_0", "K_0",
               "phi_U", "phi_E", "phi_S", "h_V",
               "p_S", "p_H_cond_S", 
               "gamma", 
               "quarantine_start_time",
               "social_distancing_start_time",
               "PCR_sens",
               "b_q", "b_a", "b_p",
               "R_0",
               "E_0", "z_0", "N_0", "C_0", "sigma_M",
               "beta_w_3", "beta_w_2",
               "beta_w_1", "beta_w_0",
               "g_0", "g_F",
               "sigma_epsilon",
               "G_w_y_scaling")

# ---- rproc ----
#Process model Csnippet with queue incorporated
rproc <- Csnippet("
                  //Declare Arrays
                  double rate_I_P[2], trans_I_P[2]; //Declare arrays for eulermultinom transitions from Exposed Compartment
                  double rate_I_S_1[2], trans_I_S_1[2]; //Declare arrays for eulermultinom transitions from symptomatic infection Compartment
                  
                  double multinom_output[8], multinom_prob[8]; //Declare arrays for multinom in Q5
                  
                  int M = (int) M_0; //Number of exposed compartments
                  int V = (int) V_0; //Number of days spent in hospital (number of cohorts in Queues 1 and NC)
                  int K = (int) K_0; //Number of days spent in quarantine (number of cohorts in Queue 3)
                  
                  int m; //Exposed compartment number
                  int v; //Queue 1/ Queue NC cohort number
                  int k; //Queue 3 cohort number
                  
                  double dE_m_E_m_1[M-1]; //Declare array for binomial transitions between Exposed Compartments
                  
                  double P_Q1_old[V]; //Declare array to store old value of P_Q1 during capacity calculations. 
                  double P_Q3_old[K]; //Declare array to store old value of P_Q3 during capacity calculations.
                  
                  double QNC_old[V]; //Declare array to store old value of QNC during capacity calculations. 
                  double Q1_old[V]; //Declare array to store old value of Q1 during capacity calculations. 
                  double Q3_old[K]; //Declare array to store old value of Q3 during capacity calculations.
                  
                  //Declare E pointer array
                  double *e=&E_1;
                  
                  //Declare Queue 1 pointer arrays
                  double *q_1=&Q_1_1;
                  double *p_q1=&P_Q1_1;
                  double *total_samples_for_PCR_Testing_q1=&total_samples_for_PCR_Testing_Q1_1;
                  double *total_samples_for_PCR_Testing_lag_1_q1=&total_samples_for_PCR_Testing_lag_1_Q1_1;
                  double *total_samples_for_PCR_Testing_lag_2_q1=&total_samples_for_PCR_Testing_lag_2_Q1_1;
                  double *y_q1=&Y_Q1_1;
                  
                  //Declare Queue NC pointer arrays
                  double *q_nc=&Q_NC_1;
                  double *total_samples_for_PCR_Testing_qnc=&total_samples_for_PCR_Testing_QNC_1;
                  double *total_samples_for_PCR_Testing_lag_1_qnc=&total_samples_for_PCR_Testing_lag_1_QNC_1;
                  double *total_samples_for_PCR_Testing_lag_2_qnc=&total_samples_for_PCR_Testing_lag_2_QNC_1;
                  double *y_qnc=&Y_QNC_1;
                  
                  //Declare Queue 3 pointer arrays
                  double *q_3=&Q_3_1;
                  double *p_q3=&P_Q3_1;
                  double *total_samples_for_PCR_Testing_q3=&total_samples_for_PCR_Testing_Q3_1;
                  double *total_samples_for_PCR_Testing_lag_1_q3=&total_samples_for_PCR_Testing_lag_1_Q3_1;
                  double *total_samples_for_PCR_Testing_lag_2_q3=&total_samples_for_PCR_Testing_lag_2_Q3_1;
                  double *y_q3=&Y_Q3_1;


                   //Error checks- Top of model

                  if(R_A < 0 || R_F < 0 || R_H < 0 || I_A < 0 || I_P < 0 || I_H < 0 || I_S_1 < 0 || I_S_2 < 0 || S < 0 || N < 0 || Backlog_Queue_1 < 0 || Q_2 < 0 ||Backlog_Queue_3 < 0 || Q_4 < 0 |Backlog_Queue_NC < 0|| First_re_test_Q1 < 0 || Second_re_test_Q1 < 0 || First_re_test_Q3 < 0 || Second_re_test_Q3 < 0 || total_samples_for_PCR_Testing_backlog_Q1 < 0 || total_samples_for_PCR_Testing_backlog_lag_1_Q1 < 0 || total_samples_for_PCR_Testing_backlog_lag_2_Q1 < 0 || total_samples_for_PCR_Testing_backlog_Q3 < 0 || total_samples_for_PCR_Testing_backlog_lag_1_Q3 < 0 || total_samples_for_PCR_Testing_backlog_lag_2_Q3 < 0 ||total_samples_for_PCR_Testing_backlog_QNC < 0 || total_samples_for_PCR_Testing_backlog_lag_1_QNC < 0 || total_samples_for_PCR_Testing_backlog_lag_2_QNC < 0 || neg_samples_Q1 < 0 || neg_samples_Q3 < 0 || total_neg_samples_all_queues < 0 ||  L_advanced_2_days < 0 || G_w_y < 0 ||  L_int < 0 || F_w_y < 0 ||  L_1 < 0 ||  w < 0 ||  L_2 < 0 ||  L_3 < 0 ||  L_4 < 0 || y < 0){
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at top of process model   t = %lg \\n\", t);

                      

                  }
                  
                  double sum_everything = R_A + R_F + R_H + I_P + I_A + I_H + I_S_1 + I_S_2 + S + N;
                  sum_everything = sum_everything + Backlog_Queue_1 + Backlog_Queue_3 + Backlog_Queue_NC + Q_2  + Q_4;
                  sum_everything = sum_everything + First_re_test_Q1 + First_re_test_Q3 +Second_re_test_Q1 + Second_re_test_Q3;
                  sum_everything = sum_everything + total_samples_for_PCR_Testing_backlog_Q1 + total_samples_for_PCR_Testing_backlog_Q3 + total_samples_for_PCR_Testing_backlog_QNC;
                  sum_everything = sum_everything + total_samples_for_PCR_Testing_backlog_lag_1_Q1 + total_samples_for_PCR_Testing_backlog_lag_1_Q3 + total_samples_for_PCR_Testing_backlog_lag_1_QNC;
                  sum_everything = sum_everything + total_samples_for_PCR_Testing_backlog_lag_2_Q1 + total_samples_for_PCR_Testing_backlog_lag_2_Q3 + total_samples_for_PCR_Testing_backlog_lag_2_QNC;
                  sum_everything = sum_everything + neg_samples_Q1 + neg_samples_Q3 + total_neg_samples_all_queues;
                  sum_everything = sum_everything + L_advanced_2_days + G_w_y + F_w_y + L_int + L_1 + L_2 + L_3 + L_4 + w + y;
                  if(isnan(sum_everything)){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at top of process model t = %lg \\n\", t);

                  }
                  
                  //Check Exposed Compartments
                  for(m=0; m<M; m++) {
                    if(isnan(e[m])){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at top of process model t = %lg \\n\", t);
                  
                    }
                    if(e[m] < 0) {
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at top of process model   t = %lg \\n\", t);
                    }
                  }
                  
                  //Check Q_1 and Q_NC arrays
                  for(v=0; v<V; v++) {
                    if(isnan(q_1[v]) || isnan(q_nc[v]) || isnan(p_q1[v]) || isnan(total_samples_for_PCR_Testing_q1[v]) || isnan(total_samples_for_PCR_Testing_qnc[v]) || isnan(total_samples_for_PCR_Testing_lag_1_q1[v]) || isnan(total_samples_for_PCR_Testing_lag_1_qnc[v]) || isnan(total_samples_for_PCR_Testing_lag_2_q1[v]) || isnan(total_samples_for_PCR_Testing_lag_2_qnc[v])){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at top of process model t = %lg \\n\", t);
                  
                    }
                    if(q_1[v] < 0 || q_nc[v] < 0 || p_q1[v] < 0 || total_samples_for_PCR_Testing_q1[v] < 0 || total_samples_for_PCR_Testing_qnc[v] < 0 || total_samples_for_PCR_Testing_lag_1_q1[v] < 0 || total_samples_for_PCR_Testing_lag_1_qnc[v] < 0 || total_samples_for_PCR_Testing_lag_2_q1[v] < 0 || total_samples_for_PCR_Testing_lag_2_qnc[v] < 0 ) {
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at top of process model   t = %lg \\n\", t);
                    }
                  }
                  
                  //Check Q_3 Arrays
                  for(k=0; k<K; k++) {
                    if(isnan(q_3[k]) || isnan(p_q3[k]) || isnan(total_samples_for_PCR_Testing_q3[k]) || isnan(total_samples_for_PCR_Testing_lag_1_q3[k]) || isnan(total_samples_for_PCR_Testing_lag_2_q3[k]) || isnan(q_3[k])){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at top of process model t = %lg \\n\", t);
                  
                    }
                    if(q_3[k] < 0 || p_q3[k] < 0 || total_samples_for_PCR_Testing_q3[k] < 0 || total_samples_for_PCR_Testing_lag_1_q3[k] < 0 || total_samples_for_PCR_Testing_lag_2_q3[k] < 0 ) {
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at top of process model   t = %lg \\n\", t);
                    }
                  }
                  
                  int print_out_top = (Error_Printing_Complete == FALSE) & (Neg_State_Value_Detected == TRUE || NAN_State_Value_Detected == TRUE);
                  if(print_out_top){
                      Rprintf(\"I_S_1 = %lg \\n\", I_S_1);
                      Rprintf(\"I_S_2 = %lg \\n\", I_S_2);
                      Rprintf(\"I_H = %lg \\n\", I_H);
                      Rprintf(\"I_P = %lg \\n\", I_P);
                      Rprintf(\"I_A = %lg \\n\", I_A);

                      Rprintf(\"Backlog_Queue_1 = %lg \\n\", Backlog_Queue_1);
                      Rprintf(\"Backlog_Queue_3 = %lg \\n\", Backlog_Queue_3);
                      Rprintf(\"Backlog_Queue_NC = %lg \\n\", Backlog_Queue_NC);
                      
                      Rprintf(\"First_re_test_Q1 = %lg \\n\", First_re_test_Q1);
                      Rprintf(\"First_re_test_Q3 = %lg \\n\", First_re_test_Q3);
                      Rprintf(\"Second_re_test_Q1 = %lg \\n\", Second_re_test_Q1);
                      Rprintf(\"Second_re_test_Q3 = %lg \\n\", Second_re_test_Q3);
                      
                      
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_Q1 = %lg \\n\", total_samples_for_PCR_Testing_backlog_Q1);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_1_Q1 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_1_Q1);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_2_Q1 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_2_Q1);
                      
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_Q3 = %lg \\n\", total_samples_for_PCR_Testing_backlog_Q3);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_1_Q3 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_1_Q3);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_2_Q3 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_2_Q3);
                      
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_QNC = %lg \\n\", total_samples_for_PCR_Testing_backlog_QNC);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_1_QNC = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_1_QNC);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_2_QNC = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_2_QNC);
                      
                      Rprintf(\"neg_samples_Q1 = %lg \\n\", neg_samples_Q1);
                      Rprintf(\"neg_samples_Q3 = %lg \\n\", neg_samples_Q3);
                      Rprintf(\"total_neg_samples_all_queues = %lg \\n\", total_neg_samples_all_queues);
                      
                      Rprintf(\"Q_2 = %lg \\n\", Q_2);
                      Rprintf(\"Q_4 = %lg \\n\", Q_4);
                      
                      Rprintf(\"F_w_y = %lg \\n\", F_w_y);
                      Rprintf(\"w = %lg \\n\", w);
                      Rprintf(\"y = %lg \\n\", y);
                      
                      Rprintf(\"E_1 = %lg \\n\", E_1);
                      Rprintf(\"E_2 = %lg \\n\", E_2);
                      Rprintf(\"E_3 = %lg \\n\", E_3);
                      Rprintf(\"E_4 = %lg \\n\", E_4);
                      Rprintf(\"E_5 = %lg \\n\", E_5);


                      Rprintf(\"R_H = %lg \\n\", R_H);

                      Rprintf(\"R_A = %lg \\n\", R_A);
                      Rprintf(\"R_F = %lg \\n\", R_F);


                      Rprintf(\"N = %lg \\n\", N);
                      Rprintf(\"S = %lg \\n\", S);
                      
                      Rprintf(\"L_advanced_2_days = %lg \\n\", L_advanced_2_days);
                      Rprintf(\"G_w_y = %lg \\n\", G_w_y);
                      Rprintf(\"L_int = %lg \\n\", L_int);
                      Rprintf(\"L_1 = %lg \\n\", L_1);
                      Rprintf(\"L_2 = %lg \\n\", L_2);
                      Rprintf(\"L_3 = %lg \\n\", L_3);
                      Rprintf(\"L_4 = %lg \\n\", L_4);
                      
                      Rprintf(\"Print out params  p_S = %lg \\n\", p_S);
                      Rprintf(\"p_H_cond_S = %lg \\n\", p_H_cond_S);
                      Rprintf(\"phi_E = %lg \\n\", phi_E);
                      Rprintf(\"phi_U = %lg \\n\", phi_U);
                      Rprintf(\"phi_S = %lg \\n\", phi_S);
                      Rprintf(\"h_V = %lg \\n\", h_V);
                      Rprintf(\"gamma = %lg \\n\", gamma);
                      Rprintf(\"R_0 = %lg \\n\", R_0);
                      Rprintf(\"b_q = %lg \\n\", b_q);
                      Rprintf(\"b_a = %lg \\n\", b_a);
                      Rprintf(\"b_p = %lg \\n\", b_p);
                      Rprintf(\"z_0 = %lg \\n\", z_0);
                      Rprintf(\"E_0 = %lg \\n\", E_0);
                      Rprintf(\"N_0 = %lg \\n\", N_0);
                      Rprintf(\"C_0 = %lg \\n\", C_0);
                      Rprintf(\"G_w_y_scaling = %lg \\n\", G_w_y_scaling);
                      
                      Rprintf(\"quarantine_start_time = %lg \\n\", quarantine_start_time);
                      Rprintf(\"PCR_sens = %lg \\n\", PCR_sens);
                      Rprintf(\"sigma_M = %lg \\n\", sigma_M);
                      
                      Rprintf(\"beta_w_3 = %lg \\n\", beta_w_3);
                      Rprintf(\"beta_w_2 = %lg \\n\", beta_w_2);
                      Rprintf(\"beta_w_1 = %lg \\n\", beta_w_1);
                      Rprintf(\"beta_w_0 = %lg \\n\", beta_w_0);
                      Rprintf(\"g_0 = %lg \\n\", g_0);
                      Rprintf(\"g_F = %lg \\n\", g_F);
                      Rprintf(\"sigma_epsilon = %lg \\n\", sigma_epsilon);
                      
                      //Print out exposed compartments
                      for(m=0; m<M; m++) {
                        Rprintf(\"m = %d \\n\", m);
                        Rprintf(\"e[m] = %lg \\n\", e[m]);
                      }
                      
                      //Print out Q_1 and Q_NC compartments
                      for(v=0; v<V; v++) {
                        Rprintf(\"v = %d \\n\", v);
                        Rprintf(\"q_1[v] = %lg \\n\", q_1[v]);
                        Rprintf(\"q_nc[v] = %lg \\n\", q_nc[v]);
                        Rprintf(\"p_q1[v] = %lg \\n\", p_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_q1[v] = %lg \\n\", total_samples_for_PCR_Testing_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_qnc[v] = %lg \\n\", total_samples_for_PCR_Testing_qnc[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_1_q1[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_1_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_1_qnc[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_1_qnc[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_2_q1[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_2_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_2_qnc[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_2_qnc[v]);
                      }
                      
                      //Print out Q_3 compartments
                      for(k=0; k<K; k++) {
                        Rprintf(\"k = %d \\n\", k);
                        Rprintf(\"q_3[k] = %lg \\n\", q_3[k]);
                        Rprintf(\"p_q3[k] = %lg \\n\", p_q3[k]);
                        Rprintf(\"total_samples_for_PCR_Testing_q3[k] = %lg \\n\", total_samples_for_PCR_Testing_q3[k]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_1_q3[k] = %lg \\n\", total_samples_for_PCR_Testing_lag_1_q3[k]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_2_q3[k] = %lg \\n\", total_samples_for_PCR_Testing_lag_2_q3[k]);
                      }
                      

                      Error_Printing_Complete = TRUE;
                  }
                  
                  
                   //Main Process Model Code Block
                   //Initialize transition array (make it NAN so it will trigger warnings if 
                   //used when unchanged)
                   for(m=0; m<M-1; m++) {
                        dE_m_E_m_1[m] = NAN;
                      }
                   
                   double total_time_infected = (1/gamma) + (1/phi_S);
                   double gamma_total = 1/total_time_infected;
                   double Beta_0 = R_0*(gamma_total);
                   double Beta_1 = b_q*Beta_0;
                   beta_t = Beta_0;
                   
                   if(t > quarantine_start_time){
                      beta_t = Beta_1;
                   }else{
                    if(t > social_distancing_start_time){
                      double m_q = (Beta_1 - Beta_0)/(quarantine_start_time - social_distancing_start_time);
                      beta_t = Beta_0 + m_q*(t-social_distancing_start_time);
                    }
                   }
                   
                   //Calculate transmssion rates in pre-symptomatic classes
                   double beta_p = b_p*beta_t;
                   
                   //Calculate transmssion rates in asymptomatic classes
                   double beta_a = b_a*beta_t;
                   
                   //Rates for Euler Multinom leaving Infected Pre-symptomatic class
                   double mu_I_P_I_A = (1-p_S)*phi_U;        //Moving from I_P to I_A
                   double mu_I_P_I_S_1 = p_S*phi_U;    //Moving from I_P to I_S_1 
                   
                   rate_I_P[0] = mu_I_P_I_A;
                   rate_I_P[1] = mu_I_P_I_S_1;
                   
                   //Rates for Euler Multinom leaving Infected Symptomatic class
                   double mu_I_S_1_I_H = p_H_cond_S*phi_S;            //Moving from I_S_1 to I_H
                   double mu_I_S_1_I_S_2 = (1-p_H_cond_S)*phi_S;        //Moving from I_S_1 to I_S_2
                   
                   rate_I_S_1[0] = mu_I_S_1_I_H;
                   rate_I_S_1[1] = mu_I_S_1_I_S_2;
                   
                   //Calculate force of infection
                   double lambda_FOI = (beta_t*(I_S_1 + I_S_2) + beta_a*(I_A) + beta_p*(I_P))/N;
                   double mu_S_E_1 = lambda_FOI;
                   
                   //Calculate rate of moving bewteen exposed compartments
                   double mu_E_m_E_m_1 = phi_E;
                   
                   //Calculate rate of moving from Infected Asymptomatic class to 
                   // Recovered Asymptomatic class
                   double mu_I_A_R_A = phi_S;
                   
                   
                   //Calculate rate of moving from last Exposed compartment to 
                   //Infected Pre-sympatomatic class
                   double mu_E_M_I_P = phi_E;

                   //Calculate rate of moving from Infected Flu-like class to Recovered Flu-Like Class
                   double mu_I_S_2_R_F = gamma;
                   
                   //Calculate rate of moving from Infected Hospitalized class to Recovered Hospitalized Class
                   double mu_I_H_R_H = h_V;
                   
                   //Binomial transitions (out of S, I_S_2, E_M, I_A, and I_H, and within E)
                   double dSE_1 = rbinom(S, 1 - exp(-mu_S_E_1*dt));
                   
                   double dI_A_R_A = rbinom(I_A, 1 - exp(-mu_I_A_R_A*dt));
                   
                   double dE_M_I_P = rbinom(e[M-1], 1 - exp(-mu_E_M_I_P*dt));
                   
                   double dI_S_2_R_F = rbinom(I_S_2, 1 - exp(-mu_I_S_2_R_F*dt));
                   double dI_H_R_H = rbinom(I_H, 1 - exp(-mu_I_H_R_H*dt));
                   
                   //Exposed class binomial transitions (start from second compartment)
                   for(m = 0; m < M-1; m++){
                    dE_m_E_m_1[m] = rbinom(e[m], 1 - exp(-mu_E_m_E_m_1*dt));
                   }
                   
                   //Euler multinomial transitions
                   
                   //Infected Pre-symptomatic compartment
                   reulermultinom(2,I_P, &rate_I_P[0], dt, &trans_I_P[0]);
                   
                   //Infected Symptomatic comparment
                   reulermultinom(2,I_S_1, &rate_I_S_1[0], dt, &trans_I_S_1[0]);
                   
                   //Get compartment transitions from Euler multinom output
                   
                   //Infected Pre-symptomatic Compartment
                   double dI_P_I_A = trans_I_P[0];
                   double dI_P_I_S_1 = trans_I_P[1];
                   
                   //Infected Symptomatic Compartment
                   double dI_S_1_I_H = trans_I_S_1[0];
                   double dI_S_1_I_S_2 = trans_I_S_1[1];
                   
                   //Update state variables using transition increments
                   
                   //Susceptible Compartment
                   S += -dSE_1;
                   
                   //Exposed Compartments
                   //First compartment
                   e[0] += dSE_1 - dE_m_E_m_1[0];
                   
                   //Inner compartments (start from second compartment,
                   // end with second to last comparment)
                   for(m = 1; m < M-1; m++){
                    e[m] += dE_m_E_m_1[m-1] - dE_m_E_m_1[m];
                   }
                   
                   //Outermost exposed compartment
                   e[M-1] += dE_m_E_m_1[M-2]  - dE_M_I_P;
                   
                   
                   
                   //Infected Compartments
                   
                   //Pre-Symptomatic Infected Compartment
                   I_P += dE_M_I_P- dI_P_I_A - dI_P_I_S_1;
                   
                   I_S_1 += dI_P_I_S_1 - dI_S_1_I_H - dI_S_1_I_S_2;
                   
                   //Hospitalized Compartment
                   I_H += dI_S_1_I_H - dI_H_R_H;
                   
                   //Flu-Like Infections Compartment
                   I_S_2 += dI_S_1_I_S_2 - dI_S_2_R_F;
                   
                   //Recovered Compartments
                   I_A += dI_P_I_A - dI_A_R_A;
                   R_A += dI_A_R_A;
                   R_F += dI_S_2_R_F;
                   R_H += dI_H_R_H;
                   
                   //Total Population (does not change)
                   N += 0;
                   
                   //Reported Cases
                   C_Q1 += dI_S_1_I_H;  //Entering Queue 1
                   C_Q3 += dI_S_1_I_S_2;  //Entering Queue 2
                   
                   //Read in testing info
                   L_int = nearbyint(L_advanced_2_days);
                   
                   double epsilon = rnorm(0,sigma_epsilon);
                   
                   //Calculate estimated non-COVID respiratory infections
                   G_w_y = g_0 + g_F*F_w_y + beta_w_3*(w*w*w) + beta_w_2*(w*w) + beta_w_1*w + beta_w_0 + epsilon;
                   if(t > quarantine_start_time){
                    G_w_y = G_w_y_scaling*G_w_y;
                   }else{
                    G_w_y = G_w_y_scaling*G_w_y;
                   }
                   
                   
                   C_QNC = nearbyint(G_w_y);
                   
                   //Queue 1
                   
                   //Initial states
                   neg_samples_Q1 = 0;
                   
                   //Add new cases to Q1
                   q_1[0] = q_1[0] + C_Q1;
                   
                   //Add new cases to QNC
                   q_nc[0] = q_nc[0] + C_QNC;
                   
                   //Determine number of samples that will tested from the backlog
                   // of Queue 1
                   double total_samples_for_PCR_Testing_backlog_always_test = Backlog_Queue_1 + Backlog_Queue_NC;
                   
                   //If backlog is greater than testing capacity
                   if(total_samples_for_PCR_Testing_backlog_always_test > L_int){
                      total_samples_for_PCR_Testing_backlog_Q1 = rhyper(Backlog_Queue_1, Backlog_Queue_NC, L_int);
                      total_samples_for_PCR_Testing_backlog_QNC = L_int - total_samples_for_PCR_Testing_backlog_Q1;
                      Backlog_Queue_1 = Backlog_Queue_1 - total_samples_for_PCR_Testing_backlog_Q1;
                      Backlog_Queue_NC = Backlog_Queue_NC - total_samples_for_PCR_Testing_backlog_QNC;
                      L_1 = 0;
                      
                    //Else (if the testing capacity is greater than the backlog)
                   }else{
                      total_samples_for_PCR_Testing_backlog_Q1 = Backlog_Queue_1;
                      total_samples_for_PCR_Testing_backlog_QNC = Backlog_Queue_NC;
                      
                      L_1 = L_int - total_samples_for_PCR_Testing_backlog_Q1 - total_samples_for_PCR_Testing_backlog_QNC;
                      
                      Backlog_Queue_1 = 0;
                      Backlog_Queue_NC = 0;
                   }
                   
                   //Simulate PCR for backlogged cases in Q1
                   Y_Q1_backlog = rbinom(total_samples_for_PCR_Testing_backlog_lag_2_Q1, PCR_sens);
                   Y_QNC_backlog = total_samples_for_PCR_Testing_backlog_lag_2_QNC;
                   neg_samples_Q1 = neg_samples_Q1 + total_samples_for_PCR_Testing_backlog_lag_2_Q1 - Y_Q1_backlog;
                   
                   //Update PCR testing compartments for Q1 backlog
                   total_samples_for_PCR_Testing_backlog_lag_2_Q1 = total_samples_for_PCR_Testing_backlog_lag_1_Q1;
                   total_samples_for_PCR_Testing_backlog_lag_2_QNC = total_samples_for_PCR_Testing_backlog_lag_1_QNC;
                   total_samples_for_PCR_Testing_backlog_lag_1_Q1 = total_samples_for_PCR_Testing_backlog_Q1;
                   total_samples_for_PCR_Testing_backlog_lag_1_QNC = total_samples_for_PCR_Testing_backlog_QNC;
                   
                   //Determine number of samples that will be tested from each sampling cohort of Queue 1
                   
                   //Loop through each cohort in Queue 1 q_1[v] (and Queue NC q_nc[v]) 
                   // starting with the oldest (q_1[V-1]/q_nc[V-1])
                   // and ending with the most recent q_1[0]/q_nc[1].
                   for(v=V-1; v>=0; v--) {
                   
                      //If L_1 is smaller than cohort
                      if(L_1 < (q_1[v] + q_nc[v])){
                        total_samples_for_PCR_Testing_q1[v] = rhyper(q_1[v], q_nc[v], L_1);
                        total_samples_for_PCR_Testing_qnc[v] = L_1 - total_samples_for_PCR_Testing_q1[v];
                        q_1[v] = q_1[v] - total_samples_for_PCR_Testing_q1[v];
                        q_nc[v] = q_nc[v] - total_samples_for_PCR_Testing_qnc[v];
                        L_1 = 0;
                       
                       //Else there is enough capacity to test cohort v
                       // in Q1/QNC 
                      }else{
                        total_samples_for_PCR_Testing_q1[v] = q_1[v];
                        total_samples_for_PCR_Testing_qnc[v] = q_nc[v];
                        q_1[v] = 0;
                        q_nc[v] = 0;
                        L_1 = L_1 - total_samples_for_PCR_Testing_q1[v] - total_samples_for_PCR_Testing_qnc[v];
                      }
                   }
                   
                   
                   //Simulate PCR Testing on each sampling cohort in Q1
                   //Loop over all cohorts v from 1:V
                   for(v=0; v<V; v++) {
                    y_q1[v] = rbinom(total_samples_for_PCR_Testing_lag_2_q1[v], PCR_sens);
                    y_qnc[v] = total_samples_for_PCR_Testing_lag_2_qnc[v];
                    neg_samples_Q1 = neg_samples_Q1 + total_samples_for_PCR_Testing_lag_2_q1[v] - y_q1[v];
                   }
                   
                   //Update lags for total samples for PCR testing for Q1
                   //Loop over all cohorts v from 1:V
                   for(v=0; v<V; v++) {
                      total_samples_for_PCR_Testing_lag_2_q1[v] = total_samples_for_PCR_Testing_lag_1_q1[v];
                      total_samples_for_PCR_Testing_lag_2_qnc[v] = total_samples_for_PCR_Testing_lag_1_qnc[v];
                      total_samples_for_PCR_Testing_lag_1_q1[v] = total_samples_for_PCR_Testing_q1[v];
                      total_samples_for_PCR_Testing_lag_1_qnc[v] = total_samples_for_PCR_Testing_qnc[v]; 
                   }
                   
                   //Update positive case matrix (p_Q1)
                   //For all daily sampling cohorts v within the last 14 days:
                   for(v=0; v<V; v++) {
                    p_q1[v] = p_q1[v] + y_q1[v];
                   }
                   
                   //Re-testing (Add lagged positive samples from Queue 1 into Queue 2)
                   //The state variable First_re_test_Q1 are samples that were 
                   // first re-sampled during the previous day. 
                   // They now need to be re-sampled a second time.
                   Second_re_test_Q1 = First_re_test_Q1;
                   C_Q2 = Second_re_test_Q1;
                   
                   //Let V-1 be the oldest cohort stored (V=13). 
                   // This cohort will have their first re-sampling conducted.
                   First_re_test_Q1 = p_q1[V-1];
                   C_Q2 = C_Q2 + First_re_test_Q1;
                   
                   //Increment P_Q1 Sampling Cohorts
                   P_Q1_old[0] = p_q1[0]; //Store first cohort 
                   //For v in 2:V: (or in C notation v in 1:V-1):
                   for(v=1; v<V; v++) {
                    P_Q1_old[v] = p_q1[v];
                    p_q1[v] = P_Q1_old[v-1];
                   }
                   
                   //For the newest cohort where v = 1 (or 0 in C notation):
                   p_q1[0] = 0; //(Making space for next cohort arrival)
                   
                   //Note that the oldest cohort (p_q1_old[V-1]) is never used
                   //since p_q1[V-1] has already been transferred to First_re_test_Q1
                   
                   //Create placeholder arrays to store
                   // current queues
                   for(v=0; v<V; v++) {
                    Q1_old[v] = q_1[v];
                    QNC_old[v] = q_nc[v];
                   }
                   
                   //Add oldest cohort to backlog 
                   // v=V (or V-1 in C notation)
                   Backlog_Queue_1 = Backlog_Queue_1 + Q1_old[V-1];
                   Backlog_Queue_NC = Backlog_Queue_NC + QNC_old[V-1];
                   
                   //Update Q1 Sampling cohorts by 1
                   //For integer v in v>1 and v<=V 
                   // (or in C notation v = 1 to v <V):
                   for(v=1; v<V; v++) {
                    q_1[v] = Q1_old[v-1];
                    q_nc[v] = QNC_old[v-1];
                   }
                   
                   //Make space for newest cohort when v=1 
                   // (in C notation v = 0) :
                   q_1[0] = 0;
                   q_nc[0] = 0;
                   
                   //Queue 2
                   
                   //Add new cases to queue
                   Q_2 = Q_2+C_Q2;
                   
                   //Determine number of samples that will be tested from Queue 2
                   //Recall that L2 is the testing capacity available at the end of Queue 2.
                   L_2 = L_1;
                   
                   //If there is not enough testing capacity to test all of Queue 2
                   if (Q_2 > L_2){
                     total_samples_for_PCR_Testing_Q2 = L_2;
                     Q_2 = Q_2 - L_2;
                     L_2 = 0;
                    //Else if there is enough capcity to test all of Queue 2
                   }else{
                      total_samples_for_PCR_Testing_Q2 = Q_2;
                      Q_2 = 0;
                      L_2 = L_2 - total_samples_for_PCR_Testing_Q2;
                   }
                   
                   //Recall that we do not keep track of the results of the PCR testing in Queue 2, 
                   //as it will not impact the count of reported cases. 
                   //We are also not worried about lags here.
                   
                   //Queue 3
                   
                   //Initial states
                   neg_samples_Q3 = 0;
                   L_3 = L_2;
                   
                   //Take into account loss rate due to recovery
                   //NOT IMPLEMENTED YET
                   double mu_2 = gamma;
                   double single_cohort_loss = 0;
                   for(k=0; k<K; k++) {
                    single_cohort_loss = rbinom(q_3[k], 1 - exp(-mu_2*dt));
                    q_3[k] = q_3[k] - single_cohort_loss;
                   }
                   
                   double backlog_loss = rbinom(Backlog_Queue_3, 1 - exp(-mu_2*dt));
                   Backlog_Queue_3 = Backlog_Queue_3 - backlog_loss;
                   
                   //Simulate loss in Q3 backlog
                   
                   //Add new cases to Q3
                   q_3[0] = q_3[0] + C_Q3; 
                   
                   //Determine number of samples that will be tested from the backlog of Queue 3
                   //If Backlog is greater than testing capacity:
                   if(Backlog_Queue_3 > L_3){
                    total_samples_for_PCR_Testing_backlog_Q3 = L_3;
                    Backlog_Queue_3 = Backlog_Queue_3 - L_3;
                    L_3 = 0;
                    
                    //There is enough capcity to test the whole Q_3 backlog
                   }else{
                    total_samples_for_PCR_Testing_backlog_Q3 = Backlog_Queue_3;
                    L_3 = L_3 - total_samples_for_PCR_Testing_backlog_Q3;
                    Backlog_Queue_3 = 0;
                   }
                   
                   //Simulate PCR for backlogged cases in Q3
                   Y_Q3_backlog = rbinom(total_samples_for_PCR_Testing_backlog_lag_2_Q3, PCR_sens);
                   neg_samples_Q3 = neg_samples_Q3 + total_samples_for_PCR_Testing_backlog_lag_2_Q3 - Y_Q3_backlog;
                   
                   //Update PCR testing compartments for Q3 backlog
                   total_samples_for_PCR_Testing_backlog_lag_2_Q3 = total_samples_for_PCR_Testing_backlog_lag_1_Q3;
                   total_samples_for_PCR_Testing_backlog_lag_1_Q3 = total_samples_for_PCR_Testing_backlog_Q3;
                   
                   //Determine number of samples that will be tested from each sampling cohort of Queue 3
                   //Loop through each cohort in Queue 3 (Q3_k) starting with the oldest (Q3_K) 
                   // and ending with the most recent (Q3_1).
                   for(k=K-1; k>=0; k--) {
                    //If L_3 is smaller than cohort (i.e.  L_3 < Q_{3_k} :
                    if(q_3[k] > L_3){
                      total_samples_for_PCR_Testing_q3[k] = L_3;
                      q_3[k] = q_3[k] - L_3;
                      L_3 = 0;
                      
                      //There is enough capacity to test cohort $k$
                    }else{
                      total_samples_for_PCR_Testing_q3[k] = q_3[k];
                      q_3[k] = 0;
                      L_3 = L_3 - total_samples_for_PCR_Testing_q3[k];
                    }
                   }
                   
                   //Simulate PCR Testing on each sampling cohort in Q3
                   //Loop over all cohorts k from 1:K
                   for(k=0; k<K; k++) {
                    y_q3[k] = rbinom(total_samples_for_PCR_Testing_lag_2_q3[k], PCR_sens);
                    neg_samples_Q3 = neg_samples_Q3 + total_samples_for_PCR_Testing_lag_2_q3[k] - y_q3[k];
                   }
                   
                   //Update lags for total samples for PCR testing for Q3
                   for(k=0; k<K; k++) {
                    total_samples_for_PCR_Testing_lag_2_q3[k] = total_samples_for_PCR_Testing_lag_1_q3[k];
                    total_samples_for_PCR_Testing_lag_1_q3[k] = total_samples_for_PCR_Testing_q3[k];
                   }
                   
                   //Update positive case matrix (P_Q3)
                   //For all daily sampling cohorts k within the last 14 days:
                   for(k=0; k<K; k++) {
                    p_q3[k] = p_q3[k] + y_q3[k];
                   }
                   
                   //Re-testing (Add lagged positive samples from Queue 3 into Queue 4)
                   //The state variable First_re_test_Q3 are samples that we first re-sampled
                   // during the previous day. They now need to be re-sampled a second time.
                   Second_re_test_Q3 = First_re_test_Q3;
                   C_Q4 = Second_re_test_Q3;
                   
                   //Let K be the oldest cohort stored (K=14).
                   //This cohort will have their first re-sampling conducted.
                   // We use K-1 for C notation.
                   First_re_test_Q3 = p_q3[K-1];
                   C_Q4 = C_Q4 + First_re_test_Q3;
                     
                   //Increment P_Q3 Sampling Cohorts
                   P_Q3_old[0] = p_q3[0]; //Store first cohort 
                   //For k in 2:K: (or in C notation k in 1:K-1):
                   for(k=1; k<K; k++) {
                    P_Q3_old[k] = p_q3[k];
                    p_q3[k] = P_Q3_old[k-1];
                   }
                   
                   //For the newest cohort where k = 1 (or k=0 in C notation):
                   p_q3[0] = 0; //(Making space for next cohort arrival)
                   
                   //Note that the oldest cohort (p_q3_old[K-1]) is never used
                   //since p_q3[K-1] has already been transferred to First_re_test_Q3
                   
                   //Create placeholder array to store
                   // current Q3
                   for(k=0; k<K; k++) {
                    Q3_old[k] = q_3[k];
                   }
                   
                   //Add oldest cohort to backlog 
                   // k=K (or K-1 in C notation)
                   Backlog_Queue_3 = Backlog_Queue_3 + Q3_old[K-1];
                   
                   //Update Q3 Sampling cohorts by 1
                   //For integer k in k>1 and k<=K 
                   // (or in C notation k = 1 to k < K):
                   for(k=1; k<K; k++) {
                    q_3[k] = Q3_old[k-1];
                   }
                   
                   //Make space for newest cohort when k=1 
                   // (in C notation k = 0) :
                   q_3[0] = 0;
                   
                   //Queue 4
                   //Add new cases to queue
                   Q_4 = Q_4 + C_Q4;
                   
                   //Determine number of samples that will be tested from Queue 4
                   //Recall that L_4 is the testing capacity available at the end of Queue 4.
                   L_4 = L_3;
                   
                   //If there is insufficent capacity (Queue 4 is bigger than L_4)
                   if(Q_4>L_4){
                    total_samples_for_PCR_Testing_Q4 = L_4;
                    Q_4 = Q_4 - L_4;
                    L_4 = 0;
                    
                    //There is enough capacity(Q_4 < L_4)
                   }else{
                    total_samples_for_PCR_Testing_Q4 = Q_4;
                    Q_4 = 0;
                    L_4 = L_4 - total_samples_for_PCR_Testing_Q4;
                   }
                   
                   //Recall that we do not keep track of the results of the PCR testing in Queue 4,
                   //as it will not impact the count of reported cases. 
                   //We are also not worried about lags here,
                   //and we assume that all samples that enter Queue 4
                   //are eventually tested (no additional loss rates).
                   
                   //Queue 5:Asymptomatic Testing
                   
                   
                   total_sample_size = S + E_1 + E_2 + E_3 + E_4 + E_5 + I_P + I_S_1 + R_A + R_F + R_H;
                   multinom_prob[0] = PCR_sens*E_1/total_sample_size;
                   multinom_prob[1] = PCR_sens*E_2/total_sample_size;
                   multinom_prob[2] = PCR_sens*E_3/total_sample_size;
                   multinom_prob[3] = PCR_sens*E_4/total_sample_size;
                   multinom_prob[4] = PCR_sens*E_5/total_sample_size;
                   multinom_prob[5] = PCR_sens*I_P/total_sample_size;
                   multinom_prob[6] = PCR_sens*I_S_1/total_sample_size;
                   double total_Q5_prob = multinom_prob[0] +  multinom_prob[1] +  multinom_prob[2];
                   total_Q5_prob = total_Q5_prob +  multinom_prob[3] +  multinom_prob[4] +  multinom_prob[5];
                   total_Q5_prob = total_Q5_prob +  multinom_prob[6] ;
                   multinom_prob[7] = 1 -  total_Q5_prob;
                   //rmultinom(1,L_4,  &multinom_prob, &multinom_output);
                   rmultinom(L_4, &multinom_prob[0], 8,  &multinom_output[0]);
                   double E_1_infected = 0;
                   double E_2_infected = 0;
                   double E_3_infected = 0;
                   double E_4_infected = 0;
                   double E_5_infected = 0;
                   double I_P_infected = 0;
                   double I_S_1_infected = 0;
                   neg_samples_Q5 = 0;
                   
                   //Update population comparmtents 
                   E_1 = E_1 - E_1_infected;
                   E_2 = E_2 - E_2_infected;
                   E_3 = E_3 - E_3_infected;
                   E_4 = E_4 - E_4_infected;
                   E_5 = E_5 - E_5_infected;
                   I_P = I_P - I_P_infected;
                   I_S_1 = I_S_1 - I_S_1_infected;
                   
                   Y_Q5 = E_1_infected + E_2_infected + E_3_infected + E_4_infected + E_5_infected +I_P_infected + I_S_1_infected;
                   
                   A_T = A_T + Y_Q5;

                   //Calculate total cases to report
                   total_neg_samples_all_queues = 0;
                   for(v=0; v<V; v++) {
                    Y_sum = Y_sum + y_q1[v];
                    total_neg_samples_all_queues = total_neg_samples_all_queues +  y_qnc[v];
                  
                   }
                   for(k=0; k<K; k++) {
                    Y_sum = Y_sum + y_q3[k];

                   }
                   
                   Y_sum = Y_sum + Y_Q1_backlog + Y_Q3_backlog;
                   total_neg_samples_all_queues = total_neg_samples_all_queues +neg_samples_Q1;
                   total_neg_samples_all_queues = total_neg_samples_all_queues +neg_samples_Q3;
                   
                   //double positive_plus_negative_tests = Y_sum + total_neg_samples_all_queues;
                   
                   
                   
                   //Prop_Positive_Tests_Track = Y_sum/positive_plus_negative_tests;
                   //Toy Reporting (Queues not yet implemented)
                   //Y_Q1 = C_Q1;
                   //Y_Q3 = C_Q3;
                   //Y_sum = C_Q1 + C_Q3;
                   if(L_orig >0){
                      Prop_Positive_Tests_Track = Y_sum/L_orig;
                   }else{
                      Prop_Positive_Tests_Track = 0; //Assign 0 if no testing yet
                   }
                      

                  
                   
                   //Error checks- Bottom of model

                  if(R_A < 0 || R_F < 0 || R_H < 0 || I_A < 0 || I_P < 0 ||  I_H < 0 || I_S_1 < 0 || I_S_2 < 0 || S < 0 || N < 0 || Backlog_Queue_1 < 0 || Q_2 < 0 ||Backlog_Queue_3 < 0 || Q_4 < 0 |Backlog_Queue_NC < 0|| First_re_test_Q1 < 0 || Second_re_test_Q1 < 0 || First_re_test_Q3 < 0 || Second_re_test_Q3 < 0 || total_samples_for_PCR_Testing_backlog_Q1 < 0 || total_samples_for_PCR_Testing_backlog_lag_1_Q1 < 0 || total_samples_for_PCR_Testing_backlog_lag_2_Q1 < 0 || total_samples_for_PCR_Testing_backlog_Q3 < 0 || total_samples_for_PCR_Testing_backlog_lag_1_Q3 < 0 || total_samples_for_PCR_Testing_backlog_lag_2_Q3 < 0 ||total_samples_for_PCR_Testing_backlog_QNC < 0 || total_samples_for_PCR_Testing_backlog_lag_1_QNC < 0 || total_samples_for_PCR_Testing_backlog_lag_2_QNC < 0 || neg_samples_Q1 < 0 || neg_samples_Q3 < 0 || total_neg_samples_all_queues < 0 ||  L_advanced_2_days < 0 || G_w_y < 0 ||  L_int < 0 || F_w_y < 0 ||  L_1 < 0 ||  w < 0 ||  L_2 < 0 ||  L_3 < 0 ||  L_4 < 0 || y < 0){
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at bottom of process model   t = %lg \\n\", t);

                      

                  }
                  
                  sum_everything = R_A + R_F + R_H + I_A + I_P + I_H + I_S_1 + I_S_2 + S + N;
                  sum_everything = sum_everything + Backlog_Queue_1 + Backlog_Queue_3 + Backlog_Queue_NC + Q_2  + Q_4;
                  sum_everything = sum_everything + First_re_test_Q1 + First_re_test_Q3 +Second_re_test_Q1 + Second_re_test_Q3;
                  sum_everything = sum_everything + total_samples_for_PCR_Testing_backlog_Q1 + total_samples_for_PCR_Testing_backlog_Q3 + total_samples_for_PCR_Testing_backlog_QNC;
                  sum_everything = sum_everything + total_samples_for_PCR_Testing_backlog_lag_1_Q1 + total_samples_for_PCR_Testing_backlog_lag_1_Q3 + total_samples_for_PCR_Testing_backlog_lag_1_QNC;
                  sum_everything = sum_everything + total_samples_for_PCR_Testing_backlog_lag_2_Q1 + total_samples_for_PCR_Testing_backlog_lag_2_Q3 + total_samples_for_PCR_Testing_backlog_lag_2_QNC;
                  sum_everything = sum_everything + neg_samples_Q1 + neg_samples_Q3 + total_neg_samples_all_queues;
                  sum_everything = sum_everything + L_advanced_2_days + G_w_y + F_w_y + L_int + L_1 + L_2 + L_3 + L_4 + w + y;
                  if(isnan(sum_everything)){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at bottom of process model t = %lg \\n\", t);

                  }
                  
                  //Check Exposed Compartments
                  for(m=0; m<M; m++) {
                    if(isnan(e[m])){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at bottom of process model t = %lg \\n\", t);
                  
                    }
                    if(e[m] < 0) {
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at bottom of process model   t = %lg \\n\", t);
                    }
                  }
                  
                  //Check Q_1 and Q_NC arrays
                  for(v=0; v<V; v++) {
                    if(isnan(q_1[v]) || isnan(q_nc[v]) || isnan(p_q1[v]) || isnan(total_samples_for_PCR_Testing_q1[v]) || isnan(total_samples_for_PCR_Testing_qnc[v]) || isnan(total_samples_for_PCR_Testing_lag_1_q1[v]) || isnan(total_samples_for_PCR_Testing_lag_1_qnc[v]) || isnan(total_samples_for_PCR_Testing_lag_2_q1[v]) || isnan(total_samples_for_PCR_Testing_lag_2_qnc[v])){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at bottom of process model t = %lg \\n\", t);
                  
                    }
                    if(q_1[v] < 0 || q_nc[v] < 0 || p_q1[v] < 0 || total_samples_for_PCR_Testing_q1[v] < 0 || total_samples_for_PCR_Testing_qnc[v] < 0 || total_samples_for_PCR_Testing_lag_1_q1[v] < 0 || total_samples_for_PCR_Testing_lag_1_qnc[v] < 0 || total_samples_for_PCR_Testing_lag_2_q1[v] < 0 || total_samples_for_PCR_Testing_lag_2_qnc[v] < 0 ) {
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at bottom of process model   t = %lg \\n\", t);
                    }
                  }
                  
                  //Check Q_3 Arrays
                  for(k=0; k<K; k++) {
                    if(isnan(q_3[k]) || isnan(p_q3[k]) || isnan(total_samples_for_PCR_Testing_q3[k]) || isnan(total_samples_for_PCR_Testing_lag_1_q3[k]) || isnan(total_samples_for_PCR_Testing_lag_2_q3[k]) || isnan(q_3[k])){
                      NAN_State_Value_Detected = TRUE;
                      Rprintf(\"nan state variable detected at bottom of process model t = %lg \\n\", t);
                  
                    }
                    if(q_3[k] < 0 || p_q3[k] < 0 || total_samples_for_PCR_Testing_q3[k] < 0 || total_samples_for_PCR_Testing_lag_1_q3[k] < 0 || total_samples_for_PCR_Testing_lag_2_q3[k] < 0 ) {
                      Neg_State_Value_Detected = TRUE;
                      Rprintf(\"Negative state variable detected at bottom of process model   t = %lg \\n\", t);
                    }
                  }
                  
                  int print_out_bottom = (Error_Printing_Complete == FALSE) & (Neg_State_Value_Detected == TRUE || NAN_State_Value_Detected == TRUE);
                  if(print_out_bottom){
                      Rprintf(\"I_S_1 = %lg \\n\", I_S_1);
                      Rprintf(\"I_S_2 = %lg \\n\", I_S_2);
                      Rprintf(\"I_H = %lg \\n\", I_H);
                      Rprintf(\"I_A = %lg \\n\", I_A);
                      Rprintf(\"I_P = %lg \\n\", I_P);

                      Rprintf(\"Backlog_Queue_1 = %lg \\n\", Backlog_Queue_1);
                      Rprintf(\"Backlog_Queue_3 = %lg \\n\", Backlog_Queue_3);
                      Rprintf(\"Backlog_Queue_NC = %lg \\n\", Backlog_Queue_NC);
                      
                      Rprintf(\"First_re_test_Q1 = %lg \\n\", First_re_test_Q1);
                      Rprintf(\"First_re_test_Q3 = %lg \\n\", First_re_test_Q3);
                      Rprintf(\"Second_re_test_Q1 = %lg \\n\", Second_re_test_Q1);
                      Rprintf(\"Second_re_test_Q3 = %lg \\n\", Second_re_test_Q3);
                      
                      
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_Q1 = %lg \\n\", total_samples_for_PCR_Testing_backlog_Q1);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_1_Q1 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_1_Q1);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_2_Q1 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_2_Q1);
                      
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_Q3 = %lg \\n\", total_samples_for_PCR_Testing_backlog_Q3);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_1_Q3 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_1_Q3);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_2_Q3 = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_2_Q3);
                      
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_QNC = %lg \\n\", total_samples_for_PCR_Testing_backlog_QNC);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_1_QNC = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_1_QNC);
                      Rprintf(\"total_samples_for_PCR_Testing_backlog_lag_2_QNC = %lg \\n\", total_samples_for_PCR_Testing_backlog_lag_2_QNC);
                      
                      Rprintf(\"neg_samples_Q1 = %lg \\n\", neg_samples_Q1);
                      Rprintf(\"neg_samples_Q3 = %lg \\n\", neg_samples_Q3);
                      Rprintf(\"total_neg_samples_all_queues = %lg \\n\", total_neg_samples_all_queues);
                      
                      Rprintf(\"Q_2 = %lg \\n\", Q_2);
                      Rprintf(\"Q_4 = %lg \\n\", Q_4);
                      
                      Rprintf(\"F_w_y = %lg \\n\", F_w_y);
                      Rprintf(\"w = %lg \\n\", w);
                      Rprintf(\"y = %lg \\n\", y);
                      
                      Rprintf(\"E_1 = %lg \\n\", E_1);
                      Rprintf(\"E_2 = %lg \\n\", E_2);
                      Rprintf(\"E_3 = %lg \\n\", E_3);
                      Rprintf(\"E_4 = %lg \\n\", E_4);
                      Rprintf(\"E_5 = %lg \\n\", E_5);
                      Rprintf(\"I_P = %lg \\n\", I_P);
                      
                      Rprintf(\"dE_m_E_m_1[M-2] = %lg \\n\", dE_m_E_m_1[M-2]);
                      Rprintf(\"lambda_FOI = %lg \\n\", lambda_FOI);
                      Rprintf(\"E_1_infected = %lg \\n\", E_1_infected);
                      Rprintf(\"dSE_1 = %lg \\n\", dSE_1);
                      Rprintf(\"dE_m_E_m_1[0] = %lg \\n\", dE_m_E_m_1[0]);
                      Rprintf(\"beta_t = %lg \\n\", beta_t);
                      Rprintf(\"beta_a = %lg \\n\", beta_a);
                      Rprintf(\"dE_m_E_m_1[0] = %lg \\n\", dE_m_E_m_1[0]);
                      Rprintf(\"dE_m_E_m_1[1] = %lg \\n\", dE_m_E_m_1[1]);
                      Rprintf(\"dE_m_E_m_1[2] = %lg \\n\", dE_m_E_m_1[2]);
                      
                      Rprintf(\"dI_P_I_A = %lg \\n\", dI_P_I_A);
                      Rprintf(\"dI_P_I_S_1 = %lg \\n\", dI_P_I_S_1);
                      Rprintf(\"dE_M_I_P = %lg \\n\", dE_M_I_P);
                      
                      

                      Rprintf(\"R_H = %lg \\n\", R_H);

                      Rprintf(\"R_A = %lg \\n\", R_A);
                      Rprintf(\"R_F = %lg \\n\", R_F);


                      Rprintf(\"N = %lg \\n\", N);
                      Rprintf(\"S = %lg \\n\", S);
                      
                      Rprintf(\"L_advanced_2_days = %lg \\n\", L_advanced_2_days);
                      Rprintf(\"G_w_y = %lg \\n\", G_w_y);
                      Rprintf(\"L_int = %lg \\n\", L_int);
                      Rprintf(\"L_1 = %lg \\n\", L_1);
                      Rprintf(\"L_2 = %lg \\n\", L_2);
                      Rprintf(\"L_3 = %lg \\n\", L_3);
                      Rprintf(\"L_4 = %lg \\n\", L_4);
                      
                      Rprintf(\"Print out params  p_S = %lg \\n\", p_S);
                      Rprintf(\"p_H_cond_S = %lg \\n\", p_H_cond_S);
                      Rprintf(\"phi_E = %lg \\n\", phi_E);
                      Rprintf(\"phi_U = %lg \\n\", phi_U);
                      Rprintf(\"phi_S = %lg \\n\", phi_S);
                      Rprintf(\"h_V = %lg \\n\", h_V);
                      Rprintf(\"gamma = %lg \\n\", gamma);
                      Rprintf(\"R_0 = %lg \\n\", R_0);
                      Rprintf(\"b_q = %lg \\n\", b_q);
                      Rprintf(\"b_a = %lg \\n\", b_a);
                      Rprintf(\"b_p = %lg \\n\", b_p);
                      Rprintf(\"z_0 = %lg \\n\", z_0);
                      Rprintf(\"E_0 = %lg \\n\", E_0);
                      Rprintf(\"N_0 = %lg \\n\", N_0);
                      Rprintf(\"C_0 = %lg \\n\", C_0);
                      Rprintf(\"G_w_y_scaling = %lg \\n\", G_w_y_scaling);
                      
                      Rprintf(\"quarantine_start_time = %lg \\n\", quarantine_start_time);
                      Rprintf(\"PCR_sens = %lg \\n\", PCR_sens);
                      Rprintf(\"sigma_M = %lg \\n\", sigma_M);
                      
                      Rprintf(\"beta_w_3 = %lg \\n\", beta_w_3);
                      Rprintf(\"beta_w_2 = %lg \\n\", beta_w_2);
                      Rprintf(\"beta_w_1 = %lg \\n\", beta_w_1);
                      Rprintf(\"beta_w_0 = %lg \\n\", beta_w_0);
                      Rprintf(\"g_0 = %lg \\n\", g_0);
                      Rprintf(\"g_F = %lg \\n\", g_F);
                      Rprintf(\"sigma_epsilon = %lg \\n\", sigma_epsilon);
                      
                      //Print out exposed compartments
                      for(m=0; m<M; m++) {
                        Rprintf(\"m = %d \\n\", m);
                        Rprintf(\"e[m] = %lg \\n\", e[m]);
                      }
                      
                      //Print out Q_1 and Q_NC compartments
                      for(v=0; v<V; v++) {
                        Rprintf(\"v = %d \\n\", v);
                        Rprintf(\"q_1[v] = %lg \\n\", q_1[v]);
                        Rprintf(\"q_nc[v] = %lg \\n\", q_nc[v]);
                        Rprintf(\"p_q1[v] = %lg \\n\", p_q1[v]);
                        Rprintf(\"y_q1[v] = %lg \\n\", y_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_q1[v] = %lg \\n\", total_samples_for_PCR_Testing_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_qnc[v] = %lg \\n\", total_samples_for_PCR_Testing_qnc[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_1_q1[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_1_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_1_qnc[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_1_qnc[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_2_q1[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_2_q1[v]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_2_qnc[v] = %lg \\n\", total_samples_for_PCR_Testing_lag_2_qnc[v]);
                      }
                      
                      //Print out Q_3 compartments
                      for(k=0; k<K; k++) {
                        Rprintf(\"k = %d \\n\", k);
                        Rprintf(\"q_3[k] = %lg \\n\", q_3[k]);
                        Rprintf(\"p_q3[k] = %lg \\n\", p_q3[k]);
                        Rprintf(\"y_q3[k] = %lg \\n\", y_q3[k]);
                        Rprintf(\"total_samples_for_PCR_Testing_q3[k] = %lg \\n\", total_samples_for_PCR_Testing_q3[k]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_1_q3[k] = %lg \\n\", total_samples_for_PCR_Testing_lag_1_q3[k]);
                        Rprintf(\"total_samples_for_PCR_Testing_lag_2_q3[k] = %lg \\n\", total_samples_for_PCR_Testing_lag_2_q3[k]);
                      }
                      

                      Error_Printing_Complete = TRUE;
                  }
                  
                   
                   ")
                  




# ---- init ----
init <- Csnippet("
                  //Rprintf(\"At top of init N_0 = %lg \\n\", N_0);
                  //Rprintf(\"At top of init E_0 = %lg \\n\", E_0);
                  
                  int M = (int) M_0; //Number of exposed compartments
                  int V = (int) V_0; //Number of days spent in hospital (number of cohorts in Queues 1 and NC)
                  int K = (int) K_0; //Number of days spent in quarantine (number of cohorts in Queue 3)
                  
                  int m; //Exposed compartment number
                  int v; //Queue 1/ Queue NC cohort number
                  int k; //Queue 3 cohort number
                  
                  //Declare E pointer array
                  double *e=&E_1;
                  
                  //Declare Queue 1 pointer arrays
                  double *q_1=&Q_1_1;
                  double *p_q1=&P_Q1_1;
                  double *total_samples_for_PCR_Testing_q1=&total_samples_for_PCR_Testing_Q1_1;
                  double *total_samples_for_PCR_Testing_lag_1_q1=&total_samples_for_PCR_Testing_lag_1_Q1_1;
                  double *total_samples_for_PCR_Testing_lag_2_q1=&total_samples_for_PCR_Testing_lag_2_Q1_1;
                  double *y_q1=&Y_Q1_1;
                  
                  //Declare Queue NC pointer arrays
                  double *q_nc=&Q_NC_1;
                  double *total_samples_for_PCR_Testing_qnc=&total_samples_for_PCR_Testing_QNC_1;
                  double *total_samples_for_PCR_Testing_lag_1_qnc=&total_samples_for_PCR_Testing_lag_1_QNC_1;
                  double *total_samples_for_PCR_Testing_lag_2_qnc=&total_samples_for_PCR_Testing_lag_2_QNC_1;
                  double *y_qnc=&Y_QNC_1;
                  
                  //Declare Queue 3 pointer arrays
                  double *q_3=&Q_3_1;
                  double *p_q3=&P_Q3_1;
                  double *total_samples_for_PCR_Testing_q3=&total_samples_for_PCR_Testing_Q3_1;
                  double *total_samples_for_PCR_Testing_lag_1_q3=&total_samples_for_PCR_Testing_lag_1_Q3_1;
                  double *total_samples_for_PCR_Testing_lag_2_q3=&total_samples_for_PCR_Testing_lag_2_Q3_1;
                  double *y_q3=&Y_Q3_1;



                  double E_init_total = 0;
                  double I_init_total = 0;
                  
                  if(z_0 > N_0){
                      I_init_total = nearbyint(N_0);
                      E_init_total = 0;
                      S = 0;
                  }else{
                      if(E_0 > N_0){
                          E_init_total = nearbyint(N_0);
                          I_init_total = 0;
                          S = 0;
                      }else{
                        E_init_total = nearbyint(E_0);
                        int extra_cap = nearbyint(N_0) - nearbyint(E_0);
                        if(extra_cap < nearbyint(z_0)){
                          I_init_total = nearbyint(extra_cap);
                          S = 0;
                        }else{
                          I_init_total = nearbyint(z_0);
                          S = nearbyint(N_0) - nearbyint(z_0) - nearbyint(E_0);
                        }
                      }
                  }
                  
                  //Assign early stage infections
                  double time_pre_symp = 1/phi_U;
                  double time_symp = 1/phi_S;
                  double prop_time_pre_symp = time_pre_symp/(time_pre_symp + time_symp);
                  
                  double total_init_I_symp = nearbyint(p_S*I_init_total);
                  I_A = nearbyint((1-p_S)*I_init_total);
                  
                  I_P = nearbyint(prop_time_pre_symp*total_init_I_symp);
                  I_S_1 = nearbyint((1-prop_time_pre_symp)*total_init_I_symp);
                  
                  //Late stage infection compartments
                  I_S_2 = 0;
                  I_H = 0;

                  //Recovered Compartments
                  R_A = 0;
                  R_F = 0;
                  R_H = 0;
  
                  //Whole Population
                  N = nearbyint(N_0);
                  
                  //Asymptomatic Testing Positive Cases (in isolation)
                  A_T = 0;
                  
                  //Transmission Rate
                  double total_time_infected = (1/gamma) + (1/phi_S);
                  double gamma_total = 1/total_time_infected;
                  double Beta_0 = R_0*(gamma_total);
                  beta_t = Beta_0;
                  //Reported Cases
                  C_Q1 = nearbyint(p_H_cond_S*C_0);
                  C_Q2 = 0;
                  C_Q3 = nearbyint((1-p_H_cond_S)*C_0);
                  C_Q4 = 0;
                  
                  
                  
                  Y_sum = 0;
                  
                  //Queue
                  L_int = nearbyint(L_advanced_2_days);
                  L_1 = 0;
                  L_2 = 0;
                  L_3 = 0;
                  L_4 = 0;
                  
                  Prop_Positive_Tests_Track = 0;
                  
                  G_w_y = 0;
                  
                  Q_2 = 0;
                  Q_4 = 0;
                  
                  
                  Backlog_Queue_1 = 0;
                  Backlog_Queue_NC = 0;
                  Backlog_Queue_3 = 0;
                  
                  First_re_test_Q1 = 0;
                  First_re_test_Q3 = 0;
                  
                  Second_re_test_Q1 = 0;
                  Second_re_test_Q3 = 0;
                  
                  total_samples_for_PCR_Testing_backlog_QNC = 0;
                  total_samples_for_PCR_Testing_backlog_Q1 = 0;
                  total_samples_for_PCR_Testing_backlog_Q3 = 0;
                  
                  total_samples_for_PCR_Testing_backlog_lag_1_Q1 = 0;
                  total_samples_for_PCR_Testing_backlog_lag_1_Q3 = 0;
                  total_samples_for_PCR_Testing_backlog_lag_1_QNC = 0;
                  
                  total_samples_for_PCR_Testing_backlog_lag_2_Q1 = 0;
                  total_samples_for_PCR_Testing_backlog_lag_2_Q3 = 0;
                  total_samples_for_PCR_Testing_backlog_lag_2_QNC = 0;
                  
                  total_samples_for_PCR_Testing_Q2 = 0;
                  total_samples_for_PCR_Testing_Q4 = 0;
                  
                  neg_samples_Q1 = 0;
                  neg_samples_Q3 = 0;
                  neg_samples_Q5 = 0;
                  total_neg_samples_all_queues = 0;
                  
                  Y_Q1_backlog = 0;
                  Y_QNC_backlog = 0;
                  Y_Q3_backlog = 0;
                  
                  //Initialize Q5 variables
                  infected_sample_size = 0;
                  total_sample_size = 0;
                  prob_infected_Q5 = 0;
                  total_samples_for_PCR_Testing_Q5 = 0;
                  Y_Q5 = 0;
                  
                  //Initialize Arrays
                  //Exposed compartment (including first one)
                  for(m = 0; m < M; m++){
                    e[m] = nearbyint(E_init_total/5);
                  }
                  
                  //Queue 1 and Queue NC arrays
                  for(v = 0; v < V; v++){
                    q_1[v] = 0;
                    q_nc[v] = 0;
                    p_q1[v] = 0;
                    total_samples_for_PCR_Testing_q1[v] = 0;
                    total_samples_for_PCR_Testing_qnc[v] = 0;
                    total_samples_for_PCR_Testing_lag_1_q1[v] = 0;
                    total_samples_for_PCR_Testing_lag_1_qnc[v] = 0;
                    total_samples_for_PCR_Testing_lag_2_q1[v] = 0;
                    total_samples_for_PCR_Testing_lag_2_qnc[v] = 0;
                    y_q1[v] = 0;
                    y_qnc[v] = 0;
                  }
                  
                  //Queue 3 arrays
                  for(k = 0; k < K; k++){
                    q_3[k] = 0;
                    p_q3[k] = 0;
                    total_samples_for_PCR_Testing_q3[k] = 0;
                    total_samples_for_PCR_Testing_lag_1_q3[k] = 0;
                    total_samples_for_PCR_Testing_lag_2_q3[k] = 0;
                    y_q3[k] = 0;
                  }
                  
                  //Set Flags
                  Neg_State_Value_Detected = FALSE;
                  NAN_State_Value_Detected = FALSE;
                  Error_Printing_Complete = FALSE;
                  
                  //Rprintf(\"At init N = %lg \\n\", N);
                  //Rprintf(\"At init S = %lg \\n\", S);
                  //Rprintf(\"At init E_1 = %lg \\n\", E_1);
                  //Rprintf(\"At init I_P = %lg \\n\", I_P);
                  //Rprintf(\"At init I_S_1 = %lg \\n\", I_S_1);
                  //Rprintf(\"At init C_Q1 = %lg \\n\", C_Q2);

                 ")

# ---- rmeas ----
rmeas <- Csnippet("
                  double size = 1.0/sigma_M/sigma_M;
                  Y = rnbinom_mu(size,Y_sum);
                  double prop_sd = sqrt(Prop_Positive_Tests_Track*(1-Prop_Positive_Tests_Track)/L_int);
                  obs_prop_positive = rnorm(Prop_Positive_Tests_Track, prop_sd);
                  ")

# ---- dmeas ----
dmeas <- Csnippet("
                  if(isnan(Y)){
                    lik = 0;
                  }else{
                    if(G_w_y_scaling > 0.33){
                      lik = -39;
                    }else{
                      double size = 1.0/sigma_M/sigma_M;
                      static double tol = 0.1;
                      double prop_sd = sqrt(Prop_Positive_Tests_Track*(1-Prop_Positive_Tests_Track)/L_int);
                      double lik_2 = dnorm(obs_prop_positive,Prop_Positive_Tests_Track,prop_sd,1);
                      lik = dnbinom_mu(Y,size,Y_sum+tol,1);
                    }
                    
                    
                  }
                      
                  
                  //Debugging Print Code
                  //Rprintf(\"t = %lg \\n\", t);
                  //Rprintf(\"I_S_1 = %lg \\n\", I_S_1);
                  //Rprintf(\"Lik = %lg \\n\", lik);
                  //Rprintf(\"Y = %lg \\n\", Y);
                  //Rprintf(\"Y_sum = %lg \\n\", Y_sum);
                  //Rprintf(\"tol = %lg \\n\", tol);
                  //Rprintf(\"size = %lg \\n\", size);

                  if (!give_log) lik = exp(lik);
                  ")

# ---- par_trans ----
par_trans = parameter_trans(log = c("R_0", "gamma", "h_V",
                                    "phi_E", "phi_U", "phi_S",
                                    "E_0", "z_0", "N_0", "C_0", "sigma_M"),
                            logit = c("PCR_sens", "p_S", "p_H_cond_S",
                                      "b_q", "b_a", "b_p", "G_w_y_scaling"))




