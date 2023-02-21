#Updating Disease history status
update_risk_factor <- function(sim_out_t, risk_factor, predictor, s) {
  # Updates risk factor values from time t to time t+1
  # References risk_factor_trend_sim for APC trends in Total_Chol, HDL, SBP/DBP, BMI, Trig, and Glucose
  risk_factor_vals <- case_when(
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHWM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHWM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHWM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHWM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHWF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHWF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHWF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHWF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHBM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHBM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHBM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHBM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHBF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHBF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHBF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHBF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "HM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "HM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "HM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "HM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "HF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "HF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "HF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "HF"),s+5]/100)), 
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "Male" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "Male"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "Male" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "Male"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "Female" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "Female"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "Female" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "Female"),s+5]/100)),
    TRUE ~ as.numeric(NA))
  
  return(risk_factor_vals)
}



run_sim <- function(s, intervention) {
  
  if (intervention == "Policy") {
    policy_effect_fruit <-   policy_effect_fruit_sim
    policy_effect_veg <-  policy_effect_veg_sim
  } else {
    policy_effect_fruit <- policy_effect_veg <- rep(0, n.sim)
  }
  
  ##$$Policy specific settings: produce RX project only
  
  if (scenario=="Hba1c&BMI" & intervention=="Policy") {
    HCE_parameter_sim[8,s]<-HCE_parameter_sim[8,s]*(1+policy_effect_HA1c_sim[s]*0.13)
  } else { 
    HCE_parameter_sim[8,s]= HCE_parameter_sim[8,s]
  }
  
  # Runs single simulation
  # Returns array with variables for each individual after n.cycle
  
  # Define inpute variables
  Input_variables <- c("SEQN","DEMO", "Age", "Total_Chol","HDL", "SBP", "DBP", "BMI",  "Obesity", "weight", "height","Trig", "Glucose", 
                       "Diabetes", "CVD_history","DM_prob", "CVD_prob", "CVD_recurrent_prob","HRQOL_scores", "HCE_predict")
  # Define output variables
  Out_variables<-c( "SEQN","Incident First CVD", "Incident Recurrent CVD","Incident CVD","Life Years","effect_disc", "HCE_disc", "Prod_cost", "Food_cost", "Food_cost_disc",  "Prod_cost_disc","Admin_costs", "Total_cost_health","Total_cost_societ" )
  
  #Set up initial values
  #replicate each individual in the data for n.loop times, to minimize stochastic error
  
  sim_out_t<-do.call("rbind", replicate(n.loop, data_for_analysis[,Input_variables], simplify = FALSE))
  sim_out_t$state=rep(data_for_analysis$initial_H, n.loop)
  sim_out_t[,"Age_cycle"]=sim_out_t[,"Age"]
  sim_out_t[,"BaseBMI"]=sim_out_t[,"BMI"]
  sim_out_t[,"HCE_disc"] <- sim_out_t$HCE_predict
  sim_out_t[,"effect_disc"] <- sim_out_t$HRQOL_scores
  #Assume productivity costs to be 0 at baseline for both intervention and control group, since no difference has occur 
  sim_out_t[,"Prod_cost"] <-0
  sim_out_t[,"Prod_cost_disc"] <-0
  #Assume productivity costs to be 0 at baseline for both intervention and control group, since no difference has occur 
  sim_out_t[,"Admin_costs"] <-0
  sim_out_t[,"Life Years"]=0
  sim_out_t[,"Incident First CVD"]=0
  sim_out_t[,"Incident Recurrent CVD"]=0
  sim_out_t[,"Incident CVD"]=0
  sim_out_t[,"Food_cost"] <- 0
  sim_out_t[,"Food_cost_disc"] <- 0
  
  for (t in 1:n.cycle) {
    print(t)
  
    #Time-varying data inputs
    sim_out_t[,"Age_cycle"] <- sim_out_t[,"Age"] + t
    
    # Updating Risk Factors over time: Applying secular trends based on age, gender, R/E
    
    ##$$Policy specific settings: produce RX project only
    # Adjust BMI based on change in fruit and vegetable intake (only applies in first 4 years)
    # Add sensivity scenario: change the BMI pathway. Change only apply once
    if (intervention == "Policy" & t <= 4 & scenario == "main") {
    weight = as.numeric(sim_out_t[,"weight"]) + t*wt_change_fruit[s]*policy_effect_fruit[s]+ t*wt_change_veg[s]*policy_effect_veg[s]
     sim_out_t[,"BMI"] = round(weight / (as.numeric(sim_out_t[,"height"])/100)^2, 1)
     sim_out_t[,"BMI"] <- update_risk_factor(sim_out_t, "BMI", "BMI", s)
    } else if (intervention == "Policy" & t <= 1  & (scenario=="BMI"|scenario=="Hba1c&BMI")){ 
     sim_out_t[,"BMI"]= as.numeric(sim_out_t[,"BMI"])+ policy_effect_BMI_sim[s]
     sim_out_t[,"BMI"] <- update_risk_factor(sim_out_t, "BMI", "BMI", s)   
    }  else {
      sim_out_t[,"BMI"] <- update_risk_factor(sim_out_t, "BMI", "BMI", s)
    }
    ##$$
    #without effect on indirect pathway
    #sim_out_t[,"BMI"] <- update_risk_factor(sim_out_t, "BMI", "BMI", s)
    # Truncate BMI to the interval 12-70 kg/m2. Source: ???	https://doi.org/10.2105/AJPH.2008.137364  
   sim_out_t[,"BMI"] <- ifelse(sim_out_t[,"BMI"] < 12, 12,
                                ifelse(sim_out_t[,"BMI"]>70, 70, sim_out_t[,"BMI"]))
   sim_out_t[,"Total_Chol"] <- update_risk_factor(sim_out_t, "Total_Chol", "Total_Chol", s)
   sim_out_t[,"HDL"] <- update_risk_factor(sim_out_t, "HDL", "HDL", s)
   sim_out_t[,"SBP"] <- update_risk_factor(sim_out_t, "SBP", "SBP", s)
   sim_out_t[,"DBP"] <- update_risk_factor(sim_out_t, "DBP", "SBP", s)
   sim_out_t[,"Trig"] <- update_risk_factor(sim_out_t, "Trig", "Trig", s)
   sim_out_t[,"Glucose"] <- update_risk_factor(sim_out_t, "Glucose", "Glucose", s)
    
    # Fruit-IHD/CHD
    RR_diff_fruit_hstk <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_fruit_hstk[1,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_fruit_hstk[2,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_fruit_hstk[3,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_fruit_hstk[4,s]*policy_effect_fruit[s]),
      TRUE ~ exp(random_logrr_fruit_hstk[5,s]*policy_effect_fruit[s])
    )
    RR_diff_fruit_ihd <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_fruit_ihd[1,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_fruit_ihd[2,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_fruit_ihd[3,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_fruit_ihd[4,s]*policy_effect_fruit[s]),
      TRUE ~ exp(random_logrr_fruit_ihd[5,s]*policy_effect_fruit[s])
    )
    RR_diff_fruit_istk <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_fruit_istk[1,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_fruit_istk[2,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_fruit_istk[3,s]*policy_effect_fruit[s]),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_fruit_istk[4,s]*policy_effect_fruit[s]),
      TRUE ~ exp(random_logrr_fruit_istk[5,s]*policy_effect_fruit[s])
    )
    
    # Veg-IHD/CHD
    RR_diff_veg_hstk <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_veg_hstk[1,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_veg_hstk[2,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_veg_hstk[3,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_veg_hstk[4,s]*policy_effect_veg[s]),
      TRUE ~ exp(random_logrr_veg_hstk[5,s]*policy_effect_veg[s])
    )
    RR_diff_veg_ihd <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_veg_ihd[1,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_veg_ihd[2,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_veg_ihd[3,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_veg_ihd[4,s]*policy_effect_veg[s]),
      TRUE ~ exp(random_logrr_veg_ihd[5,s]*policy_effect_veg[s])
    )
    RR_diff_veg_istk <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_veg_istk[1,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_veg_istk[2,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_veg_istk[3,s]*policy_effect_veg[s]),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_veg_istk[4,s]*policy_effect_veg[s]),
      TRUE ~ exp(random_logrr_veg_istk[5,s]*policy_effect_veg[s])
    )
    
    
    # BMI-IHD/CHD
    RR_diff_bmi_ihd <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_BMI_IHD[1,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5), 
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_BMI_IHD[2,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_BMI_IHD[3,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_BMI_IHD[4,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5),  
      TRUE ~ exp(random_logrr_BMI_IHD[5,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5)
    )
    
    # BMI-Total stroke
    RR_diff_bmi_tstk <- case_when(
     sim_out_t[,"Age_cycle"] < 45 ~ exp(random_logrr_BMI_TSTK[1,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5), 
     sim_out_t[,"Age_cycle"] < 55 ~ exp(random_logrr_BMI_TSTK[2,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5),
     sim_out_t[,"Age_cycle"] < 65 ~ exp(random_logrr_BMI_TSTK[3,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5),
     sim_out_t[,"Age_cycle"] < 75 ~ exp(random_logrr_BMI_TSTK[4,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5),
      TRUE ~ exp(random_logrr_BMI_TSTK[5,s]*(sim_out_t[,"BMI"] - sim_out_t[,"BaseBMI"])/5)
    )
    # Calculate weighted average of HSTK and ISTK risk ratios to get risk ratio for total stroke
    # 87% of strokes are ischemic (Source: Benjamin et al. Heart Disease and Stroke StatisticsÃ¢â,¬â???2019 Update: A Report From the American Heart Association)
    RR_diff_fruit_tstk = 0.87*RR_diff_fruit_istk + 0.13*RR_diff_fruit_hstk
    RR_diff_veg_tstk = 0.87*RR_diff_veg_istk + 0.13*RR_diff_veg_hstk
    
    # Combine effects by multiplying RRs
    RR_diff_total_CHD <- RR_diff_bmi_ihd*RR_diff_fruit_ihd*RR_diff_veg_ihd
    RR_diff_total_Stroke <- RR_diff_bmi_tstk*RR_diff_fruit_tstk*RR_diff_veg_tstk
    

    #Updating the disease risk
    
    vars<-c("SEQN", "Age_cycle","Glucose","BMI","Total_Chol","HDL","Trig","SBP","DBP","Diabetes","CVD_history","Age")
    raw.input.data<-cbind(do.call("rbind", replicate(n.loop, data_for_analysis[,c("Female", "Race", "DM_parent",  "HPT_Txt",  "Smoking")], simplify = FALSE)),sim_out_t[,vars])
    raw.input.data$Age<-raw.input.data$Age_cycle
    
    if (t%/%2 > 0 & t%%2 == 0)  {
      CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
     sim_out_t[,"CVD_recurrent_prob"] <- Multi_yr_Risk_to_annual_prob(time=2, risk=CVD_Recurrent_risk_2yr)
    } else {
     sim_out_t[,"CVD_recurrent_prob"] <-sim_out_t[,"CVD_recurrent_prob"]
    }
    
    if (t%/%8 > 0 & t%%8 == 0)  {
      DM_risk_8yr <- calc_DM_risk(raw.input.data)
     sim_out_t[,"DM_prob"] <- Multi_yr_Risk_to_annual_prob(time=8, risk=DM_risk_8yr)
    } else {
     sim_out_t[,"DM_prob"] <-sim_out_t[,"DM_prob"]
    }
    
    if (t%/%10 > 0 & t%%10 == 0) {
      ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
     sim_out_t[,"CVD_prob"] <- Multi_yr_Risk_to_annual_prob(time=10, risk=ASCVD_Risk_10yr)
    } else {
     sim_out_t[,"CVD_prob"] <-sim_out_t[,"CVD_prob"]
    }
    
    #Disaggregating ASCVD
  p.CHD_first <- rep(as.numeric(Prop_CHD[1,data_for_analysis[,"DEMO"]]), n.loop) #Prob. of CHD event (Angina, MI, Fatal MI, Fatal CHD) Source: Benjamin et al, Circulation, 2018 Table 13-1 & 18-1, and 18-2
  p.CHD_recurrent <-  rep(as.numeric(Prop_CHD[2,data_for_analysis[,"DEMO"]]), n.loop)
  
    #Defining transition probabilities
    # Subtract 19 from age because indexing starts at age 20
    p.death.DM <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ DM_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ DM_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ DM_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ DM_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ DM_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ DM_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ DM_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ DM_mortality_HM[sim_out_t[,"Age_cycle"]-19,s]
    )
    
    p.death.CHD <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ CHD_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ CHD_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ CHD_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ CHD_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ CHD_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ CHD_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ CHD_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ CHD_mortality_HM[sim_out_t[,"Age_cycle"]-19,s]
    )
    
    p.death.Stroke <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ stroke_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ stroke_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ stroke_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ stroke_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ stroke_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ stroke_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ stroke_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ stroke_mortality_HF[sim_out_t[,"Age_cycle"]-19,s]
    )
    
    p.death <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ Non_CVD_DM_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ Non_CVD_DM_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ Non_CVD_DM_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ Non_CVD_DM_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ Non_CVD_DM_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ Non_CVD_DM_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ Non_CVD_DM_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ Non_CVD_DM_mortality_HF[sim_out_t[,"Age_cycle"]-19,s]
    )
    

    #Markov State #1: "No CVD, No Diabetes"
    p.H.2.DM <- as.numeric(sim_out_t[,"DM_prob"])*rep(data_for_analysis$risk_adjustment.DM, n.loop)
    p.H.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(sim_out_t[,"CVD_prob"])
    p.H.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(sim_out_t[,"CVD_prob"])
    p.H.2.H <- ifelse((p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death) > 1, 0,
                      1 - (p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death))
    
    #Markov State #2: "No CVD, With Diabetes
    p.DM.2.death <- (1-exp(-(p.death+p.death.DM)))
    p.DM.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(sim_out_t[,"CVD_prob"])
    p.DM.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(sim_out_t[,"CVD_prob"])
    p.DM.2.DM <- as.numeric(ifelse((p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death) > 1, 0,
                                   1 - (p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death)))
    
    #Markov State #3: "First Stroke"
    p.initial_Stroke.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke)))))
    p.initial_Stroke.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.initial_Stroke.2.death))
    p.initial_Stroke.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.initial_Stroke.2.death, 0))
    
    #Markov State #4: "First CHD w/o RVSC" 
    p.initial_CHD_No_RVSC.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD)))))
    p.initial_CHD_No_RVSC.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.initial_CHD_No_RVSC.2.death))
    p.initial_CHD_No_RVSC.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.initial_CHD_No_RVSC.2.death, 0))
    
    #Markov State #5: "First CHD with RVSC"
    p.initial_CHD_RVSC.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s])))))
    p.initial_CHD_RVSC.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.initial_CHD_RVSC.2.death))
    p.initial_CHD_RVSC.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.initial_CHD_RVSC.2.death, 0))
    
    #Markov State #6: "CVD History, No Diabetes"
    p.CVD.2.DM <- as.numeric(sim_out_t[,"DM_prob"])*rep(data_for_analysis$risk_adjustment.DM, n.loop)
    p.CVD.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke)))
    p.CVD.2.CVD <- as.numeric(ifelse((p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death) > 1, 0,
                                     1 - (p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death)))
    
    #Markov State #7: "CVD History, With Diabetes"
    p.CVD_DM.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD_DM.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD_DM.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke+p.death.DM)))
    p.CVD_DM.2.CVD_DM <- as.numeric(ifelse((p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke) > 1, 0,
                                           1 - (p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke)))
    
    #Markov State #8: "Subsequent Stroke"
    p.sub_Stroke.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke)))))
    p.sub_Stroke.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.sub_Stroke.2.death))
    p.sub_Stroke.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.sub_Stroke.2.death, 0))
    
    #Markov State #9: "Subsequent CHD w/o RVSC" 
    p.sub_CHD_No_RVSC.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD)))))
    p.sub_CHD_No_RVSC.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.sub_CHD_No_RVSC.2.death))
    p.sub_CHD_No_RVSC.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.sub_CHD_No_RVSC.2.death, 0))
    
    #Markov State #10: "Subsequent CHD with RVSC"
    p.sub_CHD_RVSC.2.death <-  as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s])))))
    p.sub_CHD_RVSC.2.CVD_No_DM <-  as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.sub_CHD_RVSC.2.death))
    p.sub_CHD_RVSC.2.CVD_DM <-  as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.sub_CHD_RVSC.2.death, 0))
    
    #Assign transition probablities
    p.transition <- array(NA, dim=c(n.individual, n.health.state),
                          dimnames =  list(c(1:n.individual),name.health.state))
    p.transition[,"No CVD, No Diabetes"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.H.2.H, rep(0, n.individual))
    p.transition[,"No CVD, With Diabetes"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.H.2.DM, 
                                                     ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", p.DM.2.DM, rep(0, n.individual)))
    p.transition[,"First Stroke"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.H.2.initial_Stroke, 
                                            ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", p.DM.2.initial_Stroke, rep(0, n.individual)))
    p.transition[,"First CHD w/o RVSC"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", (1-p.RVSC)*p.H.2.initial_CHD, 
                                                  ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", (1-p.RVSC)*p.DM.2.initial_CHD, rep(0, n.individual)))
    p.transition[,"First CHD with RVSC"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.RVSC*p.H.2.initial_CHD, 
                                                   ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", p.RVSC*p.DM.2.initial_CHD, rep(0, n.individual)))
    p.transition[,"CVD History, No Diabetes"] <- case_when(
      sim_out_t[,"state"] == "First Stroke" ~ p.initial_Stroke.2.CVD_No_DM, 
      sim_out_t[,"state"] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_No_DM, 
      sim_out_t[,"state"] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_No_DM, 
      sim_out_t[,"state"] == "CVD History, No Diabetes" ~ p.CVD.2.CVD, 
      sim_out_t[,"state"] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_No_DM, 
      sim_out_t[,"state"] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_No_DM, 
      sim_out_t[,"state"] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_No_DM, 
      TRUE ~ rep(0, n.individual)
    )
    p.transition[,"CVD History, With Diabetes"] <- case_when(
      sim_out_t[,"state"] == "First Stroke" ~ p.initial_Stroke.2.CVD_DM, 
      sim_out_t[,"state"] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_DM, 
      sim_out_t[,"state"] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_DM, 
      sim_out_t[,"state"] == "CVD History, No Diabetes" ~ p.CVD.2.DM, 
      sim_out_t[,"state"] == "CVD History, With Diabetes" ~ p.CVD_DM.2.CVD_DM, 
      sim_out_t[,"state"] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_DM, 
      sim_out_t[,"state"] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_DM, 
      sim_out_t[,"state"] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_DM, 
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"Subsequent Stroke"] <- ifelse(sim_out_t[,"state"] == "CVD History, No Diabetes", p.CVD.2.Sub_Stroke,
                                                 ifelse(sim_out_t[,"state"] == "CVD History, With Diabetes", p.CVD_DM.2.Sub_Stroke, rep(0, n.individual)))
    p.transition[,"Subsequent CHD w/o RVSC"] <- ifelse(sim_out_t[,"state"] == "CVD History, No Diabetes", (1-p.RVSC)*p.CVD.2.Sub_CHD, 
                                                       ifelse(sim_out_t[,"state"] == "CVD History, With Diabetes", (1-p.RVSC)*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
    p.transition[,"Subsequent CHD with RVSC"] <- ifelse(sim_out_t[,"state"] == "CVD History, No Diabetes", p.RVSC*p.CVD.2.Sub_CHD, 
                                                        ifelse(sim_out_t[,"state"] == "CVD History, With Diabetes", p.RVSC*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
    p.transition[,"Death"] <- case_when(
      sim_out_t[,"state"] == "No CVD, No Diabetes" ~ p.death, 
      sim_out_t[,"state"] == "No CVD, With Diabetes" ~ p.DM.2.death, 
      sim_out_t[,"state"] == "First Stroke" ~ p.initial_Stroke.2.death, 
      sim_out_t[,"state"] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.death, 
      sim_out_t[,"state"] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.death, 
      sim_out_t[,"state"] == "CVD History, No Diabetes" ~ p.CVD.2.death, 
      sim_out_t[,"state"] == "CVD History, With Diabetes" ~ p.CVD_DM.2.death, 
      sim_out_t[,"state"] == "Subsequent Stroke" ~ p.sub_Stroke.2.death, 
      sim_out_t[,"state"] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.death, 
      sim_out_t[,"state"] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.death,
      TRUE ~rep(1, n.individual) )
    
    # Check that all transition probabilities add to 1 (rounded to 3 digits)
    # if (sum(round(rowSums(p.transition),3)!=1) != 0) {
    #   p_sums = round(rowSums(p.transition),3)
    #   error_out = sim_out[p_sums!=1,"state",t] # Output state of person(s) with error
    #   stop("Transition probabilities do not add to 1. ", paste("Simulation", s, ", Time", t, ". "), "State(s) with error: ", error_out)
    # }
    # Transition to the next health state 
    
   sim_out_t[,"state"]<- apply(p.transition, 1, function(x) sample(name.health.state, 1, prob = x))
    
   sim_out_t[,"Diabetes"] <- ifelse(sim_out_t[,"state"]%in% c("No CVD, With Diabetes", "CVD History, With Diabetes"), 1, 
                                     ifelse(sim_out_t[,"state"]== "Death", sim_out_t[,"Diabetes"], 
                                            ifelse(rep(data_for_analysis$Glucose, n.loop) > 126, 1,sim_out_t[,"Diabetes"])))
    
   sim_out_t[,"CVD_history"] <- ifelse(sim_out_t[,"CVD_history"]==1, 1,
                                        ifelse( sim_out_t[,"state"] %in% c("First Stroke", "First CHD w/o RVSC", "First CHD with RVSC", "CVD History, No Diabetes", "CVD History, With Diabetes",
                                                                   "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), 1, 0))
    
    
    
   sim_out_t[,"BMI"] <- ifelse(sim_out_t[,"state"]== "Death",NA,sim_out_t[,"BMI"])
   sim_out_t[,"Obesity"] <- ifelse(sim_out_t[,"BMI"] < 30, 0, 
                                    ifelse(sim_out_t[,"BMI"] >= 30, 1, sim_out_t[,"Obesity"]))
    
    
    raw.input.data[,vars]<-sim_out_t[,vars]
       
    # Update QALYs and costs
    HRQOL_scores_t1 <- calc_HRQOL(raw.input.data, HRQOL_parameter_sim[,s])
    HRQOL_scores_t1 <- ifelse(sim_out_t[,"state"]%in% c("First Stroke", "Subsequent Stroke"), HRQOL_scores_t1+u_stroke_sim[s], 
                                           ifelse(sim_out_t[,"state"]%in% c("First CHD w/o RVSC", "First CHD with RVSC", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), HRQOL_scores_t1+u_CHD_sim[s],
                                                  ifelse(sim_out_t[,"state"]== "Death", 0, HRQOL_scores_t1)))
    
    sim_out_t[,"HRQOL_scores"] =sim_out_t[,"HRQOL_scores"]+HRQOL_scores_t1
    sim_out_t[,"effect_disc"] <- sim_out_t[,"effect_disc"]+HRQOL_scores_t1/((1+beta_QALY)^(t-1))
    
    
    #ppp: for sensitivity analysis, add economic effect of HA1c
    
    
    HCE_predict_t1 <- calc_HCE(raw.input.data, HCE_parameter_sim[,s])
    HCE_predict_t1 <- ifelse(sim_out_t[,"state"]%in% c("First Stroke", "Subsequent Stroke"), HCE_predict_t1 +c_stroke_sim[s] , 
                                          ifelse(sim_out_t[,"state"]%in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC"), HCE_predict_t1+c_CHD_sim[s], 
                                                 ifelse(sim_out_t[,"state"]%in% c("First CHD with RVSC", "Subsequent CHD with RVSC"), HCE_predict_t1+c_CHD_sim[s]+c_RVSC_sim[s],
                                                        ifelse(sim_out_t[,"state"]== "Death", 0, HCE_predict_t1))))
    sim_out_t[,"HCE_predict"]=sim_out_t[,"HCE_predict"]+HCE_predict_t1
    sim_out_t[,"HCE_disc"] <- sim_out_t[,"HCE_disc"]+HCE_predict_t1/((1+beta_cost)^(t-1))
    
    
    
    #ppp Add productivity costs associated with stroke, and CHD respectively 
    Prod_cost_t1<- ifelse(sim_out_t[,"state"]%in% c("First Stroke", "Subsequent Stroke"),c_prod_stroke_sim[s], 
                                       ifelse(sim_out_t[,"state"]%in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC","First CHD with RVSC", "Subsequent CHD with RVSC"), c_prod_CHD_sim[s],
                                              ifelse(sim_out_t[,"state"]%in% c("CVD History, No Diabetes", "CVD History, With Diabetes"), c_prod_CVD_sim[s],      
                                                     ifelse(sim_out_t[,"state"]== "Death", 0, 0))))
    
    Prod_cost_t1<- ifelse(sim_out_t[,"state"]%in% c("First Stroke", "Subsequent Stroke"),c_prod_stroke_sim[s], 
                                       ifelse(sim_out_t[,"state"]%in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC","First CHD with RVSC", "Subsequent CHD with RVSC"), c_prod_CHD_sim[s],
                                              ifelse(sim_out_t[,"state"]%in% c("CVD History, No Diabetes", "CVD History, With Diabetes"), c_prod_CVD_sim[s],      
                                                     ifelse(sim_out_t[,"state"]== "Death", 0, 0))))
    
    sim_out_t[,"Prod_cost"] =as.numeric(sim_out_t[,"Prod_cost"])+Prod_cost_t1
    #ppp
    sim_out_t[,"Prod_cost_disc"] <- as.numeric(sim_out_t[,"Prod_cost_disc"])+Prod_cost_t1/((1+beta_cost)^(t-1))

  
    if (intervention == "Policy") {
      Food_cost_t1 <- ifelse(sim_out_t[,"state"]!= "Death", c_policy[s], 0)
    } else {
      Food_cost_t1 <- 0
    }
    
    sim_out_t[,"Food_cost"]<-sim_out_t[,"Food_cost"]+Food_cost_t1
    sim_out_t[,"Food_cost_disc"] <- sim_out_t[,"Food_cost_disc"]+Food_cost_t1/((1+beta_cost)^(t-1))
    
    # Discount effects and costs
    #ppp: add adiministrative costs to the model
   
    if (t == 1) {
      Admin_costs_t1 <- ifelse(sim_out_t[,"state"]!= "Death",  Food_cost_t1*0.352, 0)
    } else {
      Admin_costs_t1 <- ifelse(sim_out_t[,"state"]!= "Death",  Food_cost_t1*0.176, 0)
    }

    sim_out_t[,"Admin_costs"] <-sim_out_t[,"Admin_costs"]+ (Admin_costs_t1/((1+beta_cost)^(t-1)))
    
    sim_out_t[, "Total_cost_health"]<-sim_out_t[,"Food_cost_disc"]+ sim_out_t[,"Admin_costs"]+ sim_out_t[,"HCE_disc"]
    sim_out_t[, "Total_cost_societ"]<-sim_out_t[,"Food_cost_disc"]+ sim_out_t[,"Admin_costs"]+ sim_out_t[,"HCE_disc"]+
                                        sim_out_t[,"Prod_cost_disc"]

    
    sim_out_t[,"Life Years"]<-sim_out_t[,"Life Years"]+(sim_out_t[,"state"]!="Death")
    
    sim_out_t[,"Incident First CVD"] = sim_out_t[,"Incident First CVD"]+ ifelse(sim_out_t[,"state"]=="First Stroke"| sim_out_t[,"state"]=="First CHD with RVSC" |sim_out_t[,"state"]=="First CHD w/o RVSC",1, 0)
    
    sim_out_t[,"Incident Recurrent CVD"] = sim_out_t[,"Incident Recurrent CVD"]+ ifelse(sim_out_t[,"state"]=="Subsequent Stroke"|sim_out_t[,"state"]=="Subsequent CHD with RVSC"|sim_out_t[,"state"]=="Subsequent CHD w/o RVSC",1,0)
    
    sim_out_t[,"Incident CVD"] = sim_out_t[,"Incident Recurrent CVD"]+sim_out_t[,"Incident First CVD"]
  }
  sim_out_mean<-sim_out_t[ ,c(Out_variables)] %>% 
    group_by(SEQN)%>% 
    summarize_at(Out_variables[-1], mean, na.rm=TRUE)
  return(sim_out_mean)
}

