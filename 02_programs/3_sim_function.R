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
  if (scenario=="sensitivity2" & intervention=="Policy") {
    HCE_parameter_sim[8,s]<-HCE_parameter_sim[8,s]*(1+policy_effect_HA1c_sim[s]*0.13)
  } else { 
    
    HCE_parameter_sim[8,s]= HCE_parameter_sim[8,s]
  }
  
  # Runs single simulation
  # Returns array with variables for each individual at each cycle
  
  Input_variables <- c("SEQN", "Age",  "DEMO",
                       "Total_Chol","HDL", "SBP", "DBP", "BMI",  "Obesity", "weight", "height","Trig", "Glucose", 
                       "Diabetes", "CVD_history","DM_prob", "CVD_prob", "CVD_recurrent_prob")
  
  Output_variables<- c("SEQN", "HRQOL_scores", "HCE_predict","state", "HCE_disc", "effect_disc", "Food_cost", "Food_cost_disc", "Prod_cost", "Prod_cost_disc","Admin_costs")
  #ppp
  
  #aaa Initialize array to hold data inputs at t and t+1
  
  datain_t<-array(NA, dim=c(n.individual, length(Input_variables)+1),
                  dimnames = list(data_for_analysis$Subject_ID,
                                  c("Age_cycle", Input_variables)))
  datain_t1<-array(NA, dim=c(n.individual, length(Input_variables)+1),
                  dimnames = list(data_for_analysis$Subject_ID,
                                  c("Age_cycle", Input_variables)))
  

  ##aaa
  
  # Initialize array to hold output from single simulation
  sim_out<-array(NA, dim=c(n.individual, length(Output_variables), n.cycle+1),
                 dimnames = list(data_for_analysis$Subject_ID,
                                 Output_variables, paste("cycle", 0:n.cycle, sep = " ")))
  
  # Set initial state (cycle 0)
  
  ##aa
  datain_t[,Input_variables] <- as.matrix(data_for_analysis[,Input_variables])
  datain_t[,"Age_cycle"]=datain_t[,"Age"]
  sim_out[,c("SEQN", "HRQOL_scores", "HCE_predict"),1]<-as.matrix(data_for_analysis[,c("SEQN", "HRQOL_scores", "HCE_predict")])
  sim_out[,"state",1] <- initial_H
  sim_out[,"HCE_disc",1] <- data_for_analysis$HCE_predict
  sim_out[,"effect_disc",1] <- data_for_analysis$HRQOL_scores

  
  if (intervention == "Policy") {
    sim_out[,"Food_cost",1] <- c_policy
    sim_out[,"Food_cost_disc",1] <- c_policy
  } else {
    sim_out[,"Food_cost",1] <- 0
    sim_out[,"Food_cost_disc",1] <- 0
  }
  
  #ppp
  #Assume productivity costs to be 0 at baseline for both intervention and control group, since no difference has occur 
  sim_out[,"Prod_cost",1] <-0
  sim_out[,"Prod_cost_disc",1] <-0
  
  #Assume productivity costs to be 0 at baseline for both intervention and control group, since no difference has occur 
  sim_out[,"Admin_costs",1] <-0
  
  
  
  for (t in 1:n.cycle) {
    
    #Non-time varying data inputs: carry it over from the baseline data  
    datain_t1[,c("SEQN", "Age", "DEMO","weight","height")]<-datain_t[,c("SEQN", "Age", "DEMO","weight","height")]
    
    sim_out[,"SEQN", t]<-datain_t[,"SEQN"]
    #Time-varying data inputs
    datain_t1[,"Age_cycle"] <- as.numeric(data_for_analysis[,"Age"]) + t

# 5.4 Updating Risk Factors over time: Applying secular trends based on age, gender, R/E
# Adjust BMI based on change in fruit and vegetable intake (only applies in first 4 years)
# Add sensivity scenario: change the BMI pathway. Change only apply once


if (intervention == "Policy" & t <= 4 & scenario == "main") {
  weight = as.numeric(datain_t[,"weight"]) + t*wt_change_fruit[s]*policy_effect_fruit[s]+ t*wt_change_veg[s]*policy_effect_veg[s]
  datain_t1[,"BMI"] = round(weight / (as.numeric(datain_t[,"height"])/100)^2, 1)
  datain_t1[,"BMI"] <- update_risk_factor(datain_t1, "BMI", "BMI", s)
} else if (intervention == "Policy" & t <= 1  & (scenario=="sensitivity1"|scenario=="sensitivity2")){ 
  datain_t1[,"BMI"]= as.numeric(datain_t[,"BMI"])+ policy_effect_BMI_sim[s]
                                 datain_t1[,"BMI"] <- update_risk_factor(datain_t1, "BMI", "BMI", s)   
}  else {
  datain_t1[,"BMI"] <- update_risk_factor(datain_t, "BMI", "BMI", s)
}

# Truncate BMI to the interval 12-70 kg/m2. Source: ???	https://doi.org/10.2105/AJPH.2008.137364  
datain_t1[,"BMI"] <- ifelse(as.numeric(datain_t1[,"BMI"]) < 12, 12,
                             ifelse(as.numeric(datain_t1[,"BMI"])>70, 70, as.numeric(datain_t1[,"BMI"])))


datain_t1[,"Total_Chol"] <- update_risk_factor(datain_t, "Total_Chol", "Total_Chol", s)
datain_t1[,"HDL"] <- update_risk_factor(datain_t, "HDL", "HDL", s)
datain_t1[,"SBP"] <- update_risk_factor(datain_t, "SBP", "SBP", s)
datain_t1[,"DBP"] <- update_risk_factor(datain_t, "DBP", "SBP", s)
datain_t1[,"Trig"] <- update_risk_factor(datain_t, "Trig", "Trig", s)
datain_t1[,"Glucose"] <- update_risk_factor(datain_t, "Glucose", "Glucose", s)

# Fruit-IHD/CHD
RR_diff_fruit_hstk <- case_when(
  datain_t1[,"Age_cycle"] < 45 ~ exp(random_logrr_fruit_hstk[1,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 55 ~ exp(random_logrr_fruit_hstk[2,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 65 ~ exp(random_logrr_fruit_hstk[3,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 75 ~ exp(random_logrr_fruit_hstk[4,s]*policy_effect_fruit[s]),
  TRUE ~ exp(random_logrr_fruit_hstk[5,s]*policy_effect_fruit[s])
)
RR_diff_fruit_ihd <- case_when(
  datain_t1[,"Age_cycle"] < 45 ~ exp(random_logrr_fruit_ihd[1,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 55 ~ exp(random_logrr_fruit_ihd[2,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 65 ~ exp(random_logrr_fruit_ihd[3,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 75 ~ exp(random_logrr_fruit_ihd[4,s]*policy_effect_fruit[s]),
  TRUE ~ exp(random_logrr_fruit_ihd[5,s]*policy_effect_fruit[s])
)
RR_diff_fruit_istk <- case_when(
  datain_t1[,"Age_cycle"] < 45 ~ exp(random_logrr_fruit_istk[1,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 55 ~ exp(random_logrr_fruit_istk[2,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 65 ~ exp(random_logrr_fruit_istk[3,s]*policy_effect_fruit[s]),
  datain_t1[,"Age_cycle"] < 75 ~ exp(random_logrr_fruit_istk[4,s]*policy_effect_fruit[s]),
  TRUE ~ exp(random_logrr_fruit_istk[5,s]*policy_effect_fruit[s])
)

# Veg-IHD/CHD
RR_diff_veg_hstk <- case_when(
  datain_t1[,"Age_cycle"] < 45 ~ exp(random_logrr_veg_hstk[1,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 55 ~ exp(random_logrr_veg_hstk[2,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 65 ~ exp(random_logrr_veg_hstk[3,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 75 ~ exp(random_logrr_veg_hstk[4,s]*policy_effect_veg[s]),
  TRUE ~ exp(random_logrr_veg_hstk[5,s]*policy_effect_veg[s])
)
RR_diff_veg_ihd <- case_when(
  datain_t1[,"Age_cycle"] < 45 ~ exp(random_logrr_veg_ihd[1,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 55 ~ exp(random_logrr_veg_ihd[2,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 65 ~ exp(random_logrr_veg_ihd[3,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 75 ~ exp(random_logrr_veg_ihd[4,s]*policy_effect_veg[s]),
  TRUE ~ exp(random_logrr_veg_ihd[5,s]*policy_effect_veg[s])
)
RR_diff_veg_istk <- case_when(
  datain_t1[,"Age_cycle"] < 45 ~ exp(random_logrr_veg_istk[1,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 55 ~ exp(random_logrr_veg_istk[2,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 65 ~ exp(random_logrr_veg_istk[3,s]*policy_effect_veg[s]),
  datain_t1[,"Age_cycle"] < 75 ~ exp(random_logrr_veg_istk[4,s]*policy_effect_veg[s]),
  TRUE ~ exp(random_logrr_veg_istk[5,s]*policy_effect_veg[s])
)


# BMI-IHD/CHD
RR_diff_bmi_ihd <- case_when(
  datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_BMI_IHD[1,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5), 
  datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_BMI_IHD[2,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
  datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_BMI_IHD[3,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
  datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_BMI_IHD[4,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),  
  TRUE ~ exp(random_logrr_BMI_IHD[5,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5)
)

# BMI-Total stroke
RR_diff_bmi_tstk <- case_when(
  datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_BMI_TSTK[1,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5), 
  datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_BMI_TSTK[2,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
  datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_BMI_TSTK[3,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
  datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_BMI_TSTK[4,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
  TRUE ~ exp(random_logrr_BMI_TSTK[5,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5)
)
# Calculate weighted average of HSTK and ISTK risk ratios to get risk ratio for total stroke
# 87% of strokes are ischemic (Source: Benjamin et al. Heart Disease and Stroke Statisticsâ€”2019 Update: A Report From the American Heart Association)
RR_diff_fruit_tstk = 0.87*RR_diff_fruit_istk + 0.13*RR_diff_fruit_hstk
RR_diff_veg_tstk = 0.87*RR_diff_veg_istk + 0.13*RR_diff_veg_hstk

# Combine effects by multiplying RRs
RR_diff_total_CHD <- RR_diff_bmi_ihd*RR_diff_fruit_ihd*RR_diff_veg_ihd
RR_diff_total_Stroke <- RR_diff_bmi_tstk*RR_diff_fruit_tstk*RR_diff_veg_tstk

fixedvar<-c("Female", "Race", "DM_parent",  "HPT_Txt",  "Smoking")
keyinvars<-c("SEQN","Age_cycle","Glucose","BMI","Total_Chol","HDL","Trig","SBP","DBP","Diabetes","CVD_history")

#Updating the disease risk
raw.input.data<-matrix(as.numeric(datain_t[,keyinvars]), ncol=length(keyinvars))
colnames(raw.input.data) <- c("SEQN","Age","Glucose","BMI","Total_Chol","HDL","Trig","SBP","DBP","Diabetes","CVD_history")
raw.input.data <- as.data.frame(cbind(raw.input.data, data_for_analysis[,fixedvar]))

if (t%/%2 > 0 & t%%2 == 0)  {
  CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
  datain_t1[,"CVD_recurrent_prob"] <- Multi_yr_Risk_to_annual_prob(time=2, risk=CVD_Recurrent_risk_2yr)
} else {
  datain_t1[,"CVD_recurrent_prob"] <- datain_t[,"CVD_recurrent_prob"]
}

if (t%/%8 > 0 & t%%8 == 0)  {
  DM_risk_8yr <- calc_DM_risk(raw.input.data)
  datain_t1[,"DM_prob"] <- Multi_yr_Risk_to_annual_prob(time=8, risk=DM_risk_8yr)
} else {
  datain_t1[,"DM_prob"] <- datain_t[,"DM_prob"]
}

if (t%/%10 > 0 & t%%10 == 0) {
  ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
  datain_t1[,"CVD_prob"] <- Multi_yr_Risk_to_annual_prob(time=10, risk=ASCVD_Risk_10yr)
} else {
  datain_t1[,"CVD_prob"] <- datain_t[,"CVD_prob"]
}

#Defining transition probabilities
# Subtract 19 from age because indexing starts at age 20
p.death.DM <- case_when(
  data_for_analysis[,"DEMO"] == "Male" ~ DM_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "Female" ~ DM_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWM" ~ DM_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWF" ~ DM_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBM" ~ DM_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBF" ~ DM_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HM" ~ DM_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HF" ~ DM_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s]
)

p.death.CHD <- case_when(
  data_for_analysis[,"DEMO"] == "Male" ~ CHD_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "Female" ~ CHD_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWM" ~ CHD_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWF" ~ CHD_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBM" ~ CHD_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBF" ~ CHD_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HM" ~ CHD_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HF" ~ CHD_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s]
)

p.death.Stroke <- case_when(
  data_for_analysis[,"DEMO"] == "Male" ~ stroke_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "Female" ~ stroke_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWM" ~ stroke_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWF" ~ stroke_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBM" ~ stroke_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBF" ~ stroke_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HM" ~ stroke_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HF" ~ stroke_mortality_HF[as.numeric(datain_t[,"Age_cycle"])-19,s]
)

p.death <- case_when(
  data_for_analysis[,"DEMO"] == "Male" ~ Non_CVD_DM_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "Female" ~ Non_CVD_DM_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWM" ~ Non_CVD_DM_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHWF" ~ Non_CVD_DM_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBM" ~ Non_CVD_DM_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "NHBF" ~ Non_CVD_DM_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HM" ~ Non_CVD_DM_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
  data_for_analysis[,"DEMO"] == "HF" ~ Non_CVD_DM_mortality_HF[as.numeric(datain_t[,"Age_cycle"])-19,s]
)


#Disaggregating ASCVD
p.CHD_first <- as.numeric(Prop_CHD[1,data_for_analysis[,"DEMO"]]) #Prob. of CHD event (Angina, MI, Fatal MI, Fatal CHD) Source: Benjamin et al, Circulation, 2018 Table 13-1 & 18-1, and 18-2
p.CHD_recurrent <- as.numeric(Prop_CHD[2,data_for_analysis[,"DEMO"]])

#Markov State #1: "No CVD, No Diabetes"
p.H.2.DM <- as.numeric(datain_t[,"DM_prob"])*data_for_analysis$risk_adjustment.DM
p.H.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(datain_t[,"CVD_prob"])
p.H.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(datain_t[,"CVD_prob"])
p.H.2.H <- ifelse((p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death) > 1, 0,
                  1 - (p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death))

#Markov State #2: "No CVD, With Diabetes
p.DM.2.death <- (1-exp(-(p.death+p.death.DM)))
p.DM.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(datain_t[,"CVD_prob"])
p.DM.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(datain_t[,"CVD_prob"])
p.DM.2.DM <- ifelse((p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death) > 1, 0,
                    1 - (p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death))

#Markov State #3: "First Stroke"
p.initial_Stroke.2.death <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke))))
p.initial_Stroke.2.CVD_No_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 0, 1 - p.initial_Stroke.2.death)
p.initial_Stroke.2.CVD_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 1 - p.initial_Stroke.2.death, 0)

#Markov State #4: "First CHD w/o RVSC" 
p.initial_CHD_No_RVSC.2.death <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD))))
p.initial_CHD_No_RVSC.2.CVD_No_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 0, 1 - p.initial_CHD_No_RVSC.2.death)
p.initial_CHD_No_RVSC.2.CVD_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 1 - p.initial_CHD_No_RVSC.2.death, 0)

#Markov State #5: "First CHD with RVSC"
p.initial_CHD_RVSC.2.death <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s]))))
p.initial_CHD_RVSC.2.CVD_No_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 0, 1 - p.initial_CHD_RVSC.2.death)
p.initial_CHD_RVSC.2.CVD_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 1 - p.initial_CHD_RVSC.2.death, 0)

#Markov State #6: "CVD History, No Diabetes"
p.CVD.2.DM <- as.numeric(datain_t[,"DM_prob"])*data_for_analysis$risk_adjustment.DM
p.CVD.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(datain_t[,"CVD_recurrent_prob"])
p.CVD.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(datain_t[,"CVD_recurrent_prob"])
p.CVD.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke)))
p.CVD.2.CVD <- ifelse((p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death) > 1, 0,
                      1 - (p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death))

#Markov State #7: "CVD History, With Diabetes"
p.CVD_DM.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(datain_t[,"CVD_recurrent_prob"])
p.CVD_DM.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(datain_t[,"CVD_recurrent_prob"])
p.CVD_DM.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke+p.death.DM)))
p.CVD_DM.2.CVD_DM <- ifelse((p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke) > 1, 0,
                            1 - (p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke))

#Markov State #8: "Subsequent Stroke"
p.sub_Stroke.2.death <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke))))
p.sub_Stroke.2.CVD_No_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 0, 1 - p.sub_Stroke.2.death)
p.sub_Stroke.2.CVD_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 1 - p.sub_Stroke.2.death, 0)

#Markov State #9: "Subsequent CHD w/o RVSC" 
p.sub_CHD_No_RVSC.2.death <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD))))
p.sub_CHD_No_RVSC.2.CVD_No_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 0, 1 - p.sub_CHD_No_RVSC.2.death)
p.sub_CHD_No_RVSC.2.CVD_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 1 - p.sub_CHD_No_RVSC.2.death, 0)

#Markov State #10: "Subsequent CHD with RVSC"
p.sub_CHD_RVSC.2.death <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s]))))
p.sub_CHD_RVSC.2.CVD_No_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 0, 1 - p.sub_CHD_RVSC.2.death)
p.sub_CHD_RVSC.2.CVD_DM <- ifelse(as.numeric(datain_t[,"Diabetes"]) ==1, 1 - p.sub_CHD_RVSC.2.death, 0)

#Assign transition probablities
p.transition <- array(NA, dim=c(n.individual, n.health.state),
                      dimnames = list(data_for_analysis$Subject_ID,name.health.state))
p.transition[,"No CVD, No Diabetes"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.H.2.H, rep(0, n.individual))
p.transition[,"No CVD, With Diabetes"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.H.2.DM, 
                                                 ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", p.DM.2.DM, rep(0, n.individual)))
p.transition[,"First Stroke"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.H.2.initial_Stroke, 
                                        ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", p.DM.2.initial_Stroke, rep(0, n.individual)))
p.transition[,"First CHD w/o RVSC"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", (1-p.RVSC)*p.H.2.initial_CHD, 
                                              ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", (1-p.RVSC)*p.DM.2.initial_CHD, rep(0, n.individual)))
p.transition[,"First CHD with RVSC"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.RVSC*p.H.2.initial_CHD, 
                                               ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", p.RVSC*p.DM.2.initial_CHD, rep(0, n.individual)))
p.transition[,"CVD History, No Diabetes"] <- case_when(
  sim_out[,"state",t] == "First Stroke" ~ p.initial_Stroke.2.CVD_No_DM, 
  sim_out[,"state",t] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_No_DM, 
  sim_out[,"state",t] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_No_DM, 
  sim_out[,"state",t] == "CVD History, No Diabetes" ~ p.CVD.2.CVD, 
  sim_out[,"state",t] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_No_DM, 
  sim_out[,"state",t] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_No_DM, 
  sim_out[,"state",t] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_No_DM, 
  TRUE ~ rep(0, n.individual)
)
p.transition[,"CVD History, With Diabetes"] <- case_when(
  sim_out[,"state",t] == "First Stroke" ~ p.initial_Stroke.2.CVD_DM, 
  sim_out[,"state",t] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_DM, 
  sim_out[,"state",t] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_DM, 
  sim_out[,"state",t] == "CVD History, No Diabetes" ~ p.CVD.2.DM, 
  sim_out[,"state",t] == "CVD History, With Diabetes" ~ p.CVD_DM.2.CVD_DM, 
  sim_out[,"state",t] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_DM, 
  sim_out[,"state",t] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_DM, 
  sim_out[,"state",t] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_DM, 
  TRUE ~ rep(0, n.individual)
)

p.transition[,"Subsequent Stroke"] <- ifelse(sim_out[,"state",t] == "CVD History, No Diabetes", p.CVD.2.Sub_Stroke,
                                             ifelse(sim_out[,"state",t] == "CVD History, With Diabetes", p.CVD_DM.2.Sub_Stroke, rep(0, n.individual)))
p.transition[,"Subsequent CHD w/o RVSC"] <- ifelse(sim_out[,"state",t] == "CVD History, No Diabetes", (1-p.RVSC)*p.CVD.2.Sub_CHD, 
                                                   ifelse(sim_out[,"state",t] == "CVD History, With Diabetes", (1-p.RVSC)*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
p.transition[,"Subsequent CHD with RVSC"] <- ifelse(sim_out[,"state",t] == "CVD History, No Diabetes", p.RVSC*p.CVD.2.Sub_CHD, 
                                                    ifelse(sim_out[,"state",t] == "CVD History, With Diabetes", p.RVSC*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
p.transition[,"Death"] <- case_when(
  sim_out[,"state",t] == "No CVD, No Diabetes" ~ p.death, 
  sim_out[,"state",t] == "No CVD, With Diabetes" ~ p.DM.2.death, 
  sim_out[,"state",t] == "First Stroke" ~ p.initial_Stroke.2.death, 
  sim_out[,"state",t] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.death, 
  sim_out[,"state",t] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.death, 
  sim_out[,"state",t] == "CVD History, No Diabetes" ~ p.CVD.2.death, 
  sim_out[,"state",t] == "CVD History, With Diabetes" ~ p.CVD_DM.2.death, 
  sim_out[,"state",t] == "Subsequent Stroke" ~ p.sub_Stroke.2.death, 
  sim_out[,"state",t] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.death, 
  sim_out[,"state",t] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.death, 
  TRUE ~ rep(1, n.individual)
)

# Check that all transition probabilities add to 1 (rounded to 3 digits)
# if (sum(round(rowSums(p.transition),3)!=1) != 0) {
#   p_sums = round(rowSums(p.transition),3)
#   error_out = sim_out[p_sums!=1,"state",t] # Output state of person(s) with error
#   stop("Transition probabilities do not add to 1. ", paste("Simulation", s, ", Time", t, ". "), "State(s) with error: ", error_out)
# }
# Transition to the next health state 
sim_out[,"state",t+1] <- apply(p.transition, 1, function(x) sample(name.health.state, 1, prob = x))

datain_t1[,"Diabetes"] <- ifelse(sim_out[,"state",t+1] %in% c("No CVD, With Diabetes", "CVD History, With Diabetes"), 1, 
                                 ifelse(sim_out[,"state",t+1] == "Death", NA, 
                                        ifelse(as.numeric(data_for_analysis[,"Glucose"]) > 126, 1, datain_t[,"Diabetes"])))

datain_t1[,"CVD_history"] <- ifelse(sim_out[,"state",t] %in% c("First Stroke", "First CHD w/o RVSC", "First CHD with RVSC", "CVD History, No Diabetes", "CVD History, With Diabetes",
                                                               "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), 1, 
                                    ifelse(sim_out[,"state",t+1] == "Death", NA, datain_t[,"CVD_history"]))


datain_t1[,"BMI"] <- ifelse(sim_out[,"state",t+1] == "Death", NA, datain_t1[,"BMI"])
datain_t1[,"Obesity"] <- ifelse(as.numeric(datain_t1[,"BMI"]) < 30, 0, 
                                ifelse((as.numeric(datain_t1[,"BMI"]) >= 30), 1, NA))


keyvars2 <- c("SEQN","Age_cycle","SBP", "DBP",  "BMI","Diabetes","CVD_history")
raw.input.data<-matrix(as.numeric(datain_t1[,keyvars2]), ncol=length(keyvars2))
colnames(raw.input.data)<-c("SEQN","Age","SBP", "DBP",  "BMI","Diabetes","CVD_history")
raw.input.data <- as.data.frame(cbind(raw.input.data, data_for_analysis[,c("Female", "Race","HPT_Txt")]))

# Update QALYs and costs
sim_out[,"HRQOL_scores",t+1] <- calc_HRQOL(raw.input.data, HRQOL_parameter_sim[,s])
sim_out[,"HRQOL_scores",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "Subsequent Stroke"), as.numeric(sim_out[,"HRQOL_scores",t+1])+u_stroke_sim[s], 
                                       ifelse(sim_out[,"state",t+1] %in% c("First CHD w/o RVSC", "First CHD with RVSC", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), as.numeric(sim_out[,"HRQOL_scores",t+1])+u_CHD_sim[s],
                                              ifelse(sim_out[,"state",t+1] == "Death", NA, sim_out[,"HRQOL_scores",t+1])))

#ppp: for sensitivity analysis, add economic effect of HA1c


sim_out[,"HCE_predict",t+1] <- calc_HCE(raw.input.data, HCE_parameter_sim[,s])
sim_out[,"HCE_predict",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "Subsequent Stroke"), as.numeric(sim_out[,"HCE_predict",t+1])+c_stroke_sim[s], 
                                      ifelse(sim_out[,"state",t+1] %in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC"), as.numeric(sim_out[,"HCE_predict",t+1])+c_CHD_sim[s], 
                                             ifelse(sim_out[,"state",t+1] %in% c("First CHD with RVSC", "Subsequent CHD with RVSC"), as.numeric(sim_out[,"HCE_predict",t+1])+c_CHD_sim[s]+c_RVSC_sim[s],
                                                    ifelse(sim_out[,"state",t+1] == "Death", NA, sim_out[,"HCE_predict",t+1]))))

#ppp Add productivity costs associated with stroke, and CHD respectively 
sim_out[,"Prod_cost",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "Subsequent Stroke"),c_prod_stroke_sim[s], 
                                    ifelse(sim_out[,"state",t+1] %in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC","First CHD with RVSC", "Subsequent CHD with RVSC"), c_prod_CHD_sim[s],
                                           ifelse(sim_out[,"state",t+1] %in% c("CVD History, No Diabetes", "CVD History, With Diabetes"), c_prod_CVD_sim[s],      
                                                  ifelse(sim_out[,"state",t+1] == "Death", NA, 0))))

if (intervention == "Policy") {
  sim_out[,"Food_cost",t+1] <- ifelse(sim_out[,"state",t+1] != "Death", c_policy, 0)
} else {
  sim_out[,"Food_cost",t+1] <- 0
}

# Discount effects and costs
sim_out[,"effect_disc",t+1] <- as.numeric(sim_out[,"HRQOL_scores",t+1])/((1+beta_QALY)^(t-1))
sim_out[,"HCE_disc",t+1] <- as.numeric(sim_out[,"HCE_predict",t+1])/((1+beta_cost)^(t-1))
sim_out[,"Food_cost_disc",t+1] <- as.numeric(sim_out[,"Food_cost",t+1])/((1+beta_cost)^(t-1))
#ppp
sim_out[,"Prod_cost_disc",t+1] <- as.numeric(sim_out[,"Prod_cost",t+1])/((1+beta_cost)^(t-1))
#ppp: add adiministrative costs to the model

if (t == 1) {
  sim_out[,"Admin_costs",t+1] <- ifelse(sim_out[,"state",t+1] != "Death",  as.numeric(sim_out[,"Food_cost_disc",t+1])*0.2, 0)
} else {
  sim_out[,"Admin_costs",t+1] <- ifelse(sim_out[,"state",t+1] != "Death",  as.numeric(sim_out[,"Food_cost_disc",t+1])*0.05, 0)
}

datain_t=datain_t1
  }
  
  ###################3Summarizing output #####################################
  
  # Subset output to include variables of interest
  out_variables <- c("SEQN", "effect_disc",  "HCE_disc", "Food_cost_disc", "Prod_cost_disc", "Admin_costs",  "state")
  sim_out <- sim_out[,out_variables,]
  #ppp
  vars <- c( "effect_disc",  "HCE_disc",  "Food_cost_disc", "Prod_cost_disc", "Admin_costs",  "Life Years",  "Incident First CVD", "Incident Recurrent CVD","Incident CVD","Total_cost_health","Total_cost_societ" )
  
  sim_out_subset = array(NA,dim = c( n.individual, length(vars)),
                         dimnames = list(data_for_analysis$Subject_ID, vars ))
 
  #inflate healthcare costs 
  for (i_var in c("effect_disc",  "HCE_disc", "Food_cost_disc", "Prod_cost_disc", "Admin_costs")) {
    if (i_var== "HCE_disc") {
      temp = array(as.numeric(sim_out[,i_var,c(2:n.cycle+1)])*1.062, dim = c(n.individual,n.cycle))
    }
    else {
      temp = array(as.numeric(sim_out[,i_var,c(2:n.cycle+1)]), dim = c(n.individual,n.cycle))
    }  
    sim_out_subset[, i_var]<- rowSums(temp, na.rm=T)
  }
  sim_out_subset[, "Total_cost_health"]<-sim_out_subset[,"Food_cost_disc"]+ sim_out_subset[,"Admin_costs"]+ sim_out_subset[,"HCE_disc"]
  sim_out_subset[, "Total_cost_societ"]<-sim_out_subset[,"Food_cost_disc"]+ sim_out_subset[,"Admin_costs"]+ sim_out_subset[,"HCE_disc"]+sim_out_subset[,"Prod_cost_disc"]
  
  # We will need state-specific output, so we'll subset sim_out here.
  m.state <- sim_out[,"state",c(2:(n.cycle+1))]
  rm(sim_out)
  sim_out_subset[,"Life Years"] = rowSums(m.state!="Death", na.rm = T)
  #sim_out_subset[,"Incident First Stroke"] = rowSums(m.state=="First Stroke",na.rm = T)
  #sim_out_subset[,"Incident First CHD"] = rowSums(m.state=="First CHD with RVSC" | m.state=="First CHD w/o RVSC",na.rm = T)
  #sim_out_subset[,"Incident Recurrent Stroke"] = rowSums(m.state=="Subsequent Stroke", na.rm = T)
  #sim_out_subset[,"Incident Recurrent CHD"] = rowSums(m.state=="Subsequent CHD with RVSC"|m.state=="Subsequent CHD w/o RVSC",  na.rm = T)
  sim_out_subset[,"Incident First CVD"] = rowSums(m.state=="First Stroke" |m.state=="First CHD with RVSC" | m.state=="First CHD w/o RVSC", na.rm = T)
  sim_out_subset[,"Incident Recurrent CVD"] = rowSums(m.state=="Subsequent Stroke"|m.state=="Subsequent CHD with RVSC"|m.state=="Subsequent CHD w/o RVSC",  na.rm = T)
  sim_out_subset[,"Incident CVD"] = rowSums(m.state=="First Stroke" |m.state=="First CHD with RVSC" | m.state=="First CHD w/o RVSC"|m.state=="Subsequent Stroke"|m.state=="Subsequent CHD with RVSC"|m.state=="Subsequent CHD w/o RVSC", na.rm = T)

  return(sim_out_subset)
}


