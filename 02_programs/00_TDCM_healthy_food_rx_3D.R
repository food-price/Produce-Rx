#Tufts Diabetes-CVD Microsimulation (TDCM) Model 
#The FOOD-PRICE project (https://food-price.org/)
#Healthy Food Prescription Analysis
#Primary Author: David D. Kim, Tufts Medical Center
#Contact: DKim3@tuftsmedicalcenter.org 
#Secondary Author: Brianna Lauren
#Contact: brianna.lauren@tufts.edu
#Main Purpose: 
#1) To develop a microsimulation model to link individual risk factors with cardiometabolic diseases, including diabets and CVD
#2) To generate estimations of dietary intakes and changes in risk factors and disease risk from dietary changes

#Date: May 21, 2020

#Updates: Fixing the mortality beyond age 101 so that the model runs really long enough. 
#Updates (12/20/2019): Adding DM adjustment factors for non-whites (2.4 for NHBW & HW and 1.5 for NHBM & HM) [Source: Brancati et al, JAMA, 2000 & Narayan et al, JAMA, 2003]
#Updates (01/16/2020): Complete the full CEA Model with the dietary and policy modules
#Updates (02/12/2020): Full probablistic models
#Updates (03/12/2020): Improving Efficiency of models
#Updates (04/10/2020): Summarizing Outputs
#Updates (05/21/2020): Modified for use on HPC: seed now read from CL args, removed Windows-specific components
#                       Improved file I/O speed by replacing read.csv/write.csv with data.table::fread/fwrite
#Updates (07/16/2020): Extract race-stratified outputs

#################
# 0.Preparation #
#################

# 0.1 remove all data from memory
rm(list = ls())
#memory.size(max=T) # this is Windows-specific and will not work on the HPC running Linux

# 0.2 Start the clock
ptm <- proc.time()

# 0.3 Install/Read R Packages


library(survey)
library(svMisc)
library(psych)
library(gdata)
library(dplyr)
library(data.table)   # needed for fread/fwrite
library(foreach)
library(doParallel)
library(abind)

# 0.4 Create Functions 

# 0.4-1 To estimate gamma/beta parameters based on Mean and SD
estGammaParams <- function(mu, var) {
  beta <- var / mu
  alpha <- mu / beta
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

calc_nsims_rgamma <- function(n.sim, mu, se) {
  gamma_par <- estGammaParams(mu, se^2)
  return(rgamma(n.sim, gamma_par$alpha, 1/gamma_par$beta))
}

calc_nsims_rbeta <- function(n.sim, mu, se) {
  beta_par <- estBetaParams(mu, se^2)
  return(rbeta(n=n.sim, beta_par$alpha, beta_par$beta))
}

# 0.4-2 Converting multi-year risk (e.g., 10-year ASCVD risk) to annual probabilities
Multi_yr_Risk_to_annual_prob <- function(time, risk) {
  annual_rate <- -(1/time)*log(1-risk)
  annual_prob <- 1 - exp(-annual_rate)
}

# 0.5 Creating Working Directory 
#setwd("~/Documents/Friedman/Food-PRICE/healthy_food_rx")
setwd("/cluster/tufts/kimlab/lwang18/healthy_food_rx")

# Source other scripts
source("02_Programs/1a - Diabetes_risk_prediction_FHS (categorical points).R")
source("02_Programs/1b - ASCVD_risk_calculator.R")
source("02_Programs/1c - FHS Subsequent_CVD_risk_calculator.R")
source("02_Programs/2a - HrQOL estimator for US general population.R")
source("02_Programs/2b - HCE estimator for US general population.R")

# 0.6 Model settings (set manually or read from command line)
args <- commandArgs(trailingOnly = TRUE)  # get all arguments after script.R
# If no arguments are read from command line, set modeling choices manually. 
# Otherwise, read modeling choices from command line.
if (length(args) == 0) {
  seed <- 1234
  n.sim <- 5 #Number of probablistic samplings
  n.cycle <- 10 #Number of Cycle Length: How long does the model run (i.e., analytic time horizon)
  n.sample <- 100 #Number of individuals to be selected from the full sample; if full sample, enter "ALL"
  intervention <- "Policy" #SELECT between "Policy" or "No_Policy"
} else {
  # expecting 4 arguments
  if (length(args) != 4) {
    stop("ERROR: Incorrect number of command line arguments", call. = FALSE)
  }
  seed <- as.numeric(args[1]) # extracting first arg as seed and attempt to cast to number
  n.sim <- as.numeric(args[2])
  n.cycle <- as.numeric(args[3])
  intervention <- args[4]
  # Assume entire sample
  n.sample = "ALL"
  
  # check that modeling choices were set
  if (is.na(seed)) {
    stop("ERROR: missing seed", call. = FALSE)
  }
  if (is.na(n.sim)) {
    stop("ERROR: missing n.sim", call. = FALSE)
  }
  if (is.na(n.cycle)) {
    stop("ERROR: missing n.cycle", call. = FALSE)
  }
  if (!intervention %in% c("No_Policy", "Policy")) {
    stop("ERROR: invalid intervention", call. = FALSE)
  }
}

# 0.7 Complete Modeling Choices
beta_cost <- beta_QALY <- 0.03 #Annual discounting rate of costs and QALYs
set.seed(seed)

################################################
# 1 Importing Necessary Data Inputs + Cleaning #
################################################
print('Importing data')
# 1.1 Read in master input file
NHANES<- fread("01_Input/NHANES/NHANES_1318_Imp10_dm_fi.csv", stringsAsFactors = TRUE, data.table = FALSE)

# 1.2 Select only necessary variables from master input file
variables <- c("SEQN", "WTINT2YR", "WTMEC2YR", "WT_TOTAL", "SDMVPSU", "SDMVSTRA", 
               "Age", "Female", "Race", "CVD_history", "Diabetes",
               "Total_Chol","HDL", "SBP", "DBP", "HPT_Txt", "Smoking",
               "DM_family", "BMI", "height", "weight", "Glucose", "Trig")

NHANES <- NHANES[variables]

# 1.3 Select the target starting population to model (Age 40-79) + remove observations with missing 
NHANES_age40_79 <- na.omit(subset(NHANES, Age >= 40 & Age < 80)) 
NHANES_age40_79$Subject_ID <- c(1:nrow(NHANES_age40_79))

# 1.4 Define the analytic dataset to carry it forward
if (n.sample == "ALL"){
  data_for_analysis <- NHANES_age40_79
} else {
  random_sample <- sample(1:nrow(NHANES_age40_79),n.sample,replace=F)
  data_for_analysis <- NHANES_age40_79[random_sample,]
}
data_for_analysis <- data_for_analysis[order(data_for_analysis$Subject_ID), ]

# 1.5 Important other input data   

# 1.5-1 Policy-effect size (based on 9 intervention studies)
if (intervention == "Policy") {
  policy_effect_fruit <- calc_nsims_rbeta(n.sim, mu = 0.31, se = 0.18)
  policy_effect_veg <- calc_nsims_rbeta(n.sim, mu = 0.31, se = 0.18)
} else {
  policy_effect_fruit <- policy_effect_veg <- rep(0, n.sim)
}

# 1.5-2 Age-specific relative risk estimates between dietary intake and disease outcomes  
RR_diet_disease <- fread("01_Input/final_diet_RRs_MainAnalysis v11.csv", stringsAsFactors = TRUE, data.table = FALSE)
RR_diet_disease <- subset(RR_diet_disease, minage>=35)

# 1.5-3 Health-state/Event specific mortality data  
Non_CVD_DM_mortality <- fread("01_Input/Non_CVD_Non_DM_Cause_Mortality (Annual probability).csv", stringsAsFactors = TRUE, data.table = FALSE)
Non_CVD_DM_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Male, Non_CVD_DM_mortality$Male_SE))
Non_CVD_DM_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Female, Non_CVD_DM_mortality$Female_SE))
Non_CVD_DM_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWM, Non_CVD_DM_mortality$NHWM_SE))
Non_CVD_DM_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWF, Non_CVD_DM_mortality$NHWF_SE))
Non_CVD_DM_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBM, Non_CVD_DM_mortality$NHBM_SE))
Non_CVD_DM_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBF, Non_CVD_DM_mortality$NHBF_SE))
Non_CVD_DM_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HM, Non_CVD_DM_mortality$HM_SE))
Non_CVD_DM_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HF, Non_CVD_DM_mortality$HF_SE))

stroke_mortality <- fread("01_Input/Stroke_Cause_Mortality (Annual probability).csv", stringsAsFactors = TRUE, data.table = FALSE)
stroke_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Male, stroke_mortality$Male_SE))
stroke_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Female, stroke_mortality$Female_SE))
stroke_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWM, stroke_mortality$NHWM_SE))
stroke_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWF, stroke_mortality$NHWF_SE))
stroke_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBM, stroke_mortality$NHBM_SE))
stroke_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBF, stroke_mortality$NHBF_SE))
stroke_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HM, stroke_mortality$HM_SE))
stroke_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HF, stroke_mortality$HF_SE))
  
CHD_mortality <- fread("01_Input/IHD_Cause_Mortality (Annual probability).csv", stringsAsFactors = TRUE, data.table = FALSE)
CHD_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Male, CHD_mortality$Male_SE))
CHD_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Female, CHD_mortality$Female_SE))
CHD_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWM, CHD_mortality$NHWM_SE))
CHD_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWF, CHD_mortality$NHWF_SE))
CHD_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBM, CHD_mortality$NHBM_SE))
CHD_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBF, CHD_mortality$NHBF_SE))
CHD_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HM, CHD_mortality$HM_SE))
CHD_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HF, CHD_mortality$HF_SE))

DM_mortality <- fread("01_Input/DM_Cause_Mortality (Annual probability).csv", stringsAsFactors = TRUE, data.table = FALSE)
DM_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Male, DM_mortality$Male_SE))
DM_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Female, DM_mortality$Female_SE))
DM_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWM, DM_mortality$NHWM_SE))
DM_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWF, DM_mortality$NHWF_SE))
DM_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBM, DM_mortality$NHBM_SE))
DM_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBF, DM_mortality$NHBF_SE))
DM_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HM, DM_mortality$HM_SE))
DM_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HF, DM_mortality$HF_SE))

# 1.5-4 Secular trends in major risk factors: estimated as average annual percent change over 1999-2016 
risk_factor_trend <- fread("01_Input/Risk_Factors_Trends_Summary (Final).csv", stringsAsFactors = TRUE, data.table = FALSE)
risk_factor_trend_sim <- matrix(NA, nrow=nrow(risk_factor_trend), ncol=n.sim)
colnames(risk_factor_trend_sim) <- paste("simulation_", 1:n.sim, sep = " ")

risk_factor_trend_sim <- t(mapply(rnorm, n=n.sim, risk_factor_trend$APC, risk_factor_trend$APC_SE))

for (i in 1:nrow(risk_factor_trend)) {
  risk_factor_trend_sim[i,] <- as.matrix(rnorm(n.sim, risk_factor_trend$APC[i], risk_factor_trend$APC_SE[i]))
}

risk_factor_trend_sim <- cbind(risk_factor_trend, risk_factor_trend_sim)

# 1.5-5 Gender- and race-specific proportiona of CHD cases Among ASCVD cases (CHD + Stroke)
Prop_CHD <- fread("01_Input/Prop_CHD.csv", stringsAsFactors = FALSE, data.table = FALSE)

#1.5-6 Health Care Expenditure Model
HCE_parameter_sim <- matrix(NA, nrow=nrow(HCE_parameters), ncol=n.sim)
rownames(HCE_parameter_sim) <- rownames(HCE_parameters)
colnames(HCE_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

for (i in 1:nrow(HCE_parameters)) {
  HCE_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HCE_parameters$Beta[i], HCE_parameters$SE[i]))
}

#1.5-7 HrQOL Prediction
HRQOL_parameter_sim <- matrix(NA, nrow=nrow(HRQOL_parameters), ncol=n.sim)
rownames(HRQOL_parameter_sim) <- rownames(HRQOL_parameters)
colnames(HRQOL_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

for (i in 1:nrow(HRQOL_parameters)) {
  HRQOL_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HRQOL_parameters$Beta[i], HRQOL_parameters$SE[i]))
}

# 1.6 DATA CLEANING

data_for_analysis$initial_H[data_for_analysis$CVD_history == 0 & data_for_analysis$Diabetes == 0] <- "No CVD, No Diabetes"
data_for_analysis$initial_H[data_for_analysis$CVD_history == 0 & data_for_analysis$Diabetes == 1] <- "No CVD, With Diabetes"
data_for_analysis$initial_H[data_for_analysis$CVD_history == 1 & data_for_analysis$Diabetes == 0] <- "CVD History, No Diabetes"
data_for_analysis$initial_H[data_for_analysis$CVD_history == 1 & data_for_analysis$Diabetes == 1] <- "CVD History, With Diabetes"

data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 1] <- "NHWM"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 1] <- "NHWF"
data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 2] <- "NHBM"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 2] <- "NHBF"
data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 3] <- "HM"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 3] <- "HF"
data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 4] <- "Male"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 4] <- "Female"

data_for_analysis$Age_cycle <- NA

data_for_analysis$DM_parent <- data_for_analysis$DM_family

data_for_analysis$BMI_cat[data_for_analysis$BMI < 18.5] <- "Underweight"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 18.5 & data_for_analysis$BMI < 25] <- "Normal"
data_for_analysis$BMI_cat[data_for_analysis$BMI >=25 & data_for_analysis$BMI < 30] <- "Overweight"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 30 & data_for_analysis$BMI < 35] <- "Obesity I"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 35 & data_for_analysis$BMI < 40] <- "Obesity II"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 40] <- "Obesity III"

data_for_analysis$Obesity[data_for_analysis$BMI < 30] <- 0
data_for_analysis$Obesity[data_for_analysis$BMI >= 30] <- 1

# 1.7 DM risk adjustment for non-whites
data_for_analysis$risk_adjustment.DM <- ifelse(data_for_analysis$DEMO %in% c("NHWM", "NHWF", "Female", "Male"), 1.0,
                                               ifelse(data_for_analysis$DEMO %in% c("NHBM", "NHBM"), 1.5, 2.4))

# 1.7 Creating dupilicate (counter-factual) observations to predict the effect under policy vs. no policy
data_for_analysis$source <- intervention

# 1.8 Diet-disease data inputs 

# 1.8.1 Indirect effect on disease
# Weight change due to change in F&V intake (4-year weight change converted to yearly)
wt_change_fruit = rnorm(n.sim, -0.53, 0.046)/4
wt_change_veg = rnorm(n.sim, -0.25, 0.051)/4

# 1.8.2 Direct effect on disease
RR_fruit_hstk <- subset(RR_diet_disease, outcome == 'HSTK' & riskfactor == 'fruit', select = c(minage, logRR.perchange, se.perchange))
RR_fruit_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'fruit', select = c(minage, logRR.perchange, se.perchange))
RR_fruit_istk <- subset(RR_diet_disease, outcome == 'ISTK' & riskfactor == 'fruit', select = c(minage, logRR.perchange, se.perchange))

RR_veg_hstk <- subset(RR_diet_disease, outcome == 'HSTK' & riskfactor == 'veg', select = c(minage, logRR.perchange, se.perchange))
RR_veg_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'veg', select = c(minage, logRR.perchange, se.perchange))
RR_veg_istk <- subset(RR_diet_disease, outcome == 'ISTK' & riskfactor == 'veg', select = c(minage, logRR.perchange, se.perchange))

RR_bmi_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'bmi', select = c(minage, logRR.perchange, se.perchange))
RR_bmi_tstk <- subset(RR_diet_disease, outcome == 'TSTK' & riskfactor == 'bmi', select = c(minage, logRR.perchange, se.perchange))

# Prepare matrix
random_logrr <- matrix(NA, nrow=nrow(RR_fruit_hstk), ncol=n.sim)
rownames(random_logrr) <- c("35-44 y", "45-54 y", "55-64 y", "65-74 y","75+ y")
colnames(random_logrr) <- paste("simulation_", 1:n.sim, sep = " ")

random_logrr_fruit_hstk <- random_logrr_fruit_ihd <- random_logrr_fruit_istk <-
  random_logrr_veg_hstk <- random_logrr_veg_ihd <- random_logrr_veg_istk <- 
  random_logrr_BMI_IHD <- random_logrr_BMI_TSTK <- random_logrr

for (g in 1:nrow(random_logrr)) {
  random_logrr_fruit_hstk[g,] <- rnorm(n.sim, RR_fruit_hstk$logRR.perchange[g], RR_fruit_hstk$se.perchange[g])
  random_logrr_fruit_ihd[g,] <- rnorm(n.sim, RR_fruit_ihd$logRR.perchange[g], RR_fruit_ihd$se.perchange[g])
  random_logrr_fruit_istk[g,] <- rnorm(n.sim, RR_fruit_istk$logRR.perchange[g], RR_fruit_istk$se.perchange[g])
  random_logrr_veg_hstk[g,] <- rnorm(n.sim, RR_veg_hstk$logRR.perchange[g], RR_veg_hstk$se.perchange[g])
  random_logrr_veg_ihd[g,] <- rnorm(n.sim, RR_veg_ihd$logRR.perchange[g], RR_veg_ihd$se.perchange[g])
  random_logrr_veg_istk[g,] <- rnorm(n.sim, RR_veg_istk$logRR.perchange[g], RR_veg_istk$se.perchange[g])
  random_logrr_BMI_IHD[g,] <- rnorm(n.sim, RR_bmi_ihd$logRR.perchange[g], RR_bmi_ihd$se.perchange[g])
  random_logrr_BMI_TSTK[g,] <- rnorm(n.sim, RR_bmi_tstk$logRR.perchange[g], RR_bmi_tstk$se.perchange[g])
}

#############################################
# 2 Describing the baseline characteristics #
#############################################
  #OMITTED - SEE THE PREVIOUS VERSION

#################################################################################################################################
# 3 Estimating disease-specific risk, health-related quality of life (HrQOL), and healthcare expenditures (HCE) at the baseline #
#################################################################################################################################
print('part 3')
p3_start <- proc.time()
# 3.1 FHS 8-year Diabetes Risk Prediction (With BMI)
variable_for_raw.input <- names(data_for_analysis)
raw.input.data <- data_for_analysis
data_for_analysis$DM_risk_8yr <- calc_DM_risk(raw.input.data)
data_for_analysis$DM_prob <- Multi_yr_Risk_to_annual_prob(time=8, risk=data_for_analysis$DM_risk_8yr)

# 3.2 ACC/AHA ASCVD 10-year Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
data_for_analysis$CVD_prob <- Multi_yr_Risk_to_annual_prob(time=10, risk=data_for_analysis$ASCVD_Risk_10yr)

# 3.3 FHS 2-year CVD recurrent Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
data_for_analysis$CVD_recurrent_prob <- Multi_yr_Risk_to_annual_prob(time=2, risk=data_for_analysis$CVD_Recurrent_risk_2yr)

# 3.4 Individual HrQOL prediction
raw.input.data <- data_for_analysis
data_for_analysis$HRQOL_scores <- calc_HRQOL(raw.input.data, HRQOL_parameters[,1])

# 3.5 Individual HCE prediction
raw.input.data <- data_for_analysis
data_for_analysis$HCE_predict <- calc_HCE(raw.input.data, HCE_parameters[,1])

p3_end <- proc.time()
p3_end - p3_start

###############################################################################################################################################################
# 4 Microsimulation Model - Predicting individual-levle changes in diet, risk factors, disease outcomes, health-related quality of life, and healthcare costs # #
###############################################################################################################################################################
print('part 4')
# 4.1 Modeling Input
initial_H <- data_for_analysis$initial_H #Vector of Initial states for individuals
n.individual <- nrow(data_for_analysis) # number of individuals in the model
name.health.state <- c("No CVD, No Diabetes", "No CVD, With Diabetes", "First Stroke", "First CHD w/o RVSC", "First CHD with RVSC",
                       "CVD History, No Diabetes", "CVD History, With Diabetes", "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC", "Death") 
n.health.state <- length(name.health.state) # number of health state that individuals can transition over time

#Health care costs associated with clinical ASCVD events  
#Source: https://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/Medicare-Provider-Charge-Data/Inpatient2017.html
c_CHD <- 10034.34 # Medicare Inpatient Prospective Payment Hospitals across six DRGs for MI [280-285] (Average total payments)
c_CHD_sim <- calc_nsims_rgamma(n.sim, c_CHD, c_CHD*0.2)

c_stroke <- 15994.49 # Medicare Inpatient Prospective Payment Hospitals across three DRGs for stroke [061-063] (Average total payments)
c_stroke_sim <- calc_nsims_rgamma(n.sim, c_stroke, c_stroke*0.2)

#Disutility associated with major clinical events
#Source: Davies et al. Health and Quality of Life Outcomes (2015) 13:90.DOI 10.1186/s12955-015-0266-9

u_CHD <- -0.055 #Average disutility among MI (-0.06) and unstable angina (-0.05)
u_CHD_sim <- -calc_nsims_rbeta(n.sim, -u_CHD, -u_CHD*0.2)

u_stroke <- -0.3 
u_stroke_sim <- -calc_nsims_rbeta(n.sim, -u_stroke, -u_stroke*0.2)

#RVSC-related Input
p.RVSC <- 0.67 #Probabilities of receiving RVSC after CVD event, Source: Benjamin et al, Circulation, 2018, Table 13-1[# of Stroke Discharge] & 18-1[# of CHD Hospital Discharge] & 24-2[# of inpatient RVSC]
p.RVSC_sim <- calc_nsims_rbeta(n.sim, p.RVSC, p.RVSC*0.2)

prop_PCI <- 0.711 #Proportions of PCI (including PCI with stent) among all RVSC (PI & CABG) Source: Heart Disease and Stroke Statistics (Benjamin et al., Circulation, 2018) Table 24-2
prop_PCI_sim <- calc_nsims_rbeta(n.sim, prop_PCI, prop_PCI*0.2)

c_CABG <- 44538.01 # Medicare Inpatient Prospective Payment Hospitals across six DRGs for CABG[231-236] (Average total payments)
c_CABG_sim <- calc_nsims_rgamma(n.sim, c_CABG, c_CABG*0.2)

c_PCI <- 18476.82 # Medicare Inpatient Prospective Payment Hospitals across nine DRGs for PCI[034-036; 246-251] (Average total payments)
c_PCI_sim <- calc_nsims_rgamma(n.sim, c_PCI, c_PCI*0.2)

c_RVSC <- prop_PCI*c_PCI + (1-prop_PCI)*c_CABG
c_RVSC_sim <- prop_PCI_sim*c_PCI_sim + (1-prop_PCI_sim)*c_CABG_sim

# Policy costs
c_policy = 27 * 12 # weighted average of dollar amount in 9 studies included in the meta-analysis (annual)
# c_policy_sim = calc_nsims_rgamma(n.sim, c_policy, c_policy*0.2)

mortality_CABG <- 0.0178 #In-hospital death rate, Source: Benjamin et al, Circulation, 2018, Table 24-1
mortality_CABG_sim <- calc_nsims_rbeta(n.sim, mortality_CABG, mortality_CABG*0.2)

mortality_PCI <- 0.0207 #In-hospital death rate, Source: Benjamin et al, Circulation, 2018, Table 24-1
mortality_PCI_sim <- calc_nsims_rbeta(n.sim, mortality_PCI, mortality_PCI*0.2)

p.death.RVSC <- prop_PCI*mortality_PCI + (1-prop_PCI)*mortality_CABG
p.death.RVSC_sim <- prop_PCI_sim*mortality_PCI_sim + (1-prop_PCI_sim)*mortality_CABG_sim

#########################################################################
# 5 Main Model - PLEASE DO NOT MODIFY UNLESS YOU KNOW WHAT YOU'RE DOING #
#########################################################################
print('running model')
model_start = proc.time()

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

run_sim <- function(s) {
  # Runs single simulation
  # Returns array with variables for each individual at each cycle
  key_variables <- c("SEQN", "Age", "Age_cycle", "DEMO", "Female", "Race", "Total_Chol","HDL", "SBP", "DBP", "BMI", "BMI_cat", "Obesity", "weight", "height", "HPT_Txt", "Trig", "Glucose", 
                     "Smoking", "Diabetes", "CVD_history", "risk_adjustment.DM", "DM_parent", "DM_prob", "CVD_prob", "CVD_recurrent_prob", "HRQOL_scores", "HCE_predict")  
  
  all_variables <- c(key_variables, "state", "HCE_disc", "effect_disc", "Food_cost", "Food_cost_disc")
  # Initialize array to hold output from single simulation
  sim_out <- array(NA, dim=c(n.individual, length(all_variables), n.cycle+1),
                   dimnames = list(data_for_analysis$Subject_ID,
                                   all_variables, paste("cycle", 0:n.cycle, sep = " ")))
  
 #999 for checking only 
  RR_CHD<-array(NA, dim=c(n.individual, n.cycle+1, s))
  RR_stk<-array(NA, dim=c(n.individual, n.cycle+1, s))
  
  # Set initial state (cycle 0)
  sim_out[,1:length(key_variables),1] <- as.matrix(data_for_analysis[,key_variables])
  sim_out[,"Age_cycle",1] <- sim_out[,"Age",1]
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
  
  for (t in 1:n.cycle) {
    #Non-time varying data inputs: carry it over from the baseline data
    sim_out[,c("SEQN", "Age", "Female", "Race","DEMO", "HPT_Txt","Smoking","DM_parent", "risk_adjustment.DM"),t+1] <- sim_out[,c("SEQN", "Age", "Female", "Race","DEMO", "HPT_Txt","Smoking","DM_parent", "risk_adjustment.DM"),t]
    
    #Time-varying data inputs
    sim_out[,"Age_cycle",t+1] <- as.numeric(sim_out[,"Age",1]) + t
    
    # 5.4 Updating Risk Factors over time: Applying secular trends based on age, gender, R/E
    # Adjust BMI based on change in fruit and vegetable intake (only applies in first 4 years)
    if (intervention == "Policy" & t <= 4) {
      weight = as.numeric(sim_out[,"weight",1]) + t*wt_change_fruit[s]*policy_effect_fruit[s]
      weight = weight + t*wt_change_veg[s]*policy_effect_veg[s]
      sim_out[,"BMI",t+1] = round(weight / (as.numeric(sim_out[,"height",1])/100)^2, 1)
      sim_out[,"BMI",t+1] <- update_risk_factor(sim_out[,,t+1], "BMI", "BMI", s)
    } else {
      sim_out[,"BMI",t+1] <- update_risk_factor(sim_out[,,t], "BMI", "BMI", s)
    }
    # Truncate BMI to the interval 12-70 kg/m2. Source: ???	https://doi.org/10.2105/AJPH.2008.137364  
    sim_out[,"BMI",t+1] <- ifelse(as.numeric(sim_out[,"BMI",t+1]) < 12, 12,
                                  ifelse(as.numeric(sim_out[,"BMI",t+1])>70, 70, as.numeric(sim_out[,"BMI",t+1])))
    
    sim_out[,"BMI_cat",t+1] <- case_when(
      as.numeric(sim_out[,"BMI",t+1]) < 18.5 ~ "Underweight",
      as.numeric(sim_out[,"BMI",t+1]) < 25 ~ "Normal",
      as.numeric(sim_out[,"BMI",t+1]) < 30 ~ "Overweight",
      as.numeric(sim_out[,"BMI",t+1]) < 35 ~ "Obesity I",
      as.numeric(sim_out[,"BMI",t+1]) < 40 ~ "Obesity II",
      as.numeric(sim_out[,"BMI",t+1]) >= 40 ~ "Obesity III"
    )
    
    sim_out[,"Total_Chol",t+1] <- update_risk_factor(sim_out[,,t], "Total_Chol", "Total_Chol", s)
    sim_out[,"HDL",t+1] <- update_risk_factor(sim_out[,,t], "HDL", "HDL", s)
    
    sim_out[,"SBP",t+1] <- update_risk_factor(sim_out[,,t], "SBP", "SBP", s)
    
    sim_out[,"DBP",t+1] <- update_risk_factor(sim_out[,,t], "DBP", "SBP", s)
    sim_out[,"Trig",t+1] <- update_risk_factor(sim_out[,,t], "Trig", "Trig", s)
    sim_out[,"Glucose",t+1] <- update_risk_factor(sim_out[,,t], "Glucose", "Glucose", s)
    
    # Fruit-IHD/CHD
    RR_diff_fruit_hstk <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_fruit_hstk[1,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_fruit_hstk[2,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_fruit_hstk[3,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_fruit_hstk[4,s]*policy_effect_fruit[s]),
      TRUE ~ exp(random_logrr_fruit_hstk[5,s]*policy_effect_fruit[s])
    )
    RR_diff_fruit_ihd <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_fruit_ihd[1,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_fruit_ihd[2,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_fruit_ihd[3,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_fruit_ihd[4,s]*policy_effect_fruit[s]),
      TRUE ~ exp(random_logrr_fruit_ihd[5,s]*policy_effect_fruit[s])
    )
    RR_diff_fruit_istk <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_fruit_istk[1,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_fruit_istk[2,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_fruit_istk[3,s]*policy_effect_fruit[s]),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_fruit_istk[4,s]*policy_effect_fruit[s]),
      TRUE ~ exp(random_logrr_fruit_istk[5,s]*policy_effect_fruit[s])
    )
    
    # Veg-IHD/CHD
    RR_diff_veg_hstk <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_veg_hstk[1,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_veg_hstk[2,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_veg_hstk[3,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_veg_hstk[4,s]*policy_effect_veg[s]),
      TRUE ~ exp(random_logrr_veg_hstk[5,s]*policy_effect_veg[s])
    )
    RR_diff_veg_ihd <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_veg_ihd[1,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_veg_ihd[2,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_veg_ihd[3,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_veg_ihd[4,s]*policy_effect_veg[s]),
      TRUE ~ exp(random_logrr_veg_ihd[5,s]*policy_effect_veg[s])
    )
    RR_diff_veg_istk <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_veg_istk[1,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_veg_istk[2,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_veg_istk[3,s]*policy_effect_veg[s]),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_veg_istk[4,s]*policy_effect_veg[s]),
      TRUE ~ exp(random_logrr_veg_istk[5,s]*policy_effect_veg[s])
    )
    
    # BMI-IHD/CHD
    RR_diff_bmi_ihd <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_BMI_IHD[1,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5), 
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_BMI_IHD[2,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_BMI_IHD[3,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_BMI_IHD[4,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5),  
      TRUE ~ exp(random_logrr_BMI_IHD[5,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5)
    )
    
    # BMI-Total stroke
    RR_diff_bmi_tstk <- case_when(
      sim_out[,"Age_cycle",t] < 45 ~ exp(random_logrr_BMI_TSTK[1,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5), 
      sim_out[,"Age_cycle",t] < 55 ~ exp(random_logrr_BMI_TSTK[2,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5),
      sim_out[,"Age_cycle",t] < 65 ~ exp(random_logrr_BMI_TSTK[3,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5),
      sim_out[,"Age_cycle",t] < 75 ~ exp(random_logrr_BMI_TSTK[4,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5),
      TRUE ~ exp(random_logrr_BMI_TSTK[5,s]*(as.numeric(sim_out[,"BMI",t+1]) - as.numeric(sim_out[,"BMI",1]))/5)
    )
    
    # Calculate weighted average of HSTK and ISTK risk ratios to get risk ratio for total stroke
    # 87% of strokes are ischemic (Source: Benjamin et al. Heart Disease and Stroke Statistics—2019 Update: A Report From the American Heart Association)
    RR_diff_fruit_tstk = 0.87*RR_diff_fruit_istk + 0.13*RR_diff_fruit_hstk
    RR_diff_veg_tstk = 0.87*RR_diff_veg_istk + 0.13*RR_diff_veg_hstk
    
    # Combine effects by multiplying RRs
    RR_diff_total_CHD <- RR_diff_bmi_ihd*RR_diff_fruit_ihd*RR_diff_veg_ihd
    RR_diff_total_Stroke <- RR_diff_bmi_tstk*RR_diff_fruit_tstk*RR_diff_veg_tstk
    
  
    #ppp for checking only 
    
    RR_CHD[,t,s]<- RR_diff_total_CHD 
    RR_stk[,t,s]<-  RR_diff_total_Stroke
    

    #Updating the disease risk
    raw.input.data  <- cbind(as.numeric(sim_out[,"SEQN",t]), as.numeric(sim_out[,"Age_cycle",t]), as.numeric(sim_out[,"Female",t]), as.numeric(sim_out[,"Race",t]),
                             as.numeric(sim_out[,"Glucose",t]), as.numeric(sim_out[,"BMI",t]), as.numeric(sim_out[,"Total_Chol",t]), 
                             as.numeric(sim_out[,"HDL",t]), as.numeric(sim_out[,"Trig",t]), as.numeric(sim_out[,"SBP",t]),as.numeric(sim_out[,"DBP",t]),
                             as.numeric(sim_out[,"HPT_Txt",t]), as.numeric(sim_out[,"Smoking",t]), as.numeric(sim_out[,"Diabetes",t]),  as.numeric(sim_out[,"CVD_history",t]), as.numeric(sim_out[,"DM_parent",t]))
    colnames(raw.input.data) <- c("SEQN","Age","Female","Race","Glucose","BMI","Total_Chol","HDL","Trig","SBP","DBP","HPT_Txt","Smoking","Diabetes","CVD_history", "DM_parent" )
    raw.input.data <- as.data.frame(raw.input.data)
    
    if (t%/%2 > 0 & t%%2 == 0)  {
      CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
      sim_out[,"CVD_recurrent_prob",t+1] <- Multi_yr_Risk_to_annual_prob(time=2, risk=CVD_Recurrent_risk_2yr)
    } else {
      sim_out[,"CVD_recurrent_prob",t+1] <- sim_out[,"CVD_recurrent_prob",t]
    }
    
    if (t%/%8 > 0 & t%%8 == 0)  {
      DM_risk_8yr <- calc_DM_risk(raw.input.data)
      sim_out[,"DM_prob",t+1] <- Multi_yr_Risk_to_annual_prob(time=8, risk=DM_risk_8yr)
    } else {
      sim_out[,"DM_prob",t+1] <- sim_out[,"DM_prob",t]
    }
    
    if (t%/%10 > 0 & t%%10 == 0) {
      ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
      sim_out[,"CVD_prob",t+1] <- Multi_yr_Risk_to_annual_prob(time=10, risk=ASCVD_Risk_10yr)
    } else {
      sim_out[,"CVD_prob",t+1] <- sim_out[,"CVD_prob",t]
    }
    
    #Defining transition probabilities
    # Subtract 19 from age because indexing starts at age 20
    p.death.DM <- case_when(
      sim_out[,"DEMO",t] == "Male" ~ DM_mortality_Male[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "Female" ~ DM_mortality_Female[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWM" ~ DM_mortality_NHWM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWF" ~ DM_mortality_NHWF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBM" ~ DM_mortality_NHBM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBF" ~ DM_mortality_NHBF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HM" ~ DM_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HF" ~ DM_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s]
    )
    
    p.death.CHD <- case_when(
      sim_out[,"DEMO",t] == "Male" ~ CHD_mortality_Male[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "Female" ~ CHD_mortality_Female[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWM" ~ CHD_mortality_NHWM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWF" ~ CHD_mortality_NHWF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBM" ~ CHD_mortality_NHBM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBF" ~ CHD_mortality_NHBF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HM" ~ CHD_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HF" ~ CHD_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s]
    )
    
    p.death.Stroke <- case_when(
      sim_out[,"DEMO",t] == "Male" ~ stroke_mortality_Male[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "Female" ~ stroke_mortality_Female[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWM" ~ stroke_mortality_NHWM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWF" ~ stroke_mortality_NHWF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBM" ~ stroke_mortality_NHBM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBF" ~ stroke_mortality_NHBF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HM" ~ stroke_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HF" ~ stroke_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s]
    )
    
    p.death <- case_when(
      sim_out[,"DEMO",t] == "Male" ~ Non_CVD_DM_mortality_Male[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "Female" ~ Non_CVD_DM_mortality_Female[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWM" ~ Non_CVD_DM_mortality_NHWM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHWF" ~ Non_CVD_DM_mortality_NHWF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBM" ~ Non_CVD_DM_mortality_NHBM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "NHBF" ~ Non_CVD_DM_mortality_NHBF[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HM" ~ Non_CVD_DM_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s],
      sim_out[,"DEMO",t] == "HF" ~ Non_CVD_DM_mortality_HM[as.numeric(sim_out[,"Age_cycle",t])-19,s]
    )
    
    
    #Disaggregating ASCVD
    p.CHD_first <- as.numeric(Prop_CHD[1,sim_out[,"DEMO",t]]) #Prob. of CHD event (Angina, MI, Fatal MI, Fatal CHD) Source: Benjamin et al, Circulation, 2018 Table 13-1 & 18-1, and 18-2
    p.CHD_recurrent <- as.numeric(Prop_CHD[2,sim_out[,"DEMO",t]])
    
    #Markov State #1: "No CVD, No Diabetes"
    p.H.2.DM <- as.numeric(sim_out[,"DM_prob",t])*as.numeric(sim_out[,"risk_adjustment.DM",t])
    p.H.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(sim_out[,"CVD_prob",t])
    p.H.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(sim_out[,"CVD_prob",t])
    p.H.2.H <- ifelse((p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death) > 1, 0,
                      1 - (p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death))
    
    #Markov State #2: "No CVD, With Diabetes
    p.DM.2.death <- (1-exp(-(p.death+p.death.DM)))
    p.DM.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(sim_out[,"CVD_prob",t])
    p.DM.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(sim_out[,"CVD_prob",t])
    p.DM.2.DM <- ifelse((p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death) > 1, 0,
                        1 - (p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death))
    
    #Markov State #3: "First Stroke"
    p.initial_Stroke.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke))))
    p.initial_Stroke.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - p.initial_Stroke.2.death)
    p.initial_Stroke.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - p.initial_Stroke.2.death, 0)
    
    #Markov State #4: "First CHD w/o RVSC" 
    p.initial_CHD_No_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD))))
    p.initial_CHD_No_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - p.initial_CHD_No_RVSC.2.death)
    p.initial_CHD_No_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - p.initial_CHD_No_RVSC.2.death, 0)
    
    #Markov State #5: "First CHD with RVSC"
    p.initial_CHD_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s]))))
    p.initial_CHD_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - p.initial_CHD_RVSC.2.death)
    p.initial_CHD_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - p.initial_CHD_RVSC.2.death, 0)
    
    #Markov State #6: "CVD History, No Diabetes"
    p.CVD.2.DM <- as.numeric(sim_out[,"DM_prob",t])*as.numeric(sim_out[,"DM_prob",t])
    p.CVD.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(sim_out[,"CVD_recurrent_prob",t])
    p.CVD.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(sim_out[,"CVD_recurrent_prob",t])
    p.CVD.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke)))
    p.CVD.2.CVD <- ifelse((p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death) > 1, 0,
                          1 - (p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death))
    
    #Markov State #7: "CVD History, With Diabetes"
    p.CVD_DM.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(sim_out[,"CVD_recurrent_prob",t])
    p.CVD_DM.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(sim_out[,"CVD_recurrent_prob",t])
    p.CVD_DM.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke+p.death.DM)))
    p.CVD_DM.2.CVD_DM <- ifelse((p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke) > 1, 0,
                                1 - (p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke))
    
    #Markov State #8: "Subsequent Stroke"
    p.sub_Stroke.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke))))
    p.sub_Stroke.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - p.sub_Stroke.2.death)
    p.sub_Stroke.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - p.sub_Stroke.2.death, 0)
    
    #Markov State #9: "Subsequent CHD w/o RVSC" 
    p.sub_CHD_No_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD))))
    p.sub_CHD_No_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - p.sub_CHD_No_RVSC.2.death)
    p.sub_CHD_No_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - p.sub_CHD_No_RVSC.2.death, 0)
    
    #Markov State #10: "Subsequent CHD with RVSC"
    p.sub_CHD_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s]))))
    p.sub_CHD_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - p.sub_CHD_RVSC.2.death)
    p.sub_CHD_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - p.sub_CHD_RVSC.2.death, 0)
    
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
    
    #Updating Disease history status
    sim_out[,"Diabetes",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("No CVD, With Diabetes", "CVD History, With Diabetes"), 1, 
                                     ifelse(sim_out[,"state",t+1] == "Death", NA, 
                                            ifelse(as.numeric(sim_out[,"Glucose",t+1]) > 126, 1, sim_out[,"Diabetes",t])))
    
    sim_out[,"CVD_history",t+1] <- ifelse(sim_out[,"state",t] %in% c("First Stroke", "First CHD w/o RVSC", "First CHD with RVSC", "CVD History, No Diabetes", "CVD History, With Diabetes",
                                                                     "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), 1, 
                                        ifelse(sim_out[,"state",t+1] == "Death", NA, sim_out[,"CVD_history",t]))
    
    sim_out[,"BMI",t+1] <- ifelse(sim_out[,"state",t+1] == "Death", NA, sim_out[,"BMI",t+1])
    sim_out[,"Obesity",t+1] <- ifelse(as.numeric(sim_out[,"BMI",t+1]) < 30, 0, 
                                      ifelse((as.numeric(sim_out[,"BMI",t+1]) >= 30), 1, NA))
    
    raw.input.data  <- cbind(as.numeric(sim_out[,"SEQN",t+1]), as.numeric(sim_out[,"Age_cycle",t+1]), as.numeric(sim_out[,"Female",t+1]), as.numeric(sim_out[,"Race",t+1]),
                             as.numeric(sim_out[,"SBP",t+1]),as.numeric(sim_out[,"DBP",t+1]),as.numeric(sim_out[,"HPT_Txt",t+1]), as.numeric(sim_out[,"BMI",t+1]), 
                             as.numeric(sim_out[,"Diabetes",t+1]),  as.numeric(sim_out[,"CVD_history",t+1]))
    colnames(raw.input.data) <- c("SEQN","Age","Female","Race","SBP", "DBP", "HPT_Txt", "BMI","Diabetes","CVD_history")
    raw.input.data <- as.data.frame(raw.input.data)
    
    # Update QALYs and costs
    sim_out[,"HRQOL_scores",t+1] <- calc_HRQOL(raw.input.data, HRQOL_parameter_sim[,s])
    sim_out[,"HRQOL_scores",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "Subsequent Stroke"), as.numeric(sim_out[,"HRQOL_scores",t+1])+u_stroke_sim[s], 
                                           ifelse(sim_out[,"state",t+1] %in% c("First CHD w/o RVSC", "First CHD with RVSC", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), as.numeric(sim_out[,"HRQOL_scores",t+1])+u_CHD_sim[s],
                                                  ifelse(sim_out[,"state",t+1] == "Death", 0, sim_out[,"HRQOL_scores",t+1])))
    
    sim_out[,"HCE_predict",t+1] <- calc_HCE(raw.input.data, HCE_parameter_sim[,s])
    sim_out[,"HCE_predict",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "Subsequent Stroke"), as.numeric(sim_out[,"HCE_predict",t+1])+c_stroke_sim[s], 
                                          ifelse(sim_out[,"state",t+1] %in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC"), as.numeric(sim_out[,"HCE_predict",t+1])+c_CHD_sim[s], 
                                                 ifelse(sim_out[,"state",t+1] %in% c("First CHD with RVSC", "Subsequent CHD with RVSC"), as.numeric(sim_out[,"HCE_predict",t+1])+c_CHD_sim[s]+c_RVSC_sim[s],
                                                        ifelse(sim_out[,"state",t+1] == "Death", NA, sim_out[,"HCE_predict",t+1]))))
    
    if (intervention == "Policy") {
      sim_out[,"Food_cost",t+1] <- ifelse(sim_out[,"state",t+1] != "Death", c_policy, 0)
    } else {
      sim_out[,"Food_cost",t+1] <- 0
    }
    
    # Discount effects and costs
    sim_out[,"effect_disc",t+1] <- as.numeric(sim_out[,"HRQOL_scores",t+1])/((1+beta_QALY)^(t-1))
    sim_out[,"HCE_disc",t+1] <- as.numeric(sim_out[,"HCE_predict",t+1])/((1+beta_cost)^(t-1))
    sim_out[,"Food_cost_disc",t+1] <- as.numeric(sim_out[,"Food_cost",t+1])/((1+beta_cost)^(t-1))
  }
  # Subset output to include variables of interest
  out_variables <- c("Obesity", "Diabetes", "CVD_history", "HRQOL_scores", "HCE_predict", "Food_cost", 
                     "state", "HCE_disc", "effect_disc", "Food_cost_disc")
  sim_out_subset <- sim_out[,out_variables,]
  return (sim_out_subset)
}

# Detect system type
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
    if (os == "windows")
      os <- "windows"
  } else { ## if we still don't know
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

# Check if we are on Windows or Mac using our function.
cluster_type <- if (get_os() == "windows") {"PSOCK"} else {"FORK"}

no_cores <- detectCores() - 1
s = 1 # initialize iteration variable before forking
cl<-makeCluster(no_cores, type=cluster_type) # Make the cluster
# clusterEvalQ(cl)
registerDoParallel(cl)

acomb <- function(...) abind(..., along = 4)
sim_out <- foreach(s=1:n.sim, .combine = 'acomb', .verbose = T) %dopar% {
  set.seed(seed + n.cycle*s)
  run_sim(s)
}

stopCluster(cl)

print("Time to run model:")
proc.time() - model_start

#########################################################################
# 6 Summarize and save output                                           #
#########################################################################
print("Summarizing and saving output")
summ_start = proc.time()

# We will need state-specific output, so we'll subset sim_out here.
m.state <- array(sim_out[,"state",,], dim=c(n.individual, n.cycle+1, n.sim),
                 dimnames = list(paste("ind", 1:n.individual, sep = " "), paste("cycle", 0:n.cycle, sep = " "), paste("simulation_", 1:n.sim, sep = " ")))

# Create survey design objects
options(survey.lonely.psu = "adjust")
design <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WT_TOTAL, nest=TRUE, data=data_for_analysis)

# Define relevant variables for output
vars <- c("Obesity", "Diabetes", "CVD_history", "HRQOL_scores", "effect_disc", "HCE_predict", "HCE_disc", 
          "Food_cost", "Food_cost_disc", "Life Years", "Incident Stroke", "Incident CHD")

# First, summarize model output over time to get "lifetime" estimates.
# This gives us person-years for obesity, diabetes, and CVD prevalence, and total QALY and HCE values.
# These arrays will be used for calculating individual-level and simulation-level output
lifetime = array(NA, dim = c(length(vars), n.individual, n.sim),
                 dimnames = list(vars, data_for_analysis$Subject_ID, paste("simulation", 1:n.sim, sep = " ")))
for (i_var in c("Obesity", "Diabetes", "CVD_history", "HRQOL_scores", "HCE_predict", "HCE_disc", "effect_disc",
                "Food_cost", "Food_cost_disc")) {
  # First, create a numeric array specific to the variable
  temp = array(as.numeric(sim_out[,i_var,,]), dim = c(n.individual,n.cycle,n.sim))
  lifetime[i_var,,] = apply(temp, 3, rowSums, na.rm = T)
}
# We need to handle states differently
lifetime["Life Years",,] = apply(m.state!="Death", 3, rowSums, na.rm = T)
lifetime["Incident Stroke",,] = apply(m.state=="First Stroke", 3, rowSums, na.rm = T)
lifetime["Incident CHD",,] = apply(m.state=="First CHD with RVSC" | m.state=="First CHD w/o RVSC", 3, rowSums, na.rm = T)

# We will generate 3 different versions of the final output:
  # 1. Individual-level data (summarized over simulations and time)
  # 2. Simulation-level data (summarized over individuals and time)
  # 3. Time-level data (summarized over individuals and simulations)

# Initialize arrays to hold final output
ind_level_means = array(NA, dim = c(n.individual, length(vars)), 
                        dimnames = list(data_for_analysis$Subject_ID, vars))
ind_level_variances = array(NA, dim = c(n.individual, length(vars)), 
                       dimnames = list(data_for_analysis$Subject_ID, vars))
sim_level_means = array(NA, dim = c(n.sim, length(vars)), 
                        dimnames = list(paste("simulation_", 1:n.sim, sep = " "), vars))
sim_level_variances = array(NA, dim = c(n.sim, length(vars)), 
                            dimnames = list(paste("simulation_", 1:n.sim, sep = " "), vars))
time_level_means = array(NA, dim = c(n.cycle+1, length(vars)),
                         dimnames = list(paste("cycle", 0:n.cycle, sep = " "), vars))

# Loop through each variable and summarize
for (i_var in vars) {
  # First, pull in the lifetime array for the specific variable
  lifetime_var = lifetime[i_var,,]
  # Individual-level data: calculate the mean and variance across simulations
  ind_level_means[,i_var] = rowMeans(lifetime_var)
  ind_level_variances[,i_var] = apply(lifetime_var, 1, function(x) (sd(x)/sqrt(n.sim))^2)
  # Simulation-level data: calculate the survey-weighted mean and variance across individuals
  sim_level_means[,i_var] = apply(lifetime_var, 2, svymean, design)
  sim_level_variances[,i_var] = apply(lifetime_var, 2, function(x) SE(svymean(x, design))^2)
  }



# Create lists of the mean and variance arrays as the final output
output_ind_level = list(means = ind_level_means, vars = ind_level_variances)
output_sim_level = list(means = sim_level_means, vars = sim_level_variances)

# Save final output
saveRDS(output_ind_level, file = paste("03_Output/ind_level", intervention, "SEED", seed, n.sim, Sys.Date(), ".rda", sep = "_"))
saveRDS(output_sim_level, file = paste("03_Output/sim_level", intervention, "SEED", seed, n.sim, Sys.Date(), ".rda", sep = "_"))


saveRDS(RR_CHD, file = paste("03_Output/RR_CHD", intervention, "SEED", seed, Sys.Date(), ".rda", sep = "_"))
saveRDS(RR_stk, file = paste("03_Output/RR_stk", intervention, "SEED", seed, Sys.Date(), ".rda", sep = "_"))

###To generate 3-D time level data (i,vars,t)

time_level_mean_3d = array(NA,dim = c(n.individual, length(vars),(n.cycle+1)),
                      dimnames = list(data_for_analysis$Subject_ID, vars, c(1:(n.cycle+1))))

time_level_variance_3d = array(NA,dim = c(n.individual, length(vars),(n.cycle+1)),
                           dimnames = list(data_for_analysis$Subject_ID, vars, c(1:(n.cycle+1))))

for (i_var in c("Obesity", "Diabetes", "CVD_history", "HRQOL_scores", "HCE_predict", "HCE_disc", "effect_disc",
                "Food_cost", "Food_cost_disc")) {
temp = array(as.numeric(sim_out[,i_var,,]), dim = c(n.individual,n.cycle+1,n.sim))
#replace all NA with 0
temp[is.na(temp)] <-0

#  calculate the mean across simulations

sim_means = rowSums(temp, dims = 2)/n.sim

# Then, calculate the survey-weighted mean across individuals
time_level_mean_3d[,i_var,] = sim_means
time_level_variance_3d [,i_var, ] = apply(temp, 2, function(x) (sd(x)/sqrt(n.sim))^2)
}

temp_life = rowSums(m.state!="Death", dims=2)/n.sim
###SE for a proportion is calculated as sqrt(p/(1-p)/n)
temp_life_se=sqrt(rowSums(m.state!="Death", dims=2)*rowSums(m.state=="Death", dims=2)/n.sim^3)
time_level_mean_3d[,"Life Years",] = temp_life
time_level_variance_3d[,"Life Years",] = temp_life_se

temp_stroke = rowSums(m.state=="First Stroke", dims=2)/n.sim
temp_stroke_se = sqrt(rowSums(m.state=="First Stroke", dims=2)*rowSums(m.state!="First Stroke", dims=2)/n.sim^3)

time_level_mean_3d[,"Incident Stroke",] = temp_stroke
time_level_variance_3d[,"Incident Stroke",] = temp_stroke_se


temp_chd = rowSums(m.state=="First CHD with RVSC" | m.state=="First CHD w/o RVSC", dims=2)/n.sim
temp_chd_se = rowSums(m.state=="First CHD with RVSC" | m.state=="First CHD w/o RVSC", dims=2)*(n.sim-rowSums(m.state=="First CHD with RVSC" | m.state=="First CHD w/o RVSC", dims=2))/n.sim^3


time_level_mean_3d[,"Incident CHD",] = temp_chd
time_level_variance_3d[,"Incident CHD",] = temp_chd_se

output_time_level_3d = list(means = time_level_mean_3d, vars =time_level_variance_3d)
saveRDS(output_time_level_3d, file = paste("03_Output/time_level_3d", intervention, "SEED", seed, Sys.Date(), ".rda", sep = "_"))





print("Time to summarize and save output:")
proc.time() - summ_start

print("Total time:")
proc.time() - ptm
