#Tufts Diabetes-Obesity, CVD Microsimulation (TDOCM) Model 
#The FOOD-PRICE project (https://food-price.org/)
#Healthy Food Prescription Analysis
#Primary Author: David D. Kim, Tufts Medical Center
#Contact: david.kim@bsd.uchicago.edu
#Secondary Author: Brianna Lauren, Lu Wang
#Contact: lulu_stat@hotmail.com
#Main Purpose: 
#1) To develop a microsimulation model to link individual risk factors with cardiometabolic diseases, including diabetes and CVD
#2) To generate estimations of changes in disease risks, mortality, and  from dietary changes

#Model outcomes: 
#Incidence CVD events: first CVD events and recurrent CVD events, life-years, quality-adjusted life-years, healthcare costs, policy costs (policy specific)

#################
# 0.Preparation #
#################

# Start the clock
ptm <- proc.time()

#memory.size(max=T) # this is Windows-specific and will not work on the HPC running Linux

# 0.1 Install/Read R Packages

library(survey)
library(svMisc)
library(psych)
library(gdata)
library(dplyr)
library(data.table)   # needed for fread/fwrite
library(foreach)
library(doParallel)
library(abind)
#library(ggplot2)
library(dampack)
library(writexl)

# remove all data from memory
rm(list = ls())

# 0.2 Create Functions 

# 0.2-1 To estimate gamma/beta parameters based on Mean and SD
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

# 0.2-2 Converting multi-year risk (e.g., 10-year ASCVD risk) to annual probabilities
Multi_yr_Risk_to_annual_prob <- function(time, risk) {
  annual_rate <- -(1/time)*log(1-risk)
  annual_prob <- 1 - exp(-annual_rate)
}

# 0.3 Creating Working Directory 

#setwd("/cluster/tufts/kimlab/lwang18/healthy_food_rx")
setwd("C:/Users/lwang18/Documents/GitHub/healthy_food_rx")


# 0.4 Source other key scripts
source("02_Programs/1a - Diabetes_risk_prediction_FHS (categorical points).R")
source("02_Programs/1b - ASCVD_risk_calculator.R")
source("02_Programs/1c - FHS Subsequent_CVD_risk_calculator.R")
source("02_Programs/2a - HrQOL estimator for US general population.R")
source("02_Programs/2b - HCE estimator for US general population.R")
# source the script that defined the simulation function for running n.sim simulations 
source("02_Programs/3_sim_function_4.R")

# 0.5 Model settings (set manually or read from command line when submitted through cluster)
args <- commandArgs(trailingOnly = TRUE)  # get all arguments after script.R
# If no arguments are read from command line, set modeling choices manually. 
# Otherwise, read modeling choices from command line.
if (length(args) == 0) {
  seed <- 1234
  n.sim <- 2 #Number of probablistic samplings
  n.cycle <-5 #Number of Cycle Length: How long does the model run (i.e., analytic time horizon)
  n.sample <- 200 #Number of individuals to be selected from the full sample; if full sample, enter "ALL"
  scenario<-"Hba1c&BMI"   ##$$specific for produce RX project 
  n.loop=5
  #n.loop is the Number of replicates for each individual, each person is replicated for n.loop times and the outcomes are averaged across these replicates 
  #for the final results in each probablistic sampling. 
  bysub<-1 
  } else {
  # expecting 5 arguments
  if (length(args) != 5) {
    stop("ERROR: Incorrect number of command line arguments", call. = FALSE)
  }
  seed <- as.numeric(args[1]) # extracting first arg as seed and attempt to cast to number
  n.sim <- as.numeric(args[2])
  n.cycle <- as.numeric(args[3])
   # Assume entire sample
  n.sample = "ALL"
  scenario=args[4]
  n.loop=as.numeric(args[5])
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
}


set.seed(seed)

################################################################################################
# 1 Defining and Importing Necessary Impute parameters, creating n.sim random draws
################################################################################################

##1.0 Study specific settings: produce RX project 

# Policy-effect size, costs, and discounting rate

#get n.sim random draws on policy effect size based mean and SE of the effect, and assumed distribution
policy_effect_fruit_sim <- calc_nsims_rbeta(n.sim, mu = 0.40, se = 0.18) 
policy_effect_veg_sim <- calc_nsims_rbeta(n.sim, mu = 0.40, se = 0.18)
policy_effect_HA1c_sim <- rnorm(n.sim, -0.63,  0.18)
policy_effect_BMI_sim <- rnorm(n.sim, -0.45,  0.10)

if (intervention == "Policy") {
  policy_effect_fruit <-   policy_effect_fruit_sim
  policy_effect_veg <-  policy_effect_veg_sim
} else {
  policy_effect_fruit <- policy_effect_veg <- rep(0, n.sim)
}

# Policy costs: 31.86  per month (SE=4.4), scale to per year
c_policy = rnorm(n.sim, 382, 52.8) # weighted average of dollar amount in 13 studies included in the meta-analysis (annual, 2021 dollar)
# discounting rate
beta_cost <- beta_QALY <- 0.03 #Annual discounting rate of costs and QALYs

# 1.1 Read in master input file and cleaning

print('Importing data')

NHANES<- fread("01_Input/NHANES/NHANES_1318_Imp_dm_fi.csv", stringsAsFactors = TRUE, data.table = FALSE)

# 1.1.1 Select only necessary variables from master input file
variables <- c("SEQN", "WTINT2YR", "WTMEC2YR", "WT_TOTAL", "SDMVPSU", "SDMVSTRA", 
               "Age", "Female", "Race", "CVD_history", "Diabetes", "edu", "pir",
               "Total_Chol","HDL", "SBP", "DBP", "HPT_Txt", "Smoking",
               "DM_family", "BMI", "height", "weight", "Glucose", "Trig","Ins","Private","Medicare","Medicaid","Dual")

NHANES <- NHANES[variables]

# 1.1.2 Select the target starting population to model (Age 40-79) 
NHANES_age40_79 <- subset(NHANES, Age >= 40 & Age < 80)

NHANES_age40_79$Subject_ID <- c(1:nrow(NHANES_age40_79))

# 1.1.3 Define the analytic dataset to carry it forward

if (n.sample == "ALL"){
  data_for_analysis <- NHANES_age40_79
} else {
  random_sample <- sample(1:nrow(NHANES_age40_79),n.sample,replace=F)
  data_for_analysis <- NHANES_age40_79[random_sample,]
}

data_for_analysis <- data_for_analysis[order(data_for_analysis$Subject_ID), ]


Agesex<-ifelse(data_for_analysis$Age<65 & data_for_analysis$Female==1,1,
               ifelse(data_for_analysis$Age<65 & data_for_analysis$Female==0,2,
                      ifelse(data_for_analysis$Age>=65 & data_for_analysis$Female==1,3, 
                            ifelse(data_for_analysis$Age>=65 & data_for_analysis$Female==0,4,0))))
                      
# 1.1.4 Additional DATA CLEANING

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
data_for_analysis$pir_cat[data_for_analysis$pir<1.3]<-1
data_for_analysis$pir_cat[data_for_analysis$pir>=1.3 &data_for_analysis$pir<3]<-2
data_for_analysis$pir_cat[data_for_analysis$pir>=3]<-3
data_for_analysis$edu_cat[data_for_analysis$edu %in% c(1,2)] <-1
data_for_analysis$edu_cat[data_for_analysis$edu==3]<-2
data_for_analysis$edu_cat[data_for_analysis$edu %in% c(4,5)]<-3

# 1.1.5 DM risk adjustment for non-whites
data_for_analysis$risk_adjustment.DM <- ifelse(data_for_analysis$DEMO %in% c("NHWM", "NHWF", "Female", "Male"), 1.0,
                                               ifelse(data_for_analysis$DEMO %in% c("NHBM", "NHBM"), 1.5, 2.4))

# 1.1.6 Other inputs

##Initial health states.
name.health.state <- c("No CVD, No Diabetes", "No CVD, With Diabetes", "First Stroke", "First CHD w/o RVSC", "First CHD with RVSC",
                       "CVD History, No Diabetes", "CVD History, With Diabetes", "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC", "Death") 
n.health.state <- length(name.health.state) # number of health state that individuals can transition over time



# 1.2 Diet-disease etiologic effects data inputs : Age-specific relative risk estimates between dietary intake and disease outcomes
## Source: 
## this can be updated when updated evidence available 
RR_diet_disease <- fread("01_Input/final_diet_RRs_MainAnalysis v11.csv", stringsAsFactors = TRUE, data.table = FALSE)
RR_diet_disease <- subset(RR_diet_disease, minage>=35)

# 1.2.1 Direct effect on disease

RR_fruit_hstk <- subset(RR_diet_disease, outcome == 'HSTK' & riskfactor == 'fruit', select = c(minage, logRR.perchange, se.perchange))
RR_fruit_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'fruit', select = c(minage, logRR.perchange, se.perchange))
RR_fruit_istk <- subset(RR_diet_disease, outcome == 'ISTK' & riskfactor == 'fruit', select = c(minage, logRR.perchange, se.perchange))

RR_veg_hstk <- subset(RR_diet_disease, outcome == 'HSTK' & riskfactor == 'veg', select = c(minage, logRR.perchange, se.perchange))
RR_veg_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'veg', select = c(minage, logRR.perchange, se.perchange))
RR_veg_istk <- subset(RR_diet_disease, outcome == 'ISTK' & riskfactor == 'veg', select = c(minage, logRR.perchange, se.perchange))

RR_bmi_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'bmi', select = c(minage, logRR.perchange, se.perchange))
RR_bmi_tstk <- subset(RR_diet_disease, outcome == 'TSTK' & riskfactor == 'bmi', select = c(minage, logRR.perchange, se.perchange))

# Create n.sim random draws for logrr

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

# 1.2.2 Indirect effect on disease
# Weight change due to change in F&V intake (4-year weight change converted to yearly)
wt_change_fruit = rnorm(n.sim, -0.53, 0.046)/4
wt_change_veg = rnorm(n.sim, -0.25, 0.051)/4


# 1.3 Health-state/Event specific mortality data  
Non_CVD_DM_mortality <- fread("01_Input/Non_DM_IHD_Stroke_Cause_Mortality (Annual probability) (2001-2016 & 85+ corrected).csv", stringsAsFactors = TRUE, data.table = FALSE)
Non_CVD_DM_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Male, Non_CVD_DM_mortality$Male_SE))
Non_CVD_DM_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Female, Non_CVD_DM_mortality$Female_SE))
Non_CVD_DM_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWM, Non_CVD_DM_mortality$NHWM_SE))
Non_CVD_DM_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWF, Non_CVD_DM_mortality$NHWF_SE))
Non_CVD_DM_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBM, Non_CVD_DM_mortality$NHBM_SE))
Non_CVD_DM_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBF, Non_CVD_DM_mortality$NHBF_SE))
Non_CVD_DM_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HM, Non_CVD_DM_mortality$HM_SE))
Non_CVD_DM_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HF, Non_CVD_DM_mortality$HF_SE))

stroke_mortality <- fread("01_Input/Stroke_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
stroke_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Male, stroke_mortality$Male_SE))
stroke_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Female, stroke_mortality$Female_SE))
stroke_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWM, stroke_mortality$NHWM_SE))
stroke_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWF, stroke_mortality$NHWF_SE))
stroke_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBM, stroke_mortality$NHBM_SE))
stroke_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBF, stroke_mortality$NHBF_SE))
stroke_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HM, stroke_mortality$HM_SE))
stroke_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HF, stroke_mortality$HF_SE))

CHD_mortality <- fread("01_Input/IHD_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
CHD_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Male, CHD_mortality$Male_SE))
CHD_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Female, CHD_mortality$Female_SE))
CHD_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWM, CHD_mortality$NHWM_SE))
CHD_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWF, CHD_mortality$NHWF_SE))
CHD_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBM, CHD_mortality$NHBM_SE))
CHD_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBF, CHD_mortality$NHBF_SE))
CHD_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HM, CHD_mortality$HM_SE))
CHD_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HF, CHD_mortality$HF_SE))

DM_mortality <- fread("01_Input/DM_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
DM_mortality_Male <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Male, DM_mortality$Male_SE))
DM_mortality_Female <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Female, DM_mortality$Female_SE))
DM_mortality_NHWM <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWM, DM_mortality$NHWM_SE))
DM_mortality_NHWF <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWF, DM_mortality$NHWF_SE))
DM_mortality_NHBM <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBM, DM_mortality$NHBM_SE))
DM_mortality_NHBF <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBF, DM_mortality$NHBF_SE))
DM_mortality_HM <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HM, DM_mortality$HM_SE))
DM_mortality_HF <- t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HF, DM_mortality$HF_SE))

# 1.4 Secular trends in major risk factors: estimated as average annual percent change over 1999-2016 
risk_factor_trend <- fread("01_Input/Risk_Factors_Trends_Summary (Final).csv", stringsAsFactors = TRUE, data.table = FALSE)
risk_factor_trend_sim <- matrix(NA, nrow=nrow(risk_factor_trend), ncol=n.sim)
colnames(risk_factor_trend_sim) <- paste("simulation_", 1:n.sim, sep = " ")

risk_factor_trend_sim <- t(mapply(rnorm, n=n.sim, risk_factor_trend$APC, risk_factor_trend$APC_SE))

for (i in 1:nrow(risk_factor_trend)) {
  risk_factor_trend_sim[i,] <- as.matrix(rnorm(n.sim, risk_factor_trend$APC[i], risk_factor_trend$APC_SE[i]))
}

risk_factor_trend_sim <- cbind(risk_factor_trend, risk_factor_trend_sim)

# 1.5 Gender- and race-specific proportiona of CHD cases Among ASCVD cases (CHD + Stroke)
Prop_CHD <- fread("01_Input/Prop_CHD.csv", stringsAsFactors = FALSE, data.table = FALSE)

#1.6 Health Care Expenditure Model
HCE_parameter_sim <- matrix(NA, nrow=nrow(HCE_parameters), ncol=n.sim)
rownames(HCE_parameter_sim) <- rownames(HCE_parameters)
colnames(HCE_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

HCE_parameter_sim_diabe <- HCE_parameter_sim
HCE_parameter_sim_diabe[8,]<-HCE_parameter_sim[8,]*(1+policy_effect_HA1c_sim*0.13)

for (i in 1:nrow(HCE_parameters)) {
  HCE_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HCE_parameters$Beta[i], HCE_parameters$SE[i]))
}

#1.7 HrQOL Prediction
HRQOL_parameter_sim <- matrix(NA, nrow=nrow(HRQOL_parameters), ncol=n.sim)
rownames(HRQOL_parameter_sim) <- rownames(HRQOL_parameters)
colnames(HRQOL_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

for (i in 1:nrow(HRQOL_parameters)) {
  HRQOL_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HRQOL_parameters$Beta[i], HRQOL_parameters$SE[i]))
}

#1.8 additional health care cost parameters

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


mortality_CABG <- 0.0178 #In-hospital death rate, Source: Benjamin et al, Circulation, 2018, Table 24-1
mortality_CABG_sim <- calc_nsims_rbeta(n.sim, mortality_CABG, mortality_CABG*0.2)

mortality_PCI <- 0.0207 #In-hospital death rate, Source: Benjamin et al, Circulation, 2018, Table 24-1
mortality_PCI_sim <- calc_nsims_rbeta(n.sim, mortality_PCI, mortality_PCI*0.2)

p.death.RVSC <- prop_PCI*mortality_PCI + (1-prop_PCI)*mortality_CABG
p.death.RVSC_sim <- prop_PCI_sim*mortality_PCI_sim + (1-prop_PCI_sim)*mortality_CABG_sim


# 1.9 Productivity costs associated with CHD, stroke
#ppp Add productivity costs associated with stroke, and CHD respectively 
#average annual productivity costs for stroke and CHD cases was calculated based on total annual productivity costs and population counts of stroke and CHD in the US in 2019 respectively.
#Source: https://www.heart.org/idc/groups/heart-public/@wcm/@adv/documents/downloadable/ucm_491513.pdf
#$ were inflated from 2015 dollar to 2021 dollar

c_prod_stroke =4661
c_prod_stroke_sim = calc_nsims_rgamma(n.sim, c_prod_stroke, c_prod_stroke*0.2)

c_prod_CHD =6833
c_prod_CHD_sim = calc_nsims_rgamma(n.sim, c_prod_CHD, c_prod_CHD*0.2)

c_prod_CVD =6163
c_prod_CVD_sim = calc_nsims_rgamma(n.sim, c_prod_CVD, c_prod_CVD*0.2)

# 1.10 Other modeling Input


#################################################################################################################################
# 2 Estimating disease-specific risk, health-related quality of life (HrQOL), and healthcare expenditures (HCE) at the baseline #
#################################################################################################################################

# 2.1 FHS 8-year Diabetes Risk Prediction (With BMI)
variable_for_raw.input <- names(data_for_analysis)
raw.input.data <- data_for_analysis
data_for_analysis$DM_risk_8yr <- calc_DM_risk(raw.input.data)
data_for_analysis$DM_prob <- Multi_yr_Risk_to_annual_prob(time=8, risk=data_for_analysis$DM_risk_8yr)

# 2.2 ACC/AHA ASCVD 10-year Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
data_for_analysis$CVD_prob <- Multi_yr_Risk_to_annual_prob(time=10, risk=data_for_analysis$ASCVD_Risk_10yr)

# 2.3 FHS 2-year CVD recurrent Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
data_for_analysis$CVD_recurrent_prob <- Multi_yr_Risk_to_annual_prob(time=2, risk=data_for_analysis$CVD_Recurrent_risk_2yr)

# 2.4 Individual HrQOL prediction
raw.input.data <- data_for_analysis
data_for_analysis$HRQOL_scores <- calc_HRQOL(raw.input.data, HRQOL_parameters[,1])

# 2.5 Individual HCE prediction
raw.input.data <- data_for_analysis
data_for_analysis$HCE_predict <- calc_HCE(raw.input.data, HCE_parameters[,1])


n.individual <- nrow(data_for_analysis)*n.loop # number of individuals in the model

#########################################################################
# 3 Run the simulation model - PLEASE DO NOT MODIFY UNLESS YOU KNOW WHAT YOU'RE DOING #
#########################################################################

print('running model')
model_start = proc.time()

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

acomb <- function(...) abind(..., along = 3)

# 3.1 run n.sim times of the simulation function in parallel processes, and then combine the results in sim_out, Run the function for each arm separately 
## The output is an 3-dimensional array (n.sample, variables, n.sim )

sim_out_Policy <- foreach(s=1:n.sim, .combine = 'acomb', .verbose = T) %do% {
  set.seed(seed + n.cycle*s)
  run_sim(s, "Policy")
}
##save the output for policy arm
saveRDS(sim_out_Policy, file = paste("03_Output/sim_out_policy",  "SEED", seed, n.sim, n.cycle, n.loop, Sys.Date(), ".rda", sep = "_"))

sim_out_No_Policy <- foreach(s=1:n.sim, .combine = 'acomb', .verbose = T) %do% {
  set.seed(seed + n.cycle*s)
  run_sim(s, "No Policy")
}
##save the output for nonpolicy arm
saveRDS(sim_out_No_Policy, file = paste("03_Output/sim_out_nopolicy",  "SEED", seed, n.sim, n.cycle, n.loop, Sys.Date(), ".rda", sep = "_"))

stopCluster(cl)

print("Time to run model:")
proc.time() - model_start


#########################################################################
# 4 Summarize and save output                                           #
#########################################################################
print("Summarizing and saving output")
summ_start = proc.time()

#load the processing functions 
source("02_Programs/4_processing_function.R")

#process the data for subgroup analysis later
data_for_analysis$pir_cat[data_for_analysis$pir<1.3]<-1
data_for_analysis$pir_cat[data_for_analysis$pir>=1.3 &data_for_analysis$pir<3]<-2
data_for_analysis$pir_cat[data_for_analysis$pir>=3]<-3
data_for_analysis$edu_cat[data_for_analysis$edu %in% c(1,2)] <-1
data_for_analysis$edu_cat[data_for_analysis$edu==3]<-2
data_for_analysis$edu_cat[data_for_analysis$edu %in% c(4,5)]<-3

# Create survey design object for full population
options(survey.lonely.psu = "adjust")
design_all<- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WT_TOTAL, nest=TRUE, data=data_for_analysis)

# Create survey design object
design_sample = subset(design_all, Age >= 40 & Age < 80 )
# Pull out sample size and number of simulations from array subset
pop_count = as.numeric(dplyr::count(data_for_analysis, wt=WT_TOTAL))

################################################################################################
########### 4.1: summary and output for overall population #######################################

# Calculate difference between Policy and No Policy scenarios
diff_timehoriz = sim_out_Policy - sim_out_No_Policy
vars <- dimnames(diff_timehoriz)[[2]]
n.sample <- dim(diff_timehoriz)[1]
n.sim <- dim(diff_timehoriz)[3]

No_Policy_summary= calc_summary(sim_out_No_Policy, vars, n.sim, n.sample, design_sample)
Policy_summary = calc_summary(sim_out_Policy , vars, n.sim, n.sample, design_sample)
diff_summary1 = calc_summary(diff_timehoriz, vars, n.sim, n.sample, design_sample)

summary_table = full_join(No_Policy_summary, Policy_summary, by = "Outcome", suffix = c(".No_Policy", "")) %>%
  full_join(diff_summary1, by = "Outcome", suffix = c(".Policy", ".Diff1"))

write.csv(summary_table, paste("03_Output/01_processed/summary_table_", seed, scenario, n.cycle,  "yrs_",  Sys.Date(),  ".csv", sep = ""))

pop_summary_table = summary_table
pop_summary_table[, -c(1)]=summary_table[, -c(1)]*pop_count
write.csv(pop_summary_table, paste("03_Output/01_processed/pop_summary_table_", seed, scenario, n.cycle, "yrs_",  Sys.Date(), ".csv", sep = ""))

varcea <- c("Total_cost_health", "Total_cost_societ", "effect_disc")

diff_allsim_summary = calc_allsim_summary(diff_timehoriz ,  varcea, n.sim )*pop_count
write.csv(diff_allsim_summary, paste("03_Output/01_processed/pop_summary_allsim_", seed, scenario, n.cycle, "yrs_",Sys.Date(), ".csv", sep = ""))


# Calculate cost-effectiveness
cea_table1 = calc_ce(summary_table,"Total_cost_health") 
cea_table1$perspective<-"healthcare" 
cea_table2 = calc_ce(summary_table,"Total_cost_societ")
cea_table2$perspective<-"societ" 
cea_table=rbind(cea_table1, cea_table2)

write.csv(cea_table, paste("03_Output/01_processed/cea_table_", seed, scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))

################################################################################################
########### 4.2: summary and output by population subgroups: age, sex, race/ethnicity #######################################

###summary by age 

summarybyagepop<-NULL
summarybyage<-NULL
for (age in c(0:1)){
  # Subset model output by race/ethnicity and time horizon specified in settings
  diff_timehoriz_sub= diff_timehoriz[floor(data_for_analysis$Age/65)==age,,]
  design_subsample = subset(design_all, Age >= 40 & Age < 80 & floor(Age/65)==age)
  n.subsample <- dim(diff_timehoriz_sub)[1]
  # Calculate population count
  NHANES_sub <- subset(data_for_analysis, Age >= 40 & Age < 80 & floor(Age/65)==age) 
  pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
  
  diff_subsummary2 = calc_summary(diff_timehoriz_sub, vars, n.sim, n.subsample, design_subsample)
  diff_subsummary2pop<-diff_subsummary2
  diff_subsummary2pop[, -1]=diff_subsummary2[, -1]*pop_counts
  diff_subsummary2$Age<-age
  diff_subsummary2pop$Age<-age
  summarybyage<-rbind(summarybyage,  diff_subsummary2)
  summarybyagepop<-rbind(summarybyagepop,  diff_subsummary2pop)
}

summarybyage$ncycles=n.cycle
summarybyagepop$ncycles=n.cycle

ICERbyage<-(summarybyage[summarybyage$Outcome=="Total_cost_health", c("mean","LL","UL")])/ summarybyage[summarybyage$Outcome=="effect_disc", c("mean","LL","UL")]
ICERbyage$group<-c(1,2)

write.csv(summarybyage, paste("03_Output/01_processed/summarybyage", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
write.csv(summarybyagepop, paste("03_Output/01_processed/summarybyagepop", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))

#summarybysex

summarybysex<-NULL
summarybysexpop<-NULL
for (sex in c(0:1)){
  # Subset model output by race/ethnicity and time horizon specified in settings
  diff_timehoriz_sub= diff_timehoriz[data_for_analysis$Female==sex, ,]
  design_subsample = subset(design_all, Age >= 40 & Age < 80 & Female==sex)
  n.subsample <- dim(diff_timehoriz_sub)[1]
  # Calculate population count
  NHANES_sub <- subset(data_for_analysis, Age >= 40 & Age < 80 & Female==sex) 
  pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
  
  diff_subsummary2 = calc_summary(diff_timehoriz_sub, vars, n.sim, n.subsample, design_subsample)
  diff_subsummary2pop<-diff_subsummary2
  diff_subsummary2pop[, -1]=diff_subsummary2[, -1]*pop_counts
  diff_subsummary2$Female<-sex
  diff_subsummary2pop$Female<-sex
  
  summarybysex<-rbind(summarybysex,  diff_subsummary2)
  summarybysexpop<-rbind(summarybysexpop,  diff_subsummary2pop)
}


summarybysex$ncycles=n.cycle
summarybysexpop$ncycles=n.cycle

ICERbysex<-(summarybysex[summarybysex$Outcome=="Total_cost_health", c("mean","LL","UL")])/ summarybysex[summarybysex$Outcome=="effect_disc", c("mean","LL","UL")]
ICERbysex$group<-c(0,1)


write.csv(summarybysex, paste("03_Output/01_processed/summarybysex",scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
write.csv(summarybysexpop, paste("03_Output/01_processed/summarybysexpop",scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))

###summary by race
  summarybyrace<-NULL
  summarybyracepop<-NULL
  racecat<-sort(unique(data_for_analysis$Race))
  for (race in racecat){
    # Subset model output by race/ethnicity and time horizon specified in settings
    diff_timehoriz_sub= diff_timehoriz[data_for_analysis$Race==race,,]
    design_subsample = subset(design_all, Age >= 40 & Age < 80 & Race == race)
    n.subsample <- dim(diff_timehoriz_sub)[1]
    # Calculate population count
    NHANES_sub <- subset(data_for_analysis, Age >= 40 & Age < 80 & Race == race) 
    pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
    
    diff_subsummary2 = calc_summary(diff_timehoriz_sub, vars, n.sim, n.subsample, design_subsample)
    diff_subsummary2pop<-diff_subsummary2
    diff_subsummary2pop[, -1]=diff_subsummary2[, -1]*pop_counts
    diff_subsummary2$Race<-race
    diff_subsummary2pop$Race<-race
    summarybyrace<-rbind(summarybyrace,  diff_subsummary2)
    summarybyracepop<-rbind(summarybyracepop,  diff_subsummary2pop)
  }
  summarybyrace$ncycles=n.cycle
  summarybyracepop$ncycles=n.cycle
  
  ICERbyrace<-(summarybyrace[summarybyrace$Outcome=="Total_cost_health", c("mean","LL","UL")])/ summarybyrace[summarybyrace$Outcome=="effect_disc", c("mean","LL","UL")]
  ICERbyrace$group<-racecat
  
  write.csv(summarybyrace, paste("03_Output/01_processed/summarybyrace", scenario,n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  write.csv(summarybyracepop, paste("03_Output/01_processed/summarybyracepop", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  
 
  ###summary by education 
  summarybyedupop<-NULL
  summarybyedu<-NULL
  for (ed in c(1:3)){
    # Subset model output by race/ethnicity and time horizon specified in settings
    diff_timehoriz_sub= diff_timehoriz[data_for_analysis$edu_cat==ed,,]
    design_subsample = subset(design_all, Age >= 40 & Age < 80 & edu_cat ==ed)
    n.subsample <- dim(diff_timehoriz_sub)[1]
    # Calculate population count
    NHANES_sub <-subset(data_for_analysis, Age >= 40 & Age < 80 & edu_cat ==ed)
    pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
    
    diff_subsummary2 = calc_summary(diff_timehoriz_sub, vars, n.sim, n.subsample, design_subsample)
    
    diff_subsummary2pop<-diff_subsummary2
    diff_subsummary2pop[, -1]=diff_subsummary2[, -1]*pop_counts
    
    diff_subsummary2$edu_c<-ed
    diff_subsummary2pop$edu_c<-ed
    summarybyedu<-rbind(summarybyedu,  diff_subsummary2)
    summarybyedupop<-rbind(summarybyedupop,  diff_subsummary2pop)
  }
  
  summarybyedu$ncycles=n.cycle
  summarybyedupop$ncycles=n.cycle
  
  ICERbyedu<-(summarybyedu[summarybyedu$Outcome=="Total_cost_health", c("mean","LL","UL")])/ summarybyedu[summarybyedu$Outcome=="effect_disc", c("mean","LL","UL")]
  ICERbyedu$group<-c(1,2,3)
  write.csv(summarybyedu, paste("03_Output/01_processed/summarybyedu", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  write.csv(summarybyedupop, paste("03_Output/01_processed/summarybyedupop", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  
  ###summary by income 
  summarybypirpop<-NULL
  summarybypir<-NULL
  for (incm in c(1:3)){
    # Subset model output by race/ethnicity and time horizon specified in settings
    diff_timehoriz_sub= diff_timehoriz[data_for_analysis$pir_cat==incm,,]
    design_subsample = subset(design_all, Age >= 40 & Age < 80 & pir_cat==incm)
    n.subsample <- dim(diff_timehoriz_sub)[1]
    # Calculate population count
    NHANES_sub <- subset(data_for_analysis, Age >= 40 & Age < 80 & pir_cat==incm) 
    pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
    
    diff_subsummary2 = calc_summary(diff_timehoriz_sub, vars, n.sim, n.subsample, design_subsample)
    
    diff_subsummary2pop<-diff_subsummary2
    diff_subsummary2pop[, -1]=diff_subsummary2[, -1]*pop_counts
    
    diff_subsummary2$pirc<-incm
    diff_subsummary2pop$pirc<-incm
    summarybypir<-rbind(summarybypir,  diff_subsummary2)
    summarybypirpop<-rbind(summarybypirpop,  diff_subsummary2pop)
  }
  
  summarybypir$ncycles=n.cycle
  summarybypirpop$ncycles=n.cycle
  
  ICERbypir<-(summarybypir[summarybypir$Outcome=="Total_cost_health", c("mean","LL","UL")])/ summarybypir[summarybypir$Outcome=="effect_disc", c("mean","LL","UL")]
  ICERbypir$group<-c(1,2,3)
  
  
  write.csv(summarybypir, paste("03_Output/01_processed/summarybypir", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  write.csv(summarybypirpop, paste("03_Output/01_processed/summarybypirpop", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  
  ###summary by ins 

  summarybyins<-NULL
  summarybyinspop<-NULL
  Inscat<-c("Private","Medicare","Medicaid","Dual")
  for (ins in Inscat){
    # Subset model output by race/ethnicity and time horizon specified in settings
    diff_timehoriz_sub= diff_timehoriz[data_for_analysis[,ins]==1,,]
    design_subsample = subset(design_all, Age >= 40 & Age < 80 & get(ins) == 1)
    n.subsample <- dim(diff_timehoriz_sub)[1]
    # Calculate population count
    NHANES_sub <- subset(data_for_analysis, Age >= 40 & Age < 80 & get(ins) == 1) 
    pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
    
    diff_subsummary2 = calc_summary(diff_timehoriz_sub, vars, n.sim, n.subsample, design_subsample)
    diff_subsummary2pop<-diff_subsummary2
    diff_subsummary2pop[, -1]=diff_subsummary2[, -1]*pop_counts
    diff_subsummary2$Ins<-ins
    diff_subsummary2pop$Ins<-ins
    summarybyins<-rbind(summarybyins,  diff_subsummary2)
    summarybyinspop<-rbind(summarybyinspop,  diff_subsummary2pop)
  }
  summarybyins$ncycles=n.cycle
  summarybyinspop$ncycles=n.cycle
  
  ICERbyins<-(summarybyins[summarybyins$Outcome=="Total_cost_health", c("mean","LL","UL")])/ summarybyins[summarybyins$Outcome=="effect_disc", c("mean","LL","UL")]
  ICERbyins$group<-Inscat
  
  ##output ICER by population subgroups   
  
  ICERby<-rbind(ICERbyage,ICERbysex, ICERbyrace, ICERbyedu, ICERbypir, ICERbyins)
  
  write.csv(ICERby, paste("03_Output/01_processed/ICERby", n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))  
  
  
  write.csv(summarybyins, paste("03_Output/01_processed/summarybyins", scenario,n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  write.csv(summarybyinspop, paste("03_Output/01_processed/summarybyinspop", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  
print("Time to summarize and save output:")
proc.time() - summ_start

print("Total time:")
proc.time() - ptm
