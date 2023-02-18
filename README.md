# Produce-Rx
This is an example of using the Tufts Diabetes, Obesity, and Cardioascular disease Model (DOCM) to evaluate the health and economic impact and cost-effectiveness of a policy intervention (national produce prescription). The DOCM model forecasts US trends in the prevalence of obesity, diabetes, and cardiovascular disease, as well as all-cause mortality, quality-adjusted life years, and health care costs. Results are stratified by race/ethnicity in order to highlight health disparities.

Getting Started
Required R Packages
The following packages are required to run the code. You can install packages in R by running install.packages('package_name').
survey
svMisc
psych
gdata
dplyr
plyr
data.table
foreach
doParallel
abind

Model Settings
The 01_TDCM_health_food_rxnew.R script is the main script for the model. 
The key structure for the model is listed below 
#################
# 0.Preparation #
#################
# 0.1 Install/Read R Packages
# 0.2 Create Functions 
# 0.2-1 To estimate gamma/beta parameters based on Mean and SD
# 0.2-2 Converting multi-year risk (e.g., 10-year ASCVD risk) to annual probabilities
# 0.3 Creating Working Directory 
# 0.4 Source other key scripts
# 0.5 Model settings (set manually or read from command line when submitted through cluster)
####################################################################################
# 1 Defining and Importing Necessary Impute parameters, creating n.sim random draws
####################################################################################
##$$Policy specific settings: produce RX project 
# 1.0 Policy-effect size, costs, and discounting rate
# 1.1 Read in master input file with invidiual-level data and cleaning
The model was populated with NHANES participants in 2013-2018 cycles aged 40-80 years with diabetes and food insecurity, n=757, representing 5.7 million US adults
The incremental changes in health and economic incidence CVD events (first CVD events and recurrent CVD events are reported seperately), life-years, quality-adjusted life-years, healthcare costs, productivity costs, policy costs (policy specific)
# 1.2 Diet-disease etiologic effects data inputs : Age-specific relative risk estimates between dietary intake and disease outcomes
# 1.3 Health-state/Event specific mortality data
# 1.4 Secular trends in major risk factors: estimated as average annual percent change over 1999-2016 
# 1.5 Gender- and race-specific proportiona of CHD cases Among ASCVD cases (CHD + Stroke)
# 1.6 Health Care Expenditure Model
# 1.7 HrQOL Prediction
# 1.8 additional health care cost parameters
# 1.9 Productivity costs associated with CHD, stroke
#################################################################################################################################
# 2 Estimating disease-specific risk, health-related quality of life (HrQOL), and healthcare expenditures (HCE) at the baseline #
#################################################################################################################################
#########################################################################
# 3 Main Model - PLEASE DO NOT MODIFY UNLESS YOU KNOW WHAT YOU'RE DOING #
#########################################################################
# 3.1 source the scrit that defined the simulation function for running n.sim simulations 
# 3.2 run n.sim times of the simulation function in parallel processes, and then combine the results in sim_out, Run the function for each arm separately 
## The output is an 3-dimensional array (n.sample, variables, n.sim )
#########################################################################
# 4 Summarize and save output                                           #
#########################################################################
