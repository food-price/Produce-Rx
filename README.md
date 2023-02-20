# Evaluating the health and economic impact, and cost-effectiveness of the national Produce-Rx
This is an example of using the Tufts Diabetes, Obesity, and Cardioascular disease Model (DOCM) to evaluate the health and economic impact and cost-effectiveness of a policy intervention (national produce prescription). The DOCM model forecasts US trends in the prevalence of obesity, diabetes, and cardiovascular disease, as well as all-cause mortality, quality-adjusted life years, and health care costs. 
# Getting Started
# Required R Packages
The following packages are required to run the code. You can install packages in R by running install.packages('package_name').
* survey
* svMisc
* psych
* gdata
* dplyr
* plyr
* data.table
* foreach
* doParallel
* abind

# Model Settings and structure
The 01_DOCM_producerx.R script is the main script for the model. The Key model parts are listed below to give an overview for users to set up the parameters. 

* 0.Preparation 
* 0.1 Install/Read R Packages
* 0.2 Create Functions
* 0.3 Creating Working Directory
* 0.4 Source other key scripts

These scripts contains functions for diabetes risk prediction model, CVD risk prediction model, subsequent CVD risk prediction model, HrQOL estimation model, and healthcare cost prediciton model.

* 0.5 Model settings 

(Set n.sim for the number of simulation for parameter sampling, n.cycle for number of years of follow-up, and n.loop for number of replicate for each invidiual,these can be set manually or read from command line when submitted through cluster)
* 1 Defining and importing necessary imput parameters, creating n.sim random draws
(All inputs in .csv format are located in the folder 01_Input)
* 1.0 Study specific settings: produce RX project 
(Set Policy-effect size, costs, and discounting rate, and other parameters that is specific to your study ) 
* 1.1 Read in master input file with invidiual-level data and cleaning
(We used data from NHANES participants in 2013-2018 cycles aged 40-80 years with diabetes and food insecurity, n=757, representing 5.7 million US adults.)
* 1.2 Import diet-disease etiologic effects data inputs: Age-specific relative risk estimates between dietary intake and disease outcomes

Source: Micha et.al 2017 (https://jamanetwork.com/journals/jama/article-abstract/2608221); 
 Miller et.al 2022 (https://pubmed.ncbi.nlm.nih.gov/35113165/); 
* 1.3 Import health-state/Event specific mortality data 

Data source: CDC wonder
* 1.4 Secular trends in major risk factors: estimated as average annual percent change over 1999-2016 

Data source: NHANES 1999-2016. (Source STATA codes !!!) 
* 1.5 Gender- and race-specific proportiona of CHD cases Among ASCVD cases (CHD + Stroke) 

Source: 
* 1.6  Parameters for Health Care Expenditure Model

Parameters of a healthcare cost prediction model which estimate overall health care expenditures for each individual based on their unique characteristics, includinge age, sex, race/ethnicity, BMI, and clinical consitions including diabetes, hypertension, and hisory of CVD. The parameters represent the estimated average marginal effects to predict the incremental change in health care expenditures for each covariate. 
Source: The model was developed de nova based on 2014 to 2016 Medical Panel Expenditure Survey (MEPS) data (n=73,174). (Source STATA codes !!!)  
* 1.7 HrQOL Prediction Parameters. 

We used a previously developed HRQOL prediction model for the US nationally representative sample based on demographic, socioeconomic, and chronic disease factors. Health-related quality of life (HRQOL) was measured based on an established patient-based estimates of how changes in health status alter the quality of life, from a scale of 0.00 (death) to 1.00 (perfect health).
Source: Lubetkin EI, Jia H, Franks P, Gold MR. Relationship among sociodemographic factors, clinical conditions, and health-related quality of life: examining the EQ-5D in the U.S. general population. Qual Life Res. 2005;14(10):2187-96.)
* 1.8 additional health care cost parameters 

The weighted-average total payments for event and procedure-specific costs were estiamted based on the number of total discharges across relevant diagnosis-related groups from the Medicare Inpatient Prospective Payment data. (Source STATA codes !!!) 

* 1.9 Productivity costs associated with CHD, stroke

average annual productivity costs for stroke and CHD cases was calculated based on total annual productivity costs and population counts of stroke and CHD in the US in 2019 respectively. Source: https://www.heart.org/idc/groups/heart-public/@wcm/@adv/documents/downloadable/ucm_491513.pdf. $ were inflated from 2015 dollar to 2021 dollar

* 2 Estimating disease-specific risk, health-related quality of life (HrQOL), and healthcare expenditures (HCE) at the baseline #
* 3 Run the simulation model 
* 3.1 Run n.sim times of the simulation function in parallel processes, and then combine the results in sim_out, Run the function for each arm separately.

A model structure diagram for the simulation model is shown in seperate figure in this github folder.   
* 4 Summarize and save output                                          

Main outputs include the mean and uncertainties in policy and non-policy arm, and incremental changes on health and economic incidence CVD events (first CVD events and recurrent CVD events are reported seperately), life-years, quality-adjusted life-years, discounted healthcare costs, productivity costs, policy costs (policy specific), and the relavent uncertainties. We also generated results stratified by age, sex, race/ethnicity, education, family income, and insurance status in order to highlight the health disparities.

Contact: lulu_stat@hotmail.com
