#########################################################################
# 6 Summarize and save output                                           #
#########################################################################
setwd("/cluster/tufts/kimlab/lwang18/healthy_food_rx")
library(survey)
library(abind)
library(data.table)   # needed for fread/fwrite
library(svMisc)
library(psych)
library(gdata)
library(dplyr)
library(dampack)
library(writexl)
seed=513
scenario="sensitivity2"
n.cycle=5

print("Summarizing and saving output")
summ_start = proc.time()

NHANES<- fread("01_Input/NHANES/NHANES_1318_Imp100_dm_fi2.csv", stringsAsFactors = TRUE, data.table = FALSE)

# 1.2 Select only necessary variables from master input file
variables <- c("SEQN", "WTINT2YR", "WTMEC2YR", "WT_TOTAL", "SDMVPSU", "SDMVSTRA", 
               "Age", "Female", "Race", "CVD_history", "Diabetes", "edu", "pir",
               "Total_Chol","HDL", "SBP", "DBP", "HPT_Txt", "Smoking",
               "DM_family", "BMI", "height", "weight", "Glucose", "Trig","Private","Medicare","Medicaid","Dual")
NHANES <- NHANES[variables]
NHANES_age40_79 <- subset(NHANES, Age >= 40 & Age < 80)
NHANES_age40_79$Subject_ID <- c(1:nrow(NHANES_age40_79))
data_for_analysis <- NHANES_age40_79
data_for_analysis <- data_for_analysis[order(data_for_analysis$Subject_ID), ]
data_for_analysis$pir_cat[data_for_analysis$pir<1.3]<-1
data_for_analysis$pir_cat[data_for_analysis$pir>=1.3 &data_for_analysis$pir<3]<-2
data_for_analysis$pir_cat[data_for_analysis$pir>=3]<-3
data_for_analysis$edu_cat[data_for_analysis$edu %in% c(1,2)] <-1
data_for_analysis$edu_cat[data_for_analysis$edu==3]<-2
data_for_analysis$edu_cat[data_for_analysis$edu %in% c(4,5)]<-3


sim_out_Policy<-readRDS( paste("03_Output/sim_out_Policy", "SEED", 513, 1000, n.cycle, "2022-06-04", ".rda", sep = "_"))
sim_out_No_Policy<-readRDS( paste("03_Output/sim_out_No_Policy", "SEED", 513, 1000, n.cycle, "2022-06-04", ".rda", sep = "_"))


#load the processing functions 
source("02_Programs/4_processing_function.R")

# Create survey design object for full population
options(survey.lonely.psu = "adjust")
design_all<- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WT_TOTAL, nest=TRUE, data=data_for_analysis)

# Create survey design object
design_sample = subset(design_all, Age >= 40 & Age < 80)
# Pull out sample size and number of simulations from array subset
pop_count = as.numeric(dplyr::count(data_for_analysis, wt=WT_TOTAL)/100)

#calculate net costs from different perspective
#Total_cost_health<-sim_out_No_Policy[,"Food_cost_disc" , ]+ sim_out_No_Policy[,"Admin_costs" ,]+ sim_out_No_Policy[,"HCE_disc"  ,]
#Total_cost_societ<-sim_out_No_Policy[,"Food_cost_disc" , ]+ sim_out_No_Policy[,"Admin_costs" ,]+ sim_out_No_Policy[,"HCE_disc"  ,]+ sim_out_No_Policy[,"Prod_cost_disc" ,]
#No_Policy_timehoriz =abind(sim_out_No_Policy,    Total_cost_health,Total_cost_societ, along=2, make.names=TRUE)


#calculate net costs from different perspective
#Total_cost_health<-sim_out_Policy[,"Food_cost_disc" , ]+ sim_out_Policy[,"Admin_costs" ,]+ sim_out_Policy[,"HCE_disc"  ,]
#Total_cost_societ<-sim_out_Policy[,"Food_cost_disc" , ]+ sim_out_Policy[,"Admin_costs" ,]+ sim_out_Policy[,"HCE_disc"  ,]+ sim_out_Policy[,"Prod_cost_disc" ,]
#Policy_timehoriz =abind(sim_out_Policy,  Total_cost_health,Total_cost_societ, along=2, make.names=TRUE)


################################################################################################
########### 6.1: summary and output for overall population #######################################

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

write.csv(summary_table, paste("03_Output/01_processed/summary_table_", seed, scenario, n.cycle,  "yrs_",Sys.Date(), ".csv", sep = ""))
 
  pop_summary_table = summary_table
  pop_summary_table[, -c(1)]=summary_table[, -c(1)]*pop_count
  write.csv(pop_summary_table, paste("03_Output/01_processed/pop_summary_table_", seed, scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  
  varcea <- c("Total_cost_health", "Total_cost_societ", "effect_disc")

  diff_allsim_summary = calc_allsim_summary(diff_timehoriz ,  varcea, n.sim )*pop_count
 write.csv(diff_allsim_summary, paste("03_Output/01_processed/pop_summary_allsim_", seed, scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))


  # Calculate cost-effectiveness
  cea_table1 = calc_ce(summary_table,"Total_cost_health") 
cea_table1$perspective<-"healthcare" 
  cea_table2 = calc_ce(summary_table,"Total_cost_societ")
cea_table2$perspective<-"societ" 
 
cea_table=rbind(cea_table1, cea_table2)

  write.csv(cea_table, paste("03_Output/01_processed/cea_table_", seed, scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))

bysub=1
################################################################################################
########### 6.2: summary and output by population subgroups: age, sex, race/ethnicity #######################################

if (bysub==1) {
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

write.csv(summarybyins, paste("03_Output/01_processed/summarybyins", scenario,n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
write.csv(summarybyinspop, paste("03_Output/01_processed/summarybyinspop", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))


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
  
  ICERby<-rbind(ICERbysex, ICERbyrace)
  write.csv(summarybysex, paste("03_Output/01_processed/summarybysex",scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  write.csv(summarybysexpop, paste("03_Output/01_processed/summarybysexpop",scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
 

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


  ICERby<-rbind(ICERbyage,ICERbysex, ICERbyrace, ICERbyedu, ICERbypir)
  
  write.csv(summarybypir, paste("03_Output/01_processed/summarybypir", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  write.csv(summarybypirpop, paste("03_Output/01_processed/summarybypirpop", scenario, n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))
  

##output ICER by population subgroups   

  write.csv(ICERby, paste("03_Output/01_processed/ICERby", n.cycle, "yrs_", Sys.Date(), ".csv", sep = ""))  
}
