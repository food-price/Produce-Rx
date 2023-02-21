# Purpose: Process output from ODC-M model - Healthy food rx analysis
# Program input: 4-D array from model output (dimensions: individuals, variables, cycles, simulations)
# Program output: average per-person values (summary_table), population totals including differences
# between Policy and No Policy (pop_summary_table), and cost-effectiveness table with ICERs (cea_table)
# for a specified time horizon and/or population


##### Initialize functions to process model output #################################################


calc_summary <- function(lifetime, vars, n.sim, n.sample, design_sample) {
  # Input: output of calc_timehoriz function (dimensions: indviiduals, variables, simulations); survey design 
  # object for the sample
  # Purpose: 
  # 1. Summarize across individuals (mean and within-simulation variance)
  # 2. Summarize across simulations (survey-weighted mean and within-simulation variance)
  # 3. Apply Rubin's rule for a final calculation of variance (combining sampling and parameter uncertainty)
  # Output: Summary table with per-person mean and 95% CI for each variable
  
  # Summarize across individuals
  # Initialize arrays to hold simulation-level data
  
  sim_level_means = array(NA, dim = c(n.sim, length(vars)), dimnames = list(paste("simulation_", 1:n.sim, sep = " "), vars))
  sim_level_variances = array(NA, dim = c(n.sim, length(vars)), 
                              dimnames = list(paste("simulation_", 1:n.sim, sep = " "), vars))
  # Loop through each variable and summarize
  for (i_var in vars) {
    # First, pull in the lifetime array for the specific variable
    lifetime_var = lifetime[,i_var,]
    # Individual-level data: calculate the mean and variance across simulations
    sim_level_means[,i_var] = apply(lifetime_var, 2, svymean, design_sample)
    sim_level_variances[,i_var] = apply(lifetime_var, 2, function(x) SE(svymean(x, design_sample))^2)
  }
  
  # Summarize across individuals  
  ## Testing different method of summarizing 
  summary = data.frame(Outcome = colnames(sim_level_variances))
  ## grand mean for each variable 
  summary$mean = apply(sim_level_means, 2, mean)
  # the average NHANES sampling variance,  
  within.var = apply(sim_level_variances, 2, mean)
  ##the variance of the 1000 simulated population means. 
  
  ##according to Rubin's rule 
  ##(https://journals.sagepub.com/doi/suppl/10.1177/0272989X20916442, formula (2) on page 3 of Appendix 1) , 
  ##between.var is equal to the square of the standard deviation of the 1000 means. 
  between.var = apply(sim_level_means, 2, function(x)(sd(x)^2))
  
  final.var <- within.var +  (1 + (1/n.sample))*between.var
  final.se <- sqrt(final.var)
  
  summary$LL = apply(sim_level_means, 2, function(x)quantile(x, .025, na.rm=TRUE))
  summary$UL = apply(sim_level_means, 2, function(x)quantile(x, .975, na.rm=TRUE)) 
  summary$se =   final.se
  summary$LL2 = summary$mean - 1.96 * final.se
  # Calculate upper limit of 95% CI
  summary$UL2 = summary$mean + 1.96 * final.se  
  return (summary)
}

# Generate CEA table
calc_ce <- function(pop_summary_table, totalcosts) {
  # Input: summary table with an Outcome column containing totalcosts and "effect_disc" values,
  # and mean.No_Policy and mean.Policy columns
  temp = pop_summary_table %>%
    filter(Outcome %in% c(totalcosts, "effect_disc"))
temp$Outcome[temp$Outcome == totalcosts]<-"Total_cost_disc"

  rownames(temp) = temp$Outcome
  temp = temp %>% select(mean.No_Policy, mean.Policy)
  cost_qaly = data.frame(t(temp))
  cea_table = dampack::calculate_icers(cost_qaly$Total_cost_disc, cost_qaly$effect_disc, c("No Policy", "Policy"))
  return (cea_table)
}


calc_allsim_summary <- function(lifetime,  vars, n.sim,design_sample) {
  # Input: output of calc_timehoriz function (dimensions: indviiduals, variables, simulations); survey design 
  # object for the sample
  # Purpose: 
  # 1. Summarize across individuals (mean and within-simulation variance)
  # 2. Summarize across simulations (survey-weighted mean and within-simulation variance)
  # 3. Apply Rubin's rule for a final calculation of variance (combining sampling and parameter uncertainty)
  # Output: Summary table with per-person mean and 95% CI for each variable
  
  # Summarize across individuals
  # Initialize arrays to hold simulation-level data
  
  sim_level_means = array(NA, dim = c(n.sim, length(vars)), dimnames = list(paste("simulation_", 1:n.sim, sep = " "), vars))
  sim_level_variances = array(NA, dim = c(n.sim, length(vars)), 
                              dimnames = list(paste("simulation_", 1:n.sim, sep = " "), vars))
  # Loop through each variable and summarize
  for (i_var in vars) {
    # First, pull in the lifetime array for the specific variable
    lifetime_var = lifetime[,i_var,]
    # Individual-level data: calculate the mean and variance across simulations
    sim_level_means[,i_var] = apply(lifetime_var, 2, svymean, design_sample)
    sim_level_variances[,i_var] = apply(lifetime_var, 2, function(x) SE(svymean(x, design_sample))^2)
  }

  return ( sim_level_means)
}



calc_summary_by <- function( byvar, diff_timehoriz_sub){
summaryby<-NULL
summaryby_pop<-NULL
for (n in 1:length(byvars)){
  a=byvar[n]
  b<-unique(data_for_analysis[,a])
for (i in b){
  # Subset model output by race/ethnicity and time horizon specified in settings
  diff_timehoriz_sub= diff_timehoriz[data_for_analysis[,a]==i,,]
  design_subsample = subset(design_all, Age >= 40 & Age < 80 & get(a)==i)
  n.subsample <- dim(diff_timehoriz_sub)[1]
  # Calculate population count
  NHANES_sub <- subset(data_for_analysis, Age >= 40 & Age < 80 & get(a)==i) 
  pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
  diff_subsummary2 = calc_summary(diff_timehoriz_sub, vars, n.sim, n.subsample, design_subsample)
  diff_subsummary2pop<-diff_subsummary2
  diff_subsummary2pop[, -1]=diff_subsummary2[, -1]*pop_counts
  diff_subsummary2[,"cat"]<-i
  diff_subsummary2[,"catvar"]<-a
  diff_subsummary2pop[,"cat"]<-i
  diff_subsummary2pop[,"catvar"]<-a
  summaryby<-rbind(summaryby,  diff_subsummary2)
  summaryby_pop<-rbind(summaryby,  diff_subsummary2pop)
}

}
  return (list(summaryby, summaryby_pop))
}


calc_allsim_by <- function( byvar, vars, diff_timehoriz_sub){
  allsummaryby<-NULL
  for (n in 1:length(byvars)){
    a=byvar[n]
    b<-unique(data_for_analysis[,a])
    for (i in b){
      # Subset model output by race/ethnicity and time horizon specified in settings
      diff_timehoriz_sub= diff_timehoriz[data_for_analysis[,a]==i,,]
      design_subsample = subset(design_all, Age >= 40 & Age < 80 & get(a)==i)
      allsim_by<-calc_allsim_summary(diff_timehoriz_sub,  vars, n.sim,design_subsample ) 
      # Calculate population count
      NHANES_sub <- subset(data_for_analysis, Age >= 40 & Age < 80 & get(a)==i) 
      pop_counts = as.numeric(dplyr::count(NHANES_sub, wt=WT_TOTAL)/100)
      allsim_by<-as.data.frame( allsim_by)
      allsim_bypop<-allsim_by
      allsim_bypop[, -1]=allsim_by[, -1]*pop_counts
      allsim_by$cat<-i
      allsim_by[,"catvar"]<-a
      allsim_bypop$cat<-i
      allsim_bypop[,"catvar"]<-a
      allsummaryby<-rbind(allsummaryby,  allsim_bypop)
    }
  }
  return ( allsummaryby)
}


# Generate CEA table
calc_ce <- function(pop_summary_table, totalcosts) {
  # Input: summary table with an Outcome column containing totalcosts and "effect_disc" values,
  # and mean.No_Policy and mean.Policy columns
  temp = pop_summary_table %>%
    filter(Outcome %in% c(totalcosts, "effect_disc"))
  temp$Outcome[temp$Outcome == totalcosts]<-"Total_cost_disc"
  
  rownames(temp) = temp$Outcome
  temp = temp %>% select(mean.No_Policy, mean.Policy)
  cost_qaly = data.frame(t(temp))
  cea_table = dampack::calculate_icers(cost_qaly$Total_cost_disc, cost_qaly$effect_disc, c("No Policy", "Policy"))
  return (cea_table)
}



