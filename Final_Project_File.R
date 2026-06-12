# 1. Header
# Natalie Miller
# December 2025
# This is a QMRA model quantifying the risks of waterborne exposure to gastrointestinal illness in the Chattahoochee River - comparing before and after precipitation events 
# Exposure scenario: individuals recreating in the Chattahoochee River for 1 hour the day before a precipitation event and the day after a precipitation event. 
# Pathogen of Interest: Norovirus (as I was able to find a ratio for E. coli (indicator, measured by USGS) to Norovirus, which has a dose-reponse model completed as a case study on QMRAwiki.
# link: https://qmrawiki.org/case-studies/norovirus-drinking-water

# Data acquisition: Data was taken over the course of a year - November 2024 to November 2025 - from a USGS monitoring station, USGS-02336000 (right at the top of Atlanta on the Chattahoochee River)
# Estimated levels of E. coli concentrations (FIB), uses turbidity and other values to estimate E. coli concentrations
# Estimated by regression equation, water, colonies per 100 milliliters [-1 Std Deviation Interval]
# link: https://waterdata.usgs.gov/monitoring-location/USGS-02336000/#period=P365D&dataTypeId=continuous-99407-1906668952&showMedian=true&showFieldMeasurements=true

# Quantitative Microbial Risk Assessment as support for bathing waters profiling (Federigi, et al) gives E. coli to pathogen ratios with a uniform distribution and a min max, from which I calculated a mean value to use as the ratio value to convert E coli to Norovirus. 
# link: https://www.sciencedirect.com/science/article/abs/pii/S0025326X20304367?via%3Dihub


### Note: I use a few terms interchangably in this code, I try to make it clear but just for full documentation, here are what they are: pre = dry = day before precipitation event .... after = post = wet = the day after the precipitation event 



# 2. Install and load packages and set seed
###Kept this section the same as what you did including the added libraries and setting the seed
install.packages("ggplot2") # for plotting
install.packages("dplyr")  # for data manipulation
install.packages("mc2d")  # for Monte Carlo simulations

library(ggplot2)
library(dplyr)
library(mc2d)
set.seed(1)



# 3. Define model parameters
# first let's do this for the exposure assessment
# swimming exposure: water ingestion
###I am keeping these values for my QMRA model as I am also doing recreational water quality 
water_ingested <- 32 # in mL, from Dufour et al. 2017, for a 1-hr swim
water_ingested_stdev <- 10 #  in mL



###inputting the dataframe with the values from dry days and wet days ... measured as cfu/100mL
###these values were taken from the USGS water monitoring site (linked above) where dates with precipitation >=0.1 were selected and then the E. coli concentrations were recorded for the day prior and the day after
dry_ecoli_real <- c(68,84,73,150,75,78,17,62,60,69,83,450,66,72,73,64,70,200,140,99,95,88,160,80,130,180,120,140,270,230,140,140,120,320,310,150,270,380,82,99,76,83,120,83,89)
wet_ecoli_real <- c(850,350,150,120,1500,700,82,70,140,75,12000,1500,310,220,320,120,1900,200,710,450,180,160,160,410,3000,200,140,680,1900,8800,140,2800,600,1500,810,210,310,560,440,130,3400,94,190,90,300)


### the conversion ratio of pathogen to E. coli taken from the literature (linked above)
conversion_ratio <- 0.002416


### using the conversion ratio to alter the E. coli numbers to reflect estimated Norovirus numbers for the inputted data
dry_ecoli <- (dry_ecoli_real*conversion_ratio)/100
wet_ecoli <- (wet_ecoli_real*conversion_ratio)/100

  
#### now the mean and standard deviation for the pre (dry days) and after (wet days) from the altered dataframe 
pre_mean <- mean(dry_ecoli)
pre_stdev <- sd(dry_ecoli)
after_mean <- mean(wet_ecoli)
after_stdev <- sd(wet_ecoli)


### outputting the dry (pre) and wet (post) values for visual inspection
cat("Dry day Norovirus mean:", pre_mean, "\n")
cat("Dry day Norovirus stdev:", pre_stdev, "\n\n")

cat("Wet day Norovirus mean:", after_mean, "\n")
cat("Wet day Norovirus stdev:", after_stdev, "\n\n")





#### My model parameters for dose-response
#### QMRA wiki uses the beta-poisson model with these parameters (citation linked above):
alpha_value <- 0.136615746
N50_value <- 24153.13022






# 4. Functions developed for model
# the ways we calculate dose or probabilities need to be coded to generate distributions, estimates

# first we need to estimate dose: dose = V x C
#### keeping this the same as yours with altered variables since I am using recreational water as well
est_dose <- function(water_ingested, pre_mean){
  dose <- water_ingested*pre_mean
  return(dose)
}



#### here I have the approximate beta poisson pulled from QMRA wiki:
beta_poisson_model <- function(dose, alpha = alpha_value, N50 = N50_value){
  if(dose==0){return(0)}
  prob_infection <- 1 - (1 + (dose / N50) * (2^(1/alpha) - 1))^(-alpha)
  return(prob_infection)
}
    


### value taken from CDC guidelines stating up to 30% of infections are asymptomatic: https://www.cdc.gov/mmwr/preview/mmwrhtml/rr6003a1.htm
prob_illness_given_infection <- 0.7 # 70% of infections convert to illness 




# finally an equation to describe overall risk
### I am using yours as a base for mine but I edited it to reflect the beta-poisson model as well as removing some of the parameters you used that are not needed here (I don't think)
  est_illness_risk <- function(dose, prob_illness_given_infection){
    prob_infection <- beta_poisson_model(dose)
    risk_illness <- prob_infection*prob_illness_given_infection
    return(risk_illness)
}

  
  
  
  

# 5. first let's do a point-estimate/deterministic of risk to ensure it looks okay before moving to the more complicated analysis
#### I kept this mostly the same as yours with a few edits to reflect my variable names as well as the fact that I am using a beta-poisson model rather than exponential 
#### this whole section is for the pre-precipitation values (dry day)
mean_dose <- est_dose(water_ingested, pre_mean)
mean_prob_inf <- beta_poisson_model(mean_dose)
mean_prob_ill <- est_illness_risk(mean_dose, prob_illness_given_infection)

# Display results
###Keeping this because I like the way it depicts all the results cleanly and you said we could use parts of your code
cat("Point estimate- first order analysis")
cat("Mean Dose:", round(mean_dose, 3), "particles")
cat("Probability of Infection:", round(mean_prob_inf, 4))
cat("Probability of Illness:", round(mean_prob_ill, 4))
cat("Risk as a percentage:", round(mean_prob_ill * 100, 2),"%")

# everything looks okay- we can move on to the probabilistic analysis
##### The values here are extremely low, but I think that is okay because the level of Norovirus to E. coli conversion ratio creates very small levels of Norovirus. Additionally, these are dry days and while recreating results in some GI illness exposure risk, the percent that the risk is due to Norovirus is very low. 
##### These values I feel are more relavent when comparing to wet days as opposed to on their own for risk assessment *more important to look at the difference rather than the individual values :)






# 6. Probabilistic analysis- adding in Monte Carlo to capture uncertainties in parameter values
#### again, this is for the pre-precipitation values 
n_iterations <- 10000 # 10,0000 simulations. I already set the seed at the top of the script

mc_doses <- numeric(n_iterations) # now I'm making vectors to store the values over 10,000 iterations
mc_infection_probs <- numeric(n_iterations)
mc_illness_risks <- numeric(n_iterations)

# now I can create a for loop to randomly sample from the distributions for each parameter
### This code is mostly the smae as yours with edits to variables as needed for my values 
for (i in 1:n_iterations) {
  # Randomly sample water concentration (normal distribution)
  sampled_conc <- rnorm(1, mean = pre_mean, sd = pre_stdev)
  
  # Ensure non-negative concentration
  sampled_conc <- max(0, sampled_conc)
  
  # Calculate dose for this iteration
  sampled_dose <- est_dose(water_ingested, sampled_conc)
  mc_doses[i] <- sampled_dose
  
  # Calculate infection probability using exponential model
  mc_infection_probs[i] <- beta_poisson_model(sampled_dose)
  
  # Calculate illness risk
  mc_illness_risks[i] <- mc_infection_probs[i] * prob_illness_given_infection
}



# Calculate and report summary statistics
### again, I kept this part of your code because it is a good way to see the values and see if it makes sense. 
cat("Probabilistic results- Monte Carlo with 10,000 iterations")
cat("Mean Risk of Illness:", round(mean(mc_illness_risks), 4))
cat("Median Risk of Illness:", round(median(mc_illness_risks), 4))
cat("95th Percentile Risk:", round(quantile(mc_illness_risks, 0.95), 4))
cat("Min Risk:", round(min(mc_illness_risks), 5))
cat("Max Risk:", round(max(mc_illness_risks), 4))
cat("Risk as percentage:")
cat("Mean:", round(mean(mc_illness_risks) * 100, 2),"%")
cat("95th percentile:", round(quantile(mc_illness_risks, 0.95) * 100, 2), "%")








# 6. Inspect results visually
##### using a lot of the same code you used but with "pre" results rather than "adv"
pre_results <- data.frame(
iteration = 1:n_iterations,
dose = mc_doses,
infection_prob = mc_infection_probs,
illness_risk = mc_illness_risks,
log_dose_pre <- log10(mc_doses)
)

# Plot 1: Distribution of illness risk
plot1 <- ggplot(pre_results, aes(x = illness_risk)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(illness_risk), color = "Mean"), 
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median(illness_risk), color = "Median"), 
             linetype = "dashed", size = 1) +
  labs(title = "Distribution of Norovirus Illness Risk from Swimming Exposure - Pre-Precipitation",
       x = "Risk of Illness (probability)",
       y = "Frequency",
       color = "Statistic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

plot1


# Plot 2: Cumulative distribution (showing what % of runs fall below each risk level)
plot2 <- ggplot(pre_results, aes(x = sort(illness_risk))) +
  geom_line(aes(y = seq_along(sort(illness_risk)) / n_iterations), 
            color = "steelblue", size = 1) +
  labs(title = "Cumulative Distribution of Risk - Pre-precipitation",
       x = "Risk of Illness (probability)",
       y = "Cumulative Probability",
       subtitle = "Shows what % of iterations fall below each risk level") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

plot2


# Plot 3: Dose vs Infection Probability (showing dose-response relationship)
plot3 <- ggplot(pre_results, aes(x = dose, y = infection_prob)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Dose-Response Relationship (Beta-Poisson)",
       x = "Ingested Dose (viral particles)",
       y = "Probability of Infection") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

plot3



# 6. Sensitivity analysis- what impact does slightly altering each parameter have on the overall estimate?
# first need to make a function to do the sensitivity analysis
####kept this the same as well but changed input variables/values to the ones relavent to my use-case
run_sensitivity_model <- function(param_value, param_name, n_iter = 1000) {
  risks <- numeric(n_iter)
  for (i in 1:n_iter) {
    dose <- est_dose(water_ingested, pre_mean)
    if (param_name == "water_volume") {
      dose <- est_dose(param_value, pre_mean)
    } else if (param_name == "virus_conc") {
      dose <- est_dose(water_ingested, param_value)
    } else if (param_name == "alpha") {
      dose <- est_dose(water_ingested, pre_mean)
    }
    
    # Calculate risk using beta-poisson model (rather than exponential)
    prob_inf <- beta_poisson_model(dose, alpha = if (param_name == "alpha") param_value else alpha_value, N50 = N50_value)
    risks[i] <- prob_inf * prob_illness_given_infection
  }
  
  return(mean(risks))
}

# to test if function works- run things like this:
# risks
# length(risks)

# Test parameter ranges
###Again, this is similar to what you did to test if the function works but I added my own variables for the sensativity_df part and the sensitivity plot. 
water_volumes <- seq(0.001, 0.2, by = 0.01)
virus_concs <- seq(10, 500, by = 20)
alpha_values <- seq(0.01, 0.1, by = 0.005)

# Run sensitivity for each parameter
sensitivity_water <- sapply(water_volumes, 
                            function(x) run_sensitivity_model(x, "water_volume"))
sensitivity_conc <- sapply(virus_concs, 
                           function(x) run_sensitivity_model(x, "virus_conc"))
sensitivity_alpha <- sapply(alpha_values, 
                        function(x) run_sensitivity_model(x, "alpha"))


##to visually inspect the results 
sensitivity_df <- data.frame(
  parameter = c(rep("Water Volume (L)", length(water_volumes)),
                rep("Virus Concentration (Norovirus units/L)", length(virus_concs)),
                rep("Alpha (dose-response)", length(alpha_values))),
  value = c(water_volumes, virus_concs, alpha_values),
  risk = c(sensitivity_water, sensitivity_conc, sensitivity_alpha)
)



sensi_pre_plot <- ggplot(sensitivity_df, aes(x = value, y = risk, color = parameter)) +
  geom_line(size = 1) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(title = "Sensitivity Analysis: Impact of Parameters on Illness Risk",
       x = "Parameter Value (simulated)",
       y = "Mean Risk of Illness",
       color = "Parameter") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

sensi_pre_plot



#### The final summary values, similar to what you did but with my own problem-specific output 
cat("\nFinal Risk Summary\n")
cat("Scenario: Norovirus exposure during recreational swimming\n")
cat("Water volume ingested (mean): ", water_ingested, " liters\n")
cat("Mean Norovirus concentration (dry days): ", pre_mean, " units per 100 mL-equivalent\n")
cat("Dose-Response Model: Beta-Poisson\n")
cat("KEY FINDINGS:\n")
cat("Average risk of illness per swim: ", 
    round(mean(mc_illness_risks) * 100, 3), "%\n")
cat("In the worst 5% of scenarios, risk exceeds: ", 
    round(quantile(mc_illness_risks, 0.95) * 100, 3), "%\n")

cat("\nINTERPRETATION:\n")
cat("This model suggests that a single recreational swimming session has a ",
    round(mean(mc_illness_risks) * 100, 2),
    "% chance of resulting in Norovirus-related illness, under Pre-precipitation (dry) conditions.\n",
    "The Beta-Poisson dose-response model assumes nonconstant survival and infection probabilities.\n",
    "Following this analysis will be the wet-day (post-precipitation analysis).\n",
    sep = "")






#################### after analysis #######################


#### This section of code is basically the same as the previous expect with variable and function names changed to reflect "after" or "post-precipitation"
#### Because the code is basically the same, I will leave the comments very minimal as the prior comments from the pre section apply here

# 4. Functions developed for model (same as before)
est_dose <- function(water_ingested, mean_value){
  dose <- water_ingested * mean_value
  return(dose)
}

beta_poisson_model <- function(dose, alpha = alpha_value, N50 = N50_value){
  if(dose == 0){ return(0) }
  prob_infection <- 1 - (1 + (dose / N50) * (2^(1/alpha) - 1))^(-alpha)
  return(prob_infection)
}

prob_illness_given_infection <- 0.7  


# 5. first let's do a point-estimate/deterministic of risk to ensure it looks okay before moving to the more complicated analysis
after_mean_dose <- est_dose(water_ingested, after_mean)
after_mean_prob_inf <- beta_poisson_model(after_mean_dose)
after_mean_prob_ill <- after_mean_prob_inf * prob_illness_given_infection

# Display results
cat("Point estimate—AFTER precipitation event\n")
cat("Mean Dose:", round(after_mean_dose, 3), "particles\n")
cat("Probability Infection:", round(after_mean_prob_inf, 4), "\n")
cat("Probability Illness:", round(after_mean_prob_ill, 4), "\n")
cat("Risk (%):", round(after_mean_prob_ill * 100, 2), "%\n")





# 6. Probabilistic analysis- adding in Monte Carlo to capture uncertainties in parameter values
n_iterations <- 10000

after_mc_doses <- numeric(n_iterations)
after_mc_infection_probs <- numeric(n_iterations)
after_mc_illness_risks <- numeric(n_iterations)

for(i in 1:n_iterations){
  
  sampled_conc_after <- rnorm(1, mean = after_mean, sd = after_stdev)
  sampled_conc_after <- max(0, sampled_conc_after)
  
  sampled_dose_after <- est_dose(water_ingested, sampled_conc_after)
  after_mc_doses[i] <- sampled_dose_after
  
  after_mc_infection_probs[i] <- beta_poisson_model(sampled_dose_after)
  
  after_mc_illness_risks[i] <- after_mc_infection_probs[i] * prob_illness_given_infection
}

cat("AFTER-event Monte Carlo results (10,000 iterations)\n")
cat("Mean Risk:", round(mean(after_mc_illness_risks), 4), "\n")
cat("Median Risk:", round(median(after_mc_illness_risks), 4), "\n")
cat("95th Percentile:", round(quantile(after_mc_illness_risks, 0.95), 4), "\n")
cat("Min:", round(min(after_mc_illness_risks), 5), "\n")
cat("Max:", round(max(after_mc_illness_risks), 4), "\n")
cat("Risk (%) — Mean:", round(mean(after_mc_illness_risks)*100,2), "%\n")
cat("Risk (%) — 95th:", round(quantile(after_mc_illness_risks,0.95)*100,2), "%\n")







# 6. Inspect results visually
after_results <- data.frame(
  iteration = 1:n_iterations,
  dose = after_mc_doses,
  infection_prob = after_mc_infection_probs,
  illness_risk = after_mc_illness_risks,
  log_dose_after = log10(after_mc_doses)
)



# Plot 1: Distribution of illness risk
after_plot1 <- ggplot(after_results, aes(x = illness_risk)) +
  geom_histogram(bins = 50, fill = "firebrick", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(illness_risk), color = "Mean"),
             linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median(illness_risk), color = "Median"),
             linetype = "dashed", size = 1) +
  labs(title = "Distribution of Illness Risk — AFTER Precipitation",
       x = "Risk of Illness",
       y = "Frequency") +
  theme_minimal()

after_plot1


# Plot 2: Cumulative distribution (showing what % of runs fall below each risk level)
after_plot2 <- ggplot(after_results, aes(x = sort(illness_risk))) +
  geom_line(aes(y = seq_along(sort(illness_risk))/n_iterations),
            color = "firebrick", size = 1) +
  labs(title = "Cumulative Distribution of Illness Risk — AFTER",
       x = "Risk of Illness",
       y = "Cumulative Probability") +
  theme_minimal()

after_plot2


# Plot 3: Dose vs Infection Probability (showing dose-response relationship)
after_plot3 <- ggplot(after_results, aes(x = dose, y = infection_prob)) +
  geom_point(alpha = 0.3, color = "firebrick") +
  geom_smooth(method = "loess", color = "black", se = FALSE) +
  labs(title = "Dose–Response Relationship — AFTER",
       x = "Dose (viral particles)",
       y = "Probability of Infection") +
  theme_minimal()

after_plot3








# 6. Sensitivity analysis- what impact does slightly altering each parameter have on the overall estimate?
run_sensitivity_model_after <- function(param_value, param_name, n_iter = 1000){
  risks <- numeric(n_iter)
  for(i in 1:n_iter){
    
    if(param_name == "water_volume"){
      dose <- est_dose(param_value, after_mean)
      
    } else if(param_name == "virus_conc"){
      dose <- est_dose(water_ingested, param_value)
      
    } else if(param_name == "alpha"){
      dose <- est_dose(water_ingested, after_mean)
    }
    
    # Calculate risk using beta-poisson model (rather than exponential)
    prob_inf <- beta_poisson_model(dose, alpha = if (param_name == "alpha") param_value else alpha_value, N50 = N50_value)
    risks[i] <- prob_inf * prob_illness_given_infection
  }
  return(mean(risks))
}



# Test parameter ranges
water_volumes <- seq(0.001, 0.2, by = 0.01)
virus_concs_after <- seq(10, 500, by = 20)
alpha_values_after <- seq(0.01, 0.1, by = 0.005)

sensitivity_water_after <- sapply(water_volumes,
                                  function(x) run_sensitivity_model_after(x,"water_volume"))
sensitivity_conc_after <- sapply(virus_concs_after,
                                 function(x) run_sensitivity_model_after(x,"virus_conc"))
sensitivity_alpha_after <- sapply(alpha_values_after,
                                  function(x) run_sensitivity_model_after(x,"alpha"))

# Run sensitivity for each parameter
sensitivity_df_after <- data.frame(
  parameter = c(rep("Water Volume (L)", length(water_volumes)),
                rep("Virus Concentration (units/L)", length(virus_concs_after)),
                rep("Alpha (dose-response)", length(alpha_values_after))),
  value = c(water_volumes, virus_concs_after, alpha_values_after),
  risk = c(sensitivity_water_after, sensitivity_conc_after, sensitivity_alpha_after)
)


##to visually inspect the results 
sensi_after_plot <- ggplot(sensitivity_df_after, aes(x = value, y = risk, color = parameter)) +
  geom_line(size = 1) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(title = "Sensitivity Analysis — AFTER Precipitation",
       x = "Parameter Value",
       y = "Mean Illness Risk") +
  theme_minimal()

sensi_after_plot






#### The final summary values, similar to what you did but with my own problem-specific output (for after)
cat("\nFinal Risk Summary — AFTER precipitation\n")
cat("Scenario: Norovirus exposure during recreational swimming\n")
cat("Water volume ingested:", water_ingested, "L\n")
cat("Mean Norovirus concentration AFTER event:", after_mean, "units per 100 mL\n")
cat("Dose-Response Model: Beta-Poisson\n")
cat("Average illness risk per swim:", round(mean(after_mc_illness_risks)*100,3), "%\n")
cat("Worst 5% exceed:", round(quantile(after_mc_illness_risks,0.95)*100,3), "%\n\n")
cat("Interpretation:\n")
cat("This model suggests that a single recreational swimming session has a ",
    round(mean(after_mc_illness_risks) * 100, 2),
    "% chance of resulting in Norovirus-related illness, under Post-precipitation (wet) conditions.\n",
    "The Beta-Poisson dose-response model assumes nonconstant survival and infection probabilities.\n",
    "Following this analysis will be the wet-day (post-precipitation analysis).\n",
    sep = "")
cat("Post-event conditions show elevated risk due to increased concentrations.\n")













################# Comparing the two results!! #################



#### From the results section earlier, make sure they are appropriately named for this figure for clarity. 
pre_results$scenario   <- "Pre-Precipitation"
after_results$scenario <- "Post-Precipitation"


#### pull out just the columns needed from the resutls section and put it into a new dataframe for use in plotting
pre_just_columns   <- pre_results[, c("illness_risk", "scenario")]
after_Just_columns <- after_results[, c("illness_risk", "scenario")]


#### Combine the datasets for use in the plotting section of code
combined_results <- rbind(pre_just_columns, after_Just_columns)




###I used AI to figure out how to have pre-precipitation on the left. When I originally was making the plot, it kept putting "post precipitation" first which I attribute to it being first alphabetically. Anyway, this is my citation for this part. 
combined_results$scenario <- factor(
  combined_results$scenario,
  levels = c("Pre-Precipitation", "Post-Precipitation")
)





#### now to visualize the results using a boxplot (similar to what I did for homework 5)
#### I am using the labels Pre and Post Precipitation because I feel it makes the most sense (even though I use a lot of different things interchangably throughout the code)
comparison_boxplot <- ggplot(
  combined_results,
  aes(x = scenario, y = illness_risk, fill = scenario)
) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  theme_minimal() +
  labs(
    title = "Comparison of Illness Risk: Pre vs Post Precipitation",
    x = "",
    y = "Illness Risk (Probability)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )



comparison_boxplot
