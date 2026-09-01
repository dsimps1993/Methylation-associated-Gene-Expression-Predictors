###############################################################################
# mortality_basic_cox.r
#
# PURPOSE
#   Tests whether methylation-based senescence scores (p14, p16, p21 and their
#   mean) predict 10-year all-cause mortality in Generation Scotland, and
#   benchmarks them against a panel of established epigenetic clocks.
#
# MODEL
#   Surv(t_censor, dead) ~ scale(<score>) + age + sex
#
#   Each score is fitted in a separate Cox proportional hazards model. Scores
#   are z-scaled, so hazard ratios are per standard deviation.
#
# INPUTS   (paths relative to project root)
#   data/additional_genscot_files/2024-04-19_deaths.csv
#       One row per participant recorded as deceased (linkage extract).
#   data/additional_genscot_files/covar_temp2.csv
#       QC'd covariates and senescence scores. Requires: id, DNAm_ID, age,
#       sex, t_censor (years from baseline appointment to death or censoring).
#   data/additional_genscot_files/GS_clocks.csv
#       Published epigenetic clock estimates, keyed on DNAm_ID.
#
# OUTPUT
#   results/score_results/sen_mortality_hr_results_10Aug2026.csv
#       One row per score: HR, 95% CI, Z, P, case/sample counts, and the
#       Schoenfeld residual tests of the proportional hazards assumption
#       (local = the score term itself; global = the model as a whole).
#
# DEPENDENCIES
#   R with the survival package.
###############################################################################

### read in death data and QCd clock data ###
mortality_data<- read.csv("data/additional_genscot_files/2024-04-19_deaths.csv")
d <- read.csv("data/additional_genscot_files/covar_temp2.csv")
additional_clocks <- read.csv("data/additional_genscot_files/GS_clocks.csv")

#merge clock data with covariate data
d <- merge(d, additional_clocks, by="DNAm_ID")

### create 10-year death phenotype and tte ###
temp <- which(d$id %in% mortality_data$id)
d$dead <- 0
d$dead[temp] <- 1
temp2 <- which(d$dead==1 & d$t_censor>10)  
d$dead[temp2] <- 0
d$t_censor[d$t_censor>10] <- 10

# List of clocks to loop through
clocks <- c("p14", "p16", "p21", "mean_senesc", "Horvathv1" , "Hannum", "Lin", 
          "PhenoAge", "YingCausAge", "YingDamAge", "YingAdaptAge", "Horvathv2", "Zhang_10", 
          "DunedinPoAm38", "DunedinPACE", "DNAmGrimAge", "DNAmGrimAge.1", "DNAmTL")


# objects to store results
res.cox <- setNames(vector("list", length(clocks)), clocks)
res.cox.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.zph.global <- setNames(vector("list", length(clocks)), clocks)

library(survival)

# Cox Regression model
for(clock in clocks){
      model_formula <- as.formula(paste0("Surv(t_censor, dead) ~ scale(", clock, ") + age + sex"))
      res.cox[[clock]] <- coxph(model_formula, data = d)
      res.cox.zph.local[[clock]] <- cox.zph(res.cox[[clock]])$table[1,] 
      res.cox.zph.global[[clock]] <- cox.zph(res.cox[[clock]])$table["GLOBAL",]
}

## prop hazards (local) ##
local_p = list()

for(clock in names(res.cox.zph.local)){
local_p[[clock]] <- res.cox.zph.local[[clock]]
}

local_p = as.data.frame(do.call('rbind', local_p))
local_p$Clock <- as.character(rownames(local_p))
names(local_p)[1:3] <- c("chisq_local", "df_local", "local_p")


## prop hazards (global) ##
global_p = list()

for(clock in names(res.cox.zph.global)){
global_p[[clock]] <- res.cox.zph.global[[clock]]
}

global_p = as.data.frame(do.call('rbind', global_p))
global_p$Clock <- as.character(rownames(global_p))
names(global_p)[1:3] <- c("chisq_global", "df_global", "global_p")


## Extract summary output from Cox models for plotting 
extract_coefficients <- function(models_list) {
  # Initialize an empty list to store results
  results <- list()
  mort_coefficients <- list()
   
    # Iterate over each epigenetic clock model
    for (clock in names(models_list)) {
      model <- models_list[[clock]]
      
        # Extract the first row of coefficients
        coeff <- summary(model)$coefficients[1,]
	      n_event <- summary(model)$nevent
	      n_total <- summary(model)$n
        mort_coefficients[[clock]] <- c(coeff, n_event, n_total)
    }
    
    # Combine the coefficients into a data frame for the current disease
    results <- do.call(rbind, mort_coefficients)  

  return(as.data.frame(results))
}

out <- extract_coefficients(res.cox)
names(out) <- c("logHR","HR", "SE","Z","P","N_cases","N_total")
out$Clock <- rownames(out)

out$LCI <- exp(out$logHR - 1.96*out$SE)
out$UCI <- exp(out$logHR + 1.96*out$SE)

out1 <- out[,c("Clock","N_cases","N_total","HR","LCI","UCI","Z","P")]

out2 <- merge(out1, local_p, by=c("Clock")) 
out3 <- merge(out2, global_p, by=c("Clock")) 


write.csv(out3, "results/score_results/sen_mortality_hr_results_10Aug2026.csv", row.names=F)
