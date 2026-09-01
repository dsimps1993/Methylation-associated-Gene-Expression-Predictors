###############################################################################
# alc_pulm_disease_only_fully_adjusted_cox.r
#
# PURPOSE
#   Tests whether the p16 and p21 senescence scores predict incident
#   idiopathic pulmonary fibrosis and alcoholic liver disease independently
#   of the established clinical risk factors for those conditions.
#
#   This is the companion to alc_pulm_disease_only_basic_cox.r. That script
#   fits each risk factor separately and compares effect sizes; this one puts
#   the risk factors into the model as covariates alongside the senescence
#   score.
#
# MODELS   (disease-specific adjustment sets)
#   pulm_fibrosis   Surv(t_event, event) ~ <score> + age + sex + height
#                                          + FEV + bmi + pack_years
#   liver_alc       Surv(t_event, event) ~ <score> + age + sex + bmi
#                                          + alc_units + pack_years
#
#   Each disease is adjusted for its own clinically relevant risk factors:
#   lung function and smoking for fibrosis, alcohol intake for liver disease.
#   Height accompanies FEV because lung volume scales with body size. All
#   predictors are z-scaled, so hazard ratios are per standard deviation.
#
# COHORT CONSTRUCTION
#   Identical to all_disease_incidence_basic_cox.r -- incident cases only,
#   one row per participant, administratively censored at 10 years, minimum
#   30 events -- but restricted to the two diseases named above.
#
# INPUTS   (project-root-relative unless noted)
#   data/additional_genscot_files/covar_temp2.csv
#       QC'd covariates and senescence scores.
#   data/additional_genscot_files/GS_additional_variables.csv
#       Additional phenotypes, including FEV and height.
#   data/additional_genscot_files/2026-02-17_disease_codes_combined.csv
#       Long-format linkage extract of diagnoses.
#   /exports/igmm/eddie/tchandra-lab/scrofts/GenScot_raw/covariates.csv
#       BMI and alcohol units. 
#
# OUTPUTS
#   results/score_results/alc_liver_pulm_disease_adjusted_cox_results.csv
#       Tidy per-score effect estimates with CIs and PH diagnostics.
#   results/score_results/alc_liver_pulm_disease_adjusted_cox_full_results.txt
#       Full printed coxph summaries, including the covariate coefficients
#       that the tidy export discards.
#
# DEPENDENCIES
#   R with dplyr and survival.
###############################################################################
 
# Load Required Packages
packages <- c("dplyr", "survival")
lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
  }
  library(x, character.only = TRUE)
})

# Read in clock data
covs <- read.csv("data/additional_genscot_files/covar_temp2.csv")
covs$sex <- ifelse(covs$sex == "F", 1, 0)

additional_covs <- read.csv("data/additional_genscot_files/GS_additional_variables.csv")

#merge with covs based on "DNAm_ID"
covs <- merge(covs, additional_covs, by="id", all.x=TRUE)

# Read in additional covariates (/exports/igmm/eddie/tchandra-lab/scrofts/GenScot_raw/covariates.csv)
add_covs <- read.csv("/exports/igmm/eddie/tchandra-lab/scrofts/GenScot_raw/covariates.csv")

#merge with covs based on "id"
covs <- merge(covs, add_covs[,c("id", "bmi", "alc_units")], by="id", all.x=TRUE)

#let's log transform alc_units to reduce skew
covs$alc_units <- log(covs$alc_units + 1) 

#also log transform pack_years
covs$pack_years <- log(covs$pack_years + 1)

#let's scale all
covs$bmi <- scale(covs$bmi)
covs$alc_units <- scale(covs$alc_units)
covs$pack_years <- scale(covs$pack_years)
covs$FEV <- scale(covs$FEV)
covs$height <- scale(covs$height)

#let's scale the additional variables too
covs$p14 <- scale(covs$p14)
covs$p16 <- scale(covs$p16)
covs$p21 <- scale(covs$p21)
covs$mean_senesc <- scale(covs$mean_senesc)

#plot histogram of both to check
#open window for plot
# hist(covs$bmi, breaks=50, main="Histogram of BMI", xlab="BMI")
# hist(covs$alc_units, breaks=50, main="Histogram of Alc Units", xlab="Alc Units")
# hist(covs$pack_years, breaks=50, main="Histogram of Pack Years", xlab="Pack Years", xlim=c(0,3), ylim=c(0,500))
# hist(covs$FEV, breaks=50, main="Histogram of FEV", xlab="FEV")
# hist(covs$height, breaks=50, main="Histogram of Height", xlab="Height")

# Read in diseases
diseases <- read.csv("data/additional_genscot_files/2026-02-17_disease_codes_combined.csv", stringsAsFactors = FALSE)

# List of clocks to loop through
vars <- c("p16", "p21")

#get unique disease names
dis_list <- unique(diseases$Disease)

#let's just focus on "pulm_fibrosis" and "alc_liver"
dis_list <- c("pulm_fibrosis", "liver_alc")

# Create storage lists for results
res.cox.basic <- setNames(vector("list", length(vars)), vars)

res.cox.basic.zph.local <- setNames(vector("list", length(vars)), vars)
res.cox.basic.zph.global <- setNames(vector("list", length(vars)), vars)

# Disease summary table
disease_summary <- list()

# Loop through each disease
for (disease in dis_list) {
  disease_data <- diseases %>% filter(Disease == disease)
  cat(sprintf("Processing disease %s\n", disease))
  
  # Calculate time to incident
  disease_data <- disease_data %>%
    mutate(yod = as.numeric(substr(dt1_ym, 1, 4)),
           mod = as.numeric(substr(dt1_ym, 5, 6)),
           yoa = as.numeric(substr(gs_appt, 1, 4)),
           moa = as.numeric(substr(gs_appt, 5, 6)),
           t_disease = (yod - yoa) + ((mod - moa) / 12)) %>%
    select(-mod, -yoa, -moa)  # Remove excess columns

  # Merge with clock data
  merged_data <- merge(disease_data, covs, by = "id", all.y = TRUE)
  
  # Filter for incident cases
  merged_data_2 <- merged_data %>% filter(incident != 0 | is.na(incident))
  merged_data_2 <- merged_data_2 %>% filter(t_disease > 0 | is.na(t_disease))
  merged_data_2 <- merged_data_2 %>%
    mutate(DP = ifelse(t_disease > 0, 1, 0)) %>%
    filter(is.na(yod) | yod > 1980) 
 
  # Order data and remove duplicates
  merged_data_ord <- merged_data_2[order(merged_data_2$dt1_ym), ]
  merged_data_ord <- merged_data_ord[!duplicated(na.omit(merged_data_ord$id)), ]
  
  # Create event column and time to event
  merged_data_ord <- merged_data_ord %>%
    mutate(event = ifelse(is.na(incident), 0, 1),
           t_event = ifelse(event == 1, t_disease, t_censor))

  # Trim data
  final_data <- merged_data_ord %>% select(-yod, -dt1_ym, -gs_appt)

  # Skip disease if no events
  if (sum(final_data$event == 1) == 0) next

  # Gender filtering
  prop_F <-  length(which(final_data$sex == 1 & final_data$event ==1))/ length(which(final_data$event ==1 ))
  if (prop_F > 0.9) {
    final_data <- final_data %>% filter(sex == 1)
  } else if (prop_F < 0.1) {
    final_data <- final_data %>% filter(sex == 0)
  }

  # Filter to 10 year analyses
  final_data$event[final_data$t_event > 10 & final_data$event == 1] <- 0
  final_data$t_event[final_data$t_event > 10] <- 10 

  # Count events
  n_event <- sum(final_data$event == 1)
  final_data$n_event <- n_event

  # Skip if n_event < 30
  if (n_event < 30) next

  ## Run Cox regression model
  if (prop_F > 0.1 & prop_F < 0.9) {
    # Case 1: prop_F between 0.1 and 0.9
    for (var in vars) {

      #if variable is FEV, then we also want to adjust for height
      if (disease == "pulm_fibrosis") {
        model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ ", var,"+ age + sex + height + FEV + bmi + pack_years"))
      } else {
        model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ ", var," + age + sex + bmi + alc_units + pack_years"))
      }
      
      res.cox.basic[[var]][[disease]] <- coxph(model_formula1, data = final_data)

    }
    
  } else if (prop_F <= 0.1) {
    # Case 2: prop_F less than 0.1
    for (var in vars) {

      if (disease == "pulm_fibrosis") {
        model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ ", var, " + age + height + FEV + bmi + pack_years"))
      } else {
        model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ ", var, " + age + bmi + alc_units + pack_years"))
      }
      
      res.cox.basic[[var]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==0,])

    }
  } else {
    # Case 3: prop_F greater than or equal to 0.9
    for (var in vars) {

      if (disease == "pulm_fibrosis") {
        model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ ", var, " + age + height + FEV + bmi + pack_years"))
      } else {
        model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ ", var, " + age + bmi + alc_units + pack_years"))
      } 
      
      res.cox.basic[[var]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==1,])
    
    }
  }

cat(sprintf("Running zph for %s\n", disease))

if (prop_F > 0.1 & prop_F < 0.9) {

for (var in vars) {
  # Basic model
  res.cox.basic.zph.local[[var]][[disease]]        <- cox.zph(res.cox.basic[[var]][[disease]])$table[1, ]
  res.cox.basic.zph.global[[var]][[disease]]       <- cox.zph(res.cox.basic[[var]][[disease]])$table["GLOBAL", ]

}

  } else {

    for (var in vars) {
      # Basic model
      res.cox.basic.zph.local[[var]][[disease]]        <- cox.zph(res.cox.basic[[var]][[disease]])$table[1, ]
      res.cox.basic.zph.global[[var]][[disease]]       <- cox.zph(res.cox.basic[[var]][[disease]])$table["GLOBAL", ]

    }

}

# Cleanup
rm(final_data)
gc() 

}

process_model_results <- function(model_list, zph_local, zph_global, tag = "model") {
  
  # --- Extract local p-values
  local_p <- lapply(names(zph_local), function(var) {
    tmp <- zph_local[[var]]
    df <- data.frame(do.call('rbind', tmp))
    df$Var <- var
    df$disease <- rownames(df)
    df
  })
  local_p <- do.call(rbind, local_p)
  names(local_p)[1:3] <- c("chisq_local", "df_local", "local_p")
  
  # --- Extract global p-values
  global_p <- lapply(names(zph_global), function(var) {
    tmp <- zph_global[[var]]
    df <- data.frame(do.call('rbind', tmp))
    df$Var <- var
    df$disease <- rownames(df)
    df
  })
  global_p <- do.call(rbind, global_p)
  names(global_p)[1:3] <- c("chisq_global", "df_global", "global_p")
  
  # --- Extract summary output from Cox models
  extract_coefficients <- function(models_list) {
    results <- list()
    for (disease in names(models_list)) {
      disease_coefficients <- list()
      for (var in names(models_list[[disease]])) {
        model <- models_list[[disease]][[var]]
        coeff <- summary(model)$coefficients[1, ]
        n_event <- summary(model)$nevent
        n_total <- summary(model)$n
        disease_coefficients[[var]] <- c(coeff, n_event, n_total)
      }
      results[[disease]] <- do.call(rbind, disease_coefficients)
    }
    final_results <- do.call(rbind, lapply(names(results), function(disease) {
      cbind(Disease = disease, results[[disease]])
    }))
    return(as.data.frame(final_results))
  }
  
  out <- extract_coefficients(model_list)
  names(out)[c(1,3,4,5,6,7,8)] <- c("Var","HR", "SE","Z", "P","N_cases","N_total")
  out$disease <- sub("\\..*", "", rownames(out))
  
  out <- out %>%
    mutate(across(c("HR", "coef", "SE", "P", "Z", "N_cases", "N_total"), as.numeric))
  
  out$LCI <- exp(out$coef - 1.96*out$SE)
  out$UCI <- exp(out$coef + 1.96*out$SE)
  
  out1 <- out[,c("Var","disease","N_cases","N_total","HR","LCI","UCI","P","Z")]
  out1$df <- 1
  
  # --- Merge all results
  out2 <- merge(out1, local_p, by=c("Var", "disease")) 
  out3 <- merge(out2, global_p, by=c("Var", "disease")) 

  
  # --- Optional: add tag
  out3$model <- tag
  return(out3)
}

results_list <- list(
  basic   = list(model=res.cox.basic, zph_local=res.cox.basic.zph.local, zph_global=res.cox.basic.zph.global)
)

all_results <- do.call(rbind, lapply(names(results_list), function(tag) {
  model_info <- results_list[[tag]]
  process_model_results(model_info$model, model_info$zph_local, model_info$zph_global, tag)
}))

#make into a data frame
all_results <- as.data.frame(all_results)

#Change "model" column to be "fully_adjusted"
all_results$model <- "fully_adjusted"

#save as csv 
write.csv(all_results, "results/score_results/alc_liver_pulm_disease_adjusted_cox_results.csv", row.names = FALSE)

#Let's also export the full regression results from res.cox.basic
#we want the whole set of results, so let's export to a text file (it can't be coerced to a data frame easily)
sink("results/score_results/alc_liver_pulm_disease_adjusted_cox_full_results.txt")
print(res.cox.basic)
sink()  
