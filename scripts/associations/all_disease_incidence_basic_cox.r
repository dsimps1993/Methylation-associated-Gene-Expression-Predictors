###############################################################################
# all_disease_incidence_basic_cox.r
#
# PURPOSE
#   Phenome-wide Cox survival scan: tests each methylation-based senescence
#   score (p14, p16, p21, mean_senesc) against 10-year incidence of every
#   disease available in the Generation Scotland linkage data, and repeats
#   the scan within strata defined by sex and by smoking status.
#
# MODELS
#   Per score x disease, the following are fitted (subject to the sex rules
#   described below):
#     basic            Surv(t_event, event) ~ scale(score) + age + sex
#     basic_male       same, restricted to men,   sex term dropped
#     basic_female     same, restricted to women, sex term dropped
#     basic_smoke      basic model, ever-smokers only
#     basic_nosmoke    basic model, never-smokers only
#     basic_sex_int    Surv(...) ~ scale(score) * sex   + age
#     basic_smoke_int  Surv(...) ~ scale(score) * smoke + age + sex
#
#   Scores are z-scaled so hazard ratios are per standard deviation. These
#   are "basic" (minimally adjusted) models: age and sex only. The stratified
#   fits estimate the effect separately in each group; the interaction models
#   formally test whether the effect differs between groups, and only the
#   interaction coefficient is retained from those.
#
# COHORT CONSTRUCTION (per disease, inside the main loop)
#   1. Convert the YYYYMM diagnosis and appointment dates to a decimal
#      time-to-diagnosis in years.
#   2. Right-join onto the full covariate file so that participants with no
#      record of the disease are retained as controls.
#   3. Keep only incident cases (diagnosis after baseline); prevalent cases
#      are excluded rather than treated as controls.
#   4. Retain one row per participant, keeping the earliest diagnosis.
#   5. Administratively censor at 10 years.
#
# INCLUSION RULES
#   - Sex-specific diseases are detected empirically rather than hard-coded:
#     if >90% or <10% of cases are female, the cohort is restricted to that
#     sex and the sex term is dropped from the model.
#   - Diseases with fewer than 30 incident events are skipped, as hazard
#     ratios from smaller case counts are unstable.
#
# INPUTS   (paths relative to project root)
#   data/additional_genscot_files/covar_temp2.csv
#       QC'd covariates and senescence scores. Requires: id, age, sex,
#       pack_years, t_censor.
#   data/additional_genscot_files/2026-02-17_disease_codes_combined.csv
#       Long-format linkage extract; one row per participant-diagnosis.
#       Requires: id, Disease, dt1_ym (diagnosis YYYYMM), gs_appt (baseline
#       appointment YYYYMM), incident.
#
# OUTPUTS
#   results/score_results/Sen_Results_10Aug2026.csv
#       Main and stratified effects: HR, 95% CI, Z, P, counts, proportional
#       hazards diagnostics, and a `model` column identifying the stratum.
#   results/score_results/Sen_Interaction_Results_10Aug2026.csv
#       Score x sex and score x smoking interaction terms only.
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

# create binary smoking (NB 395 missing pack year info)
covs$smoke <- ifelse(covs$pack_years == 0, 0, 1)

# Read in diseases
diseases <- read.csv("data/additional_genscot_files/2026-02-17_disease_codes_combined.csv", stringsAsFactors = FALSE)

# List of clocks to loop through
clocks <- c("p14", "p16", "p21", "mean_senesc")

# Get unique diseases
dis_list <- unique(diseases$Disease)

# Check number of unique diseases
cat("Number of unique diseases:", length(dis_list), "\n")

# Create storage lists for results
res.cox.basic <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.male <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.female <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.sex.int <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.smoker <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.nosmoke <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.smoke.int <- setNames(vector("list", length(clocks)), clocks)

res.cox.basic.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.male.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.male.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.female.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.female.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.smoker.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.smoker.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.nosmoke.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.basic.nosmoke.zph.global <- setNames(vector("list", length(clocks)), clocks)



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
    for (clock in clocks) {
      model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + sex"))
      model_formula2 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age"))
      model_formula3 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") * sex + age"))
      model_formula4 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") * smoke + age + sex"))

      res.cox.basic[[clock]][[disease]] <- coxph(model_formula1, data = final_data)
      res.cox.basic.male[[clock]][[disease]] <- coxph(model_formula2, data = final_data[final_data$sex==0,])
      res.cox.basic.female[[clock]][[disease]] <- coxph(model_formula2, data = final_data[final_data$sex==1,])     
      res.cox.basic.sex.int[[clock]][[disease]] <- coxph(model_formula3, data = final_data)
      res.cox.basic.smoker[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$smoke==1,])
      res.cox.basic.nosmoke[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$smoke==0,])
      res.cox.basic.smoke.int[[clock]][[disease]] <- coxph(model_formula4, data = final_data)

     

    }
  } else if (prop_F <= 0.1) {
    # Case 2: prop_F less than 0.1
    for (clock in clocks) {

      model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age"))
      model_formula2 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") * smoke + age"))
      
      res.cox.basic[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==0,])
      res.cox.basic.smoker[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==0 & final_data$smoke==1,])
      res.cox.basic.nosmoke[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==0 & final_data$smoke==0,])
      res.cox.basic.smoke.int[[clock]][[disease]] <- coxph(model_formula2, data = final_data[final_data$sex==0,])



    }
  } else {
    # Case 3: prop_F greater than or equal to 0.9
    for (clock in clocks) {
      model_formula1 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age"))
      model_formula2 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") * smoke + age"))
      
      res.cox.basic[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==1,])
      res.cox.basic.smoker[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==1 & final_data$smoke==1,])
      res.cox.basic.nosmoke[[clock]][[disease]] <- coxph(model_formula1, data = final_data[final_data$sex==1 & final_data$smoke==0,])
      res.cox.basic.smoke.int[[clock]][[disease]] <- coxph(model_formula2, data = final_data[final_data$sex==1,])



    }
  }

cat(sprintf("Running zph for %s\n", disease))


if (prop_F > 0.1 & prop_F < 0.9) {

for (clock in clocks) {
  # Basic model
  res.cox.basic.zph.local[[clock]][[disease]]        <- cox.zph(res.cox.basic[[clock]][[disease]])$table[1, ]
  res.cox.basic.zph.global[[clock]][[disease]]       <- cox.zph(res.cox.basic[[clock]][[disease]])$table["GLOBAL", ]

  # Basic male/female/smoker/nosmoke
  res.cox.basic.male.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.basic.male[[clock]][[disease]])$table[1, ]
  res.cox.basic.male.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.basic.male[[clock]][[disease]])$table["GLOBAL", ]
  
  res.cox.basic.female.zph.local[[clock]][[disease]] <- cox.zph(res.cox.basic.female[[clock]][[disease]])$table[1, ]
  res.cox.basic.female.zph.global[[clock]][[disease]]<- cox.zph(res.cox.basic.female[[clock]][[disease]])$table["GLOBAL", ]
  
  res.cox.basic.smoker.zph.local[[clock]][[disease]] <- cox.zph(res.cox.basic.smoker[[clock]][[disease]])$table[1, ]
  res.cox.basic.smoker.zph.global[[clock]][[disease]]<- cox.zph(res.cox.basic.smoker[[clock]][[disease]])$table["GLOBAL", ]
  
  res.cox.basic.nosmoke.zph.local[[clock]][[disease]]<- cox.zph(res.cox.basic.nosmoke[[clock]][[disease]])$table[1, ]
  res.cox.basic.nosmoke.zph.global[[clock]][[disease]]<- cox.zph(res.cox.basic.nosmoke[[clock]][[disease]])$table["GLOBAL", ]


}

  } else {

    for (clock in clocks) {
  # Basic model
  res.cox.basic.zph.local[[clock]][[disease]]        <- cox.zph(res.cox.basic[[clock]][[disease]])$table[1, ]
  res.cox.basic.zph.global[[clock]][[disease]]       <- cox.zph(res.cox.basic[[clock]][[disease]])$table["GLOBAL", ]

  res.cox.basic.smoker.zph.local[[clock]][[disease]] <- cox.zph(res.cox.basic.smoker[[clock]][[disease]])$table[1, ]
  res.cox.basic.smoker.zph.global[[clock]][[disease]]<- cox.zph(res.cox.basic.smoker[[clock]][[disease]])$table["GLOBAL", ]
  
  res.cox.basic.nosmoke.zph.local[[clock]][[disease]]<- cox.zph(res.cox.basic.nosmoke[[clock]][[disease]])$table[1, ]
  res.cox.basic.nosmoke.zph.global[[clock]][[disease]]<- cox.zph(res.cox.basic.nosmoke[[clock]][[disease]])$table["GLOBAL", ]

}

}


  # Cleanup
  rm(final_data)
  gc() 
}


process_model_results <- function(model_list, zph_local, zph_global, tag = "model") {
  
  # --- Extract local p-values
  local_p <- lapply(names(zph_local), function(clock) {
    tmp <- zph_local[[clock]]
    df <- data.frame(do.call('rbind', tmp))
    df$Clock <- clock
    df$disease <- rownames(df)
    df
  })
  local_p <- do.call(rbind, local_p)
  names(local_p)[1:3] <- c("chisq_local", "df_local", "local_p")
  
  # --- Extract global p-values
  global_p <- lapply(names(zph_global), function(clock) {
    tmp <- zph_global[[clock]]
    df <- data.frame(do.call('rbind', tmp))
    df$Clock <- clock
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
      for (clock in names(models_list[[disease]])) {
        model <- models_list[[disease]][[clock]]
        coeff <- summary(model)$coefficients[1, ]
        n_event <- summary(model)$nevent
        n_total <- summary(model)$n
        disease_coefficients[[clock]] <- c(coeff, n_event, n_total)
      }
      results[[disease]] <- do.call(rbind, disease_coefficients)
    }
    final_results <- do.call(rbind, lapply(names(results), function(disease) {
      cbind(Disease = disease, results[[disease]])
    }))
    return(as.data.frame(final_results))
  }
  
  out <- extract_coefficients(model_list)
  names(out)[c(1,3,4,5,6,7,8)] <- c("Clock","HR", "SE","Z", "P","N_cases","N_total")
  out$disease <- sub("\\..*", "", rownames(out))
  
  out <- out %>%
    mutate(across(c("HR", "coef", "SE", "P", "Z", "N_cases", "N_total"), as.numeric))
  
  out$LCI <- exp(out$coef - 1.96*out$SE)
  out$UCI <- exp(out$coef + 1.96*out$SE)
  
  out1 <- out[,c("Clock","disease","N_cases","N_total","HR","LCI","UCI","P","Z")]
  out1$df <- 1
  
  # --- Merge all results
  out2 <- merge(out1, local_p, by=c("Clock", "disease")) 
  out3 <- merge(out2, global_p, by=c("Clock", "disease")) 

  
  # --- Optional: add tag
  out3$model <- tag
  return(out3)
}




results_list <- list(
  basic   = list(model=res.cox.basic, zph_local=res.cox.basic.zph.local, zph_global=res.cox.basic.zph.global),
  basic_male   = list(model=res.cox.basic.male, zph_local=res.cox.basic.male.zph.local, zph_global=res.cox.basic.male.zph.global),
  basic_female = list(model=res.cox.basic.female, zph_local=res.cox.basic.female.zph.local, zph_global=res.cox.basic.female.zph.global),
  basic_smoke   = list(model=res.cox.basic.smoker, zph_local=res.cox.basic.smoker.zph.local, zph_global=res.cox.basic.smoker.zph.global),
  basic_nosmoke = list(model=res.cox.basic.nosmoke, zph_local=res.cox.basic.nosmoke.zph.local, zph_global=res.cox.basic.nosmoke.zph.global)

)

all_results <- do.call(rbind, lapply(names(results_list), function(tag) {
  model_info <- results_list[[tag]]
  process_model_results(model_info$model, model_info$zph_local, model_info$zph_global, tag)
}))






### now to extract the interaction terms for sex and smoking

extract_interaction_term_only <- function(models_list, model_tag = NA, interaction_pattern = ":") {
  results <- list()

  for (clock in names(models_list)) {
    clock_coeffs <- list()
    
    for (disease in names(models_list[[clock]])) {
      model <- models_list[[clock]][[disease]]
      coef_table <- summary(model)$coefficients
      
      # Find the interaction term row
      interaction_row <- coef_table[grep(interaction_pattern, rownames(coef_table)), , drop = FALSE]
      
      if (nrow(interaction_row) == 1) {
        coef_val  <- as.numeric(interaction_row[1, "coef"])
        HR        <- as.numeric(interaction_row[1, "exp(coef)"])
        SE        <- as.numeric(interaction_row[1, "se(coef)"])
        P         <- as.numeric(interaction_row[1, "Pr(>|z|)"])
        Z         <- as.numeric(interaction_row[1, "z"])
	   df         <- 1
        n_event   <- summary(model)$nevent
        n_total   <- summary(model)$n
        
        clock_coeffs[[disease]] <- data.frame(
          Clock = clock,
          disease = disease,
          N_cases = n_event,
          N_total = n_total,
          HR = HR,
          coef = coef_val,
          SE = SE,
          P = P,
          Z = Z,
          df = df
        )
      }
    }
    
    if (length(clock_coeffs) > 0) {
      results[[clock]] <- do.call(rbind, clock_coeffs)
    }
  }

  # Combine all into one dataframe
  final <- do.call(rbind, results)

  # Add confidence intervals
  final$LCI <- exp(final$coef - 1.96 * final$SE)
  final$UCI <- exp(final$coef + 1.96 * final$SE)
  final$model <- model_tag

  # Rearrange columns
  final <- final[, c("Clock", "disease", "N_cases", "N_total", "HR", "LCI", "UCI", "P","Z","df", "model")]

  return(final)
}


interaction_models <- list(
  basic_sex_int = res.cox.basic.sex.int,
  basic_smoke_int = res.cox.basic.smoke.int
)

# Run extraction for each
interaction_results <- do.call(rbind, lapply(names(interaction_models), function(tag) {
  extract_interaction_term_only(interaction_models[[tag]], model_tag = tag)
}))

write.csv(all_results, "results/score_results/Sen_Results_10Aug2026.csv", row.names=F)
write.csv(interaction_results, "results/score_results/Sen_Interaction_Results_10Aug2026.csv", row.names=F)
