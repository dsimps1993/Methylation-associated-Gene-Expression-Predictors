
##########  Applying SAGE Predictors to external validation datasets 02/28/26 ########## 

##### GSE179325 - COVID-19 Barturen et al. 2022 #####

library(GEOquery)
library(data.table)

#### Processing GSE179325 methylation and phenodata ####

#Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179325 
GSE179325_meth <- fread("/mnt/qnapnas/data/methylation/GSE179325_COVID_dsimps93/GSE179325_B.Methylation.txt.gz")

#Tidying data
GSE179325_meth <- as.data.frame(GSE179325_meth)

rownames(GSE179325_meth) <- GSE179325_meth$ID
GSE179325_meth <- GSE179325_meth[,-1]

#Downloading metadata
GSE179325_meta <- getGEO("GSE179325",GSEMatrix = T)
GSE179325_meta <- GSE179325_meta$GSE179325_series_matrix.txt.gz

#Table of covid status
table(GSE179325_meta@phenoData@data$`subject status:ch1`)

#Exctracting clean names from metadata to match meth object column names

library(Biobase)
library(stringr)

extracted_ids <- str_replace_all(extracted_ids, "\\[|\\]", "")

identical(extracted_ids, colnames(GSE179325_meth))
#Names and IDs are identical, order is the same, replacing colnames with metadata sample names
colnames(GSE179325_meth) <- rownames(pData(GSE179325_meta))


# construct new ExpressionSet with phenodata and meth in same object
GSE179325_meta_fixed <- ExpressionSet(
  assayData = t(scale(t(as.matrix(GSE179325_meth)))), #scaling meth object as normalization. Rows=CpGs, cols=samples. Scale operates on columns so need to Transmute
  phenoData = phenoData(GSE179325_meta),
  experimentData = experimentData(GSE179325_meta),
  annotation = annotation(GSE179325_meta)
)


#Removing redundant objects
rm(GSE179325_meth,GSE179325_meta)

#Cleaning up pheno variables
GSE179325_meta_fixed$age <- as.numeric(GSE179325_meta_fixed$`age:ch1`)
GSE179325_meta_fixed$gender <- as.factor(GSE179325_meta_fixed$`gender:ch1`)

GSE179325_meta_fixed$disease <- GSE179325_meta_fixed$`disease state:ch1`
GSE179325_meta_fixed$disease[GSE179325_meta_fixed$disease == "SEVERE"] <- "Severe"
GSE179325_meta_fixed$disease[GSE179325_meta_fixed$disease == "MILD"] <- "Mild"
GSE179325_meta_fixed$disease[GSE179325_meta_fixed$disease == "NEGATIVE"] <- "Negative"

GSE179325_meta_fixed$disease <- factor(GSE179325_meta_fixed$disease,  levels = c("Negative", "Mild", "Severe"))


##### looping through files to get SAGE predictors #####

cpg_files <- c("cGAS_tcellrna_coefficients.csv",
               "p14Arf_tcellrna_coefficients.csv",
               "p16Ink4a_tcellrna_coefficients.csv",
               "p21_tcellrna_coefficients.csv")

# Define the base path
base_path <- "/mnt/qnapnas/data/Models/RNA_Predictors/weights/Full_EPICV1_seed65_feb2026/"

# Create full paths
coefficient_paths <- file.path(base_path, cpg_files)

# Optional: load the CSVs into a named list
coefficient_data <- lapply(coefficient_paths, read.csv)
names(coefficient_data) <- gsub("_coefficients.csv", "", cpg_files)


library(dplyr)
library(ggplot2)
library(ggpubr)


#Making separate object for storing predictions and plotting later
plot_data <- GSE179325_meta_fixed@phenoData@data


##### Loop applying all SAGE predictors #####
for (clock_name in names(coefficient_data)) {
  coef_df <- coefficient_data[[clock_name]]
  
  #Intercept standard deviation, the first big value that all the other values are subracted from.
  sumest <- coef_df$coefficient[length(coef_df$coefficient)]
  
  print(clock_name)
  
  #filtering for only CpGs that are present in the dataset (order for CpGs will still be the same too for the pred)
  coef_df <- coef_df[coef_df$X %in% rownames(GSE179325_meta_fixed),]
  
  #list of coef without the intercept value, intercept is already filtered earlier step
  conew <- as.numeric(coef_df$coefficient)
  
  #Extracting methylation of clock CpGs
  cpgdat <- GSE179325_meta_fixed[coef_df$X,]
  
  cpgdat <- cpgdat@assayData$exprs
  
  predout=0
  for(x in 1:ncol(cpgdat)){
    #print(x)
    predout[x] <- sumest + sum(conew*as.numeric(cpgdat[,x]))
  }
  
  #print(predout)
  
  # Add prediction to phenoData
  pred_varname <- paste0(clock_name, "_pred")
  plot_data[[pred_varname]] <- predout
  
  
}


##### Plotting difference between COVID states with multiple testing correction #####


##### Scaling prediction values as Z Score ####




library(ggpubr)
library(dplyr)
library(broom)



senPreds <- c("p14Arf",    "p16Ink4a", "p21")


# Creating scaled (z score) data frame (scaling the output of each predictor)
scaled_df <- scale(plot_data[, paste(senPreds, "tcellrna_pred", sep = "_")])

colnames(scaled_df) <- gsub("_tcellrna_pred", "_Z_score", colnames(scaled_df) )

#Adding scaled predictions to plot_data

plot_data <- cbind(plot_data, scaled_df)

#Creating aggregate senescence z-score
plot_data$senescence_zmean <- rowMeans(scaled_df)


#Creating Dataframe to store pvals
all_pvals <- data.frame(
  Clock = character(),
  P = numeric(),
  stringsAsFactors = FALSE
)

format_p <- function(p) {
  ifelse(p < 0.0001, "****",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*",
                              "ns"))))
}


# Loop getting all pvals - comparing each group to Negative reference
for (clock_name in senPreds) {
  
  pred_varname <- paste0(clock_name, "_Z_score")
  
  # Set Negative as reference level
  plot_data$disease <- factor(plot_data$disease, levels = c("Negative", "Mild", "Severe"))
  
  # Fit the linear model adjusting for sex and age
  model <- lm(as.formula(paste(pred_varname, "~ disease + gender + age")), data = plot_data)
  
  # Get tidy summary to extract coefficients
  tidy_model <- broom::tidy(model)
  
  # Extract p-values for comparisons against Negative
  p_neg_vs_mild <- tidy_model %>% filter(term == "diseaseMild") %>% pull(p.value)
  p_neg_vs_severe <- tidy_model %>% filter(term == "diseaseSevere") %>% pull(p.value)
  
  
  # Prepare a data frame for stat_pvalue_manual
  pval_df <- data.frame(
    group1 = c("Negative", "Negative"),
    group2 = c("Mild", "Severe"),
    Clock = clock_name,
    p = c(p_neg_vs_mild, p_neg_vs_severe),
    p.signif = c(format_p(p_neg_vs_mild), format_p(p_neg_vs_severe)),
    y.position = max(plot_data[[pred_varname]], na.rm = TRUE) * c(1.0, 1.25)
  )
  
  all_pvals <- bind_rows(all_pvals, pval_df)
  
}

# Apply multiple testing correction within each clock separately
all_pvals <- all_pvals %>%
  group_by(Clock) %>%
  mutate(P_BH = p.adjust(p, method = "BH"),
         P_Bonf = p.adjust(p, method = "bonferroni")) %>%
  ungroup()

all_pvals$P_BH_sig <- format_p(all_pvals$P_BH)
all_pvals$P_Bonf_sig <- format_p(all_pvals$P_Bonf)


plot_list <- list()

for (clock_name in senPreds) {
  
  
  covid_colors <- c(
    "Severe" = "#ffb7b2ff",                
    "Mild" = "#f6d6adff",            
    "Negative" = "#fff5baff"  
  )
  
  pred_varname <- paste0(clock_name, "_Z_score")
  
  # Plot
  p1 <- ggboxplot(plot_data,
                  x = "disease",
                  y = pred_varname,
                  fill = "disease",
                  outlier.shape = NA) +
    scale_fill_manual(values = covid_colors) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
    labs(x = "", y = paste(clock_name, " Predicted Expression (Z-Score)")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    stat_pvalue_manual(all_pvals[all_pvals$Clock == clock_name,], label = "P_Bonf_sig", tip.length = 0.01)
  
  # Store the plot in the list
  plot_list[[clock_name]] <- p1
}

# Arrange all plots into one figure
combined_plot <- ggarrange(plotlist = plot_list, 
                           ncol = 1, nrow = 3, 
                           common.legend = TRUE, legend = "bottom")


# Save to PDF
#ggsave("senclocks_MultiAdjP_COVID_Bart2022.pdf", combined_plot,
#       width = 4, height = 10)




##### Merged Senescence Z Score ####

# Creating scaled (z score) data frame 
#scaled_df <- scale(plot_data[, paste(senPreds, "pred", sep = "_")])

plot_data$senescence_zmean <- rowMeans(scaled_df)


all_pvals <- data.frame(
  Clock = character(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# Set Negative as reference level
plot_data$disease <- factor(plot_data$disease, levels = c("Negative", "Mild", "Severe"))

# Fit the linear model adjusting for sex and age
model <- lm(as.formula(paste("senescence_zmean ~ disease + gender + age")), data = plot_data)

# Get tidy summary to extract coefficients
tidy_model <- broom::tidy(model)

# Extract p-values for comparisons against Negative
p_neg_vs_mild <- tidy_model %>% filter(term == "diseaseMild") %>% pull(p.value)
p_neg_vs_severe <- tidy_model %>% filter(term == "diseaseSevere") %>% pull(p.value)


# Prepare a data frame for stat_pvalue_manual
pval_df <- data.frame(
  group1 = c("Negative", "Negative"),
  group2 = c("Mild", "Severe"),
  p = c(p_neg_vs_mild, p_neg_vs_severe),
  p.signif = c(format_p(p_neg_vs_mild), format_p(p_neg_vs_severe)),
  y.position = max(plot_data$senescence_zmean, na.rm = TRUE) * c(1.0, 1.25)
)

all_pvals <- bind_rows(all_pvals, pval_df)


# Apply multiple testing correction
all_pvals$P_BH <- p.adjust(all_pvals$p, method = "BH")
all_pvals$P_Bonf <- p.adjust(all_pvals$p, method = "bonferroni")

all_pvals$P_BH_sig <- format_p(all_pvals$P_BH)
all_pvals$P_Bonf_sig <- format_p(all_pvals$P_Bonf)

#old ones I prefered
#covid_colors <- c(
#  "Severe" = "#e3705f",                
#  "Mild" = "#e6be65",            
#  "Negative" = "#6db9e3"  
#)

covid_colors <- c(
  "Severe" = "#ffb7b2ff",                
  "Mild" = "#f6d6adff",            
  "Negative" = "#fff5baff"  
)


# Plot
p1 <- ggboxplot(plot_data,
                x = "disease",
                y = "senescence_zmean",
                fill = "disease",
                outlier.shape = NA) +
  scale_fill_manual(values = covid_colors) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  #labs(x = "Health Status", y = "Senescence Score (Mean Z-Score)") +
  labs(x = "", y = "Senescence Score (Mean Z-Score)") +
  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_pvalue_manual(all_pvals, label = "P_Bonf_sig", tip.length = 0.01)

p1

plot_list[["MeanZScore"]] <- p1

# Arrange all plots into one figure
combined_plot <- ggarrange(plotlist = plot_list, 
                           ncol = 1, nrow = 4, 
                           common.legend = TRUE, legend = "bottom")


# Save to PDF
ggsave("senclocksFebver_scaled_MultiAdjP_Zscore_COVID_Bart2022_BF.pdf", combined_plot,
       width = 4, height = 14)



##### GSE167202 - COVID-19 Konigsberg et al. 2021 #####


library(GEOquery)
library(data.table)

GSE167202_meth <- fread("/mnt/qnapnas/data/methylation/GSE167202_COVID_dsimps93/GSE167202_ProcessedBetaValues.txt.gz")
#GSE167202_meth <- read.table("/mnt/qnapnas/data/methylation/GSE167202_COVID_dsimps93/GSE167202_B.Methylation.txt.gz", sep = "\t", row.names = 1)

GSE167202_meth <- as.data.frame(GSE167202_meth)

rownames(GSE167202_meth) <- GSE167202_meth$ID_REF

GSE167202_anno <- GSE167202_meth[,1:22]

GSE167202_meth <- GSE167202_meth[,-c(1:22)]


GSE167202_meta <- getGEO("GSE167202",GSEMatrix = T)
GSE167202_meta <- GSE167202_meta$GSE167202_series_matrix.txt.gz

hist(as.numeric(GSE167202_meta@phenoData@data$`age:ch1`))




## getting sample order correct

colnames(GSE167202_meth) %in% GSE167202_meta$description.1
#all true
identical(colnames(GSE167202_meth), GSE167202_meta$description.1)
#wrong order
GSE167202_meth <- GSE167202_meth[,GSE167202_meta$description.1]

identical(colnames(GSE167202_meth), GSE167202_meta$description.1)

colnames(GSE167202_meth) <- GSE167202_meta$geo_accession

# construct new ExpressionSet
GSE167202_meta_fixed <- ExpressionSet(
  assayData = t(scale(t(as.matrix(GSE167202_meth)))), #scaling meth object as normalization. Rows=CpGs, cols=samples. Scale operates on columns so need to Transmute
  #assayData = as.matrix(GSE167202_meth),
  phenoData = phenoData(GSE167202_meta),
  experimentData = experimentData(GSE167202_meta),
  annotation = annotation(GSE167202_meta)
)

GSE167202_meta_fixed$`age:ch1` <- as.numeric(GSE167202_meta_fixed$`age:ch1`)



plot_data2 <- GSE167202_meta_fixed@phenoData@data



#### Loop applying all clocks ###########
for (clock_name in names(coefficient_data)) {
  coef_df <- coefficient_data[[clock_name]]
  
  #Intercept standard deviation, the first big value that all the other values are subracted from.
  sumest <- coef_df$coefficient[length(coef_df$coefficient)]
  
  print(clock_name)
  
  #filtering for only CpGs that are present in the dataset (order for CpGs will still be the same too for the pred)
  coef_df <- coef_df[coef_df$X %in% rownames(GSE167202_meta_fixed),]
  
  #list of coef without the intercept value, intercept is already filtered earlier step
  conew <- as.numeric(coef_df$coefficient)
  
  
  
  cpgdat <- GSE167202_meta_fixed[coef_df$X,]
  
  cpgdat <- cpgdat@assayData$exprs
  
  predout=0
  for(x in 1:ncol(cpgdat)){
    #print(x)
    predout[x] <- sumest + sum(conew*as.numeric(cpgdat[,x]))
  }
  
  #print(predout)
  
  # Add prediction to phenoData
  pred_varname <- paste0(clock_name, "_pred")
  plot_data2[[pred_varname]] <- predout
  
  
}




### Cleaning up COVID terms - merging 

plot_data2$sex <- plot_data2$`Sex:ch1`
plot_data2$age <-plot_data2$`age:ch1`

plot_data2$Covid_Status <- plot_data2$`severity score:ch1`

plot_data2$Covid_Status[plot_data2$`covid_status:ch1` == "negative"] <- "negative"
plot_data2$Covid_Status[plot_data2$`covid_status:ch1` == "other infection"] <- "other infection"



# taking out the ICU intensive negatives
plot_data2 <- plot_data2[!plot_data2$`covid_status:ch1` == "negative" | plot_data2$`admitted_to_icu:ch1` == "No",]



##### Grouping Covid groups by severity #########

plot_data2$Covid_Status2 <- as.character(plot_data2$Covid_Status)

plot_data2$Covid_Status2[plot_data2$Covid_Status2 == "1" | plot_data2$Covid_Status2 == "2"] <- "Mild"
plot_data2$Covid_Status2[plot_data2$Covid_Status2 == "3" | plot_data2$Covid_Status2 == "4"] <- "Severe"

plot_data2$Covid_Status2[plot_data2$Covid_Status2 == "negative"] = "Negative"
plot_data2$Covid_Status2[plot_data2$Covid_Status2 == "other infection"] = "Other Infection"



##### Scaling prediction values as Z Score (without cGAS) ####

library(ggpubr)
library(dplyr)
library(broom)



senPreds <- c("p14Arf",    "p16Ink4a", "p21")


# Creating scaled (z score) data frame 
scaled_df <- scale(plot_data2[, paste(senPreds, "tcellrna_pred", sep = "_")])

colnames(scaled_df) <- gsub("_tcellrna_pred", "_Z_score", colnames(scaled_df) )

#Adding scaled predictions to plot_data2

plot_data2 <- cbind(plot_data2, scaled_df)

#Creating aggregate senescence z-score
plot_data2$senescence_zmean2 <- rowMeans(scaled_df)



all_pvals <- data.frame(
  Clock = character(),
  #Comparison = character(),
  p = numeric(),
  stringsAsFactors = FALSE
)


#double checked that everything operates as relative to first categoric term in covid status, appears correct https://www.sthda.com/english/articles/40-regression-analysis/163-regression-with-categorical-variables-dummy-coding-essentials-in-r/

# Loop getting all pvals
for (clock_name in senPreds) {
  
  pred_varname <- paste0(clock_name, "_Z_score")
  
  # Having negative compared to others, makes it easier and its our main comparison
  plot_data2$Covid_Status2 <- factor(plot_data2$Covid_Status2, levels = c("Negative", "Mild", "Severe", "Other Infection")) 
  
  # Fit the linear model adjusting for sex
  model <- lm(as.formula(paste(pred_varname, "~ Covid_Status2 + sex + age")), data = plot_data2)
  
  # Get tidy summary to extract coefficients
  tidy_model <- broom::tidy(model) ### go back here and double check how youd get the third comparison
  
  # Create a named vector of p-values for each contrast
  p_neg_vs_mild <- tidy_model %>% filter(term == "Covid_Status2Mild") %>% pull(p.value)
  p_neg_vs_severe <- tidy_model %>% filter(term == "Covid_Status2Severe") %>% pull(p.value)
  p_neg_vs_other <- tidy_model %>% filter(term == "Covid_Status2Other Infection") %>% pull(p.value)
  
  
  
  # Prepare a data frame for stat_pvalue_manual
  pval_df <- data.frame(
    group1 = c("Negative", "Negative", "Negative"),
    group2 = c("Mild", "Severe", "Other Infection"),
    Clock = clock_name,
    p = c(p_neg_vs_mild, p_neg_vs_severe, p_neg_vs_other), #correcting names from p adj to p because technically its raw p being extracted (which is correct for BH to be applied to, not adj)
    p.signif = c(format_p(p_neg_vs_mild), format_p(p_neg_vs_severe), format_p(p_neg_vs_other)),
    y.position = max(plot_data2[[pred_varname]], na.rm = TRUE) * c(1.0, 1.4, 2)
  )
  
  all_pvals <- bind_rows(all_pvals, pval_df)
  
  
  # Restore default factor order for plotting
  #plot_data2$Covid_Status <- factor(plot_data2$Covid_Status, levels = c("positive", "other infection", "negative"))
  
}

#note on the pvals resulting from the lm model
#For each coefficient in your model (e.g., for Covid_Statusother infection, age, etc.), R performs a t-test to assess whether that coefficient is significantly different from zero:
#t = estimate/std error then pval is calc from t distribution with degrees of freedom = no of obs - no of predictors - 1

#This is not a t-test between group means directly (like you’d get from t.test()), but rather a test of the partial effect of a predictor in a multiple regression model, adjusting for other covariates.
#But mathematically, it’s still a t-test on the estimated coefficient, just in a multivariable context.


# Apply multiple testing correction within each clock separately
all_pvals <- all_pvals %>%
  group_by(Clock) %>%
  mutate(P_BH = p.adjust(p, method = "BH"),
         P_Bonf = p.adjust(p, method = "bonferroni")) %>%
  ungroup()

all_pvals$P_BH_sig <- format_p(all_pvals$P_BH)
all_pvals$P_Bonf_sig <- format_p(all_pvals$P_Bonf)


plot_list <- list()

for (clock_name in senPreds) {
  
  #covid_colors2 <- c(
  #  "Mild" = "#e6be65",                 # red
  #  "Severe" = "#e3705f",
  #  "Other Infection" = "#4DAF4A",              # blue
  #  "Negative" = "#6db9e3"  # green
  #)
  
  covid_colors2 <- c(
    "Mild" = "#f6d6adff",                 # red
    "Severe" = "#ffb7b2ff",
    "Other Infection" = "#b2e2f2ff",              # blue
    "Negative" = "#fff5baff"  # green
  )
  pred_varname <- paste0(clock_name, "_Z_score")
  
  # Plot
  p2 <- ggboxplot(plot_data2,
                  x = "Covid_Status2",
                  y = pred_varname,
                  fill = "Covid_Status2",
                  outlier.shape = NA) +
    scale_fill_manual(values = covid_colors2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
    #labs(x = "COVID Status", y = paste(clock_name, "\n Predicted Expression (Z-Score)")) +
    labs(x = "", y = "") +
    
    # mult = c(bottom_expansion, top_expansion). 
    # 0.2 adds 20% extra space to the top, giving room for p-values.
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
    
    
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    stat_pvalue_manual(all_pvals[all_pvals$Clock == clock_name,], label = "P_Bonf_sig", tip.length = 0.01)
  
  # Store the plot in the list
  plot_list[[clock_name]] <- p2
}

# Arrange all plots into one figure
combined_plot2 <- ggarrange(plotlist = plot_list, 
                            ncol = 1, nrow = 3, 
                            common.legend = TRUE, legend = "bottom")

# Save to PDF
#ggsave("senclocks_MultiAdjP_COVID3_levs2_june2025Ver3.pdf", combined_plot,
#       width = 7, height = 7)





#### Merged Zmean Sen Score - Levs2 ####

all_pvals <- data.frame(
  Clock = character(),
  #Comparison = character(),
  p = numeric(),
  stringsAsFactors = FALSE
)



# Having negative compared to others, makes it easier and its our main comparison
plot_data2$Covid_Status2 <- factor(plot_data2$Covid_Status2, levels = c("Negative", "Mild", "Severe", "Other Infection")) 

# Fit the linear model adjusting for sex
model <- lm(as.formula(paste("senescence_zmean2 ~ Covid_Status2 + sex + age")), data = plot_data2)

# Get tidy summary to extract coefficients
tidy_model <- broom::tidy(model) ### go back here and double check how youd get the third comparison

# Create a named vector of p-values for each contrast
p_neg_vs_mild <- tidy_model %>% filter(term == "Covid_Status2Mild") %>% pull(p.value)
p_neg_vs_severe <- tidy_model %>% filter(term == "Covid_Status2Severe") %>% pull(p.value)
p_neg_vs_other <- tidy_model %>% filter(term == "Covid_Status2Other Infection") %>% pull(p.value)



# Prepare a data frame for stat_pvalue_manual
pval_df <- data.frame(
  group1 = c("Negative", "Negative", "Negative"),
  group2 = c("Mild", "Severe", "Other Infection"),
  #Clock = clock_name,
  p = c(p_neg_vs_mild, p_neg_vs_severe, p_neg_vs_other), #correcting names from p adj to p because technically its raw p being extracted (which is correct for BH to be applied to, not adj)
  p.signif = c(format_p(p_neg_vs_mild), format_p(p_neg_vs_severe), format_p(p_neg_vs_other)),
  y.position = max(plot_data2$senescence_zmean2, na.rm = TRUE) * c(1.0, 1.4, 2)
)

all_pvals <- bind_rows(all_pvals, pval_df) 

#decided its better to include BH as standard (recommended for more than 1 comparison)
# 3. Apply corrections across all tests (doing both just in case)
all_pvals$P_BH <- p.adjust(all_pvals$p, method = "BH")
all_pvals$P_Bonf <- p.adjust(all_pvals$p, method = "bonferroni")

all_pvals$P_BH_sig <- format_p(all_pvals$P_BH)
all_pvals$P_Bonf_sig <- format_p(all_pvals$P_Bonf)

# Restore default factor order for plotting
plot_data2$Covid_Status <- factor(plot_data2$Covid_Status, levels = c("positive", "other infection", "negative"))


p1 <- ggboxplot(plot_data2, x = "Covid_Status2", y = "senescence_zmean2",
                fill = "Covid_Status2", outlier.shape = NA) + scale_fill_manual(values = covid_colors2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  #labs(x = "COVID Status", y = "Senescence Z-score") +
  labs(x = "", y = "") +
  
  # mult = c(bottom_expansion, top_expansion). 
  # 0.2 adds 20% extra space to the top, giving room for p-values.
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  
  
  theme(axis.text.x = element_blank()) +
  stat_pvalue_manual(all_pvals, label = "P_Bonf_sig", tip.length = 0.01)



plot_list[["MeanZScore"]] <- p1

# Arrange all plots into one figure
combined_plot2 <- ggarrange(plotlist = plot_list, 
                            ncol = 1, nrow = 4, 
                            common.legend = TRUE, legend = "bottom")


# Save to PDF
ggsave("senclocksFebver_Scaled_MultiAdjP_Zscore_COVID_Konig2021_BF.pdf", combined_plot2,
       width = 5, height = 12.5)




##### GSE112611 - Crohn's disease Somineni et al. 2019 #####


library(GEOquery)
library(data.table)

GSE112611_meth <- fread("/mnt/qnapnas/data/methylation/GSE135905_Chrons_DJS/GSE112611_beta_values.txt")
#GSE112611_meth <- read.table("/mnt/qnapnas/data/methylation/GSE112611_COVID_dsimps93/GSE112611_B.Methylation.txt.gz", sep = "\t", row.names = 1)

GSE112611_meth <- as.data.frame(GSE112611_meth)

GSE112611_meth <- GSE112611_meth[,-grep("Detection", colnames(GSE112611_meth))]

rownames(GSE112611_meth) <- GSE112611_meth$ID_REF
GSE112611_meth <- GSE112611_meth[,-1]


GSE112611_meta <- getGEO("GSE112611",GSEMatrix = T)

#epic data in GPL21145 
GSE112611_meta <- GSE112611_meta$GSE112611_series_matrix.txt.gz

hist(as.numeric(GSE112611_meta$`age:ch1`))

View(GSE112611_meta@phenoData@data)

identical(GSE112611_meta$description, colnames(GSE112611_meth))
#true

colnames(GSE112611_meth) <- rownames(pData(GSE112611_meta))

table(GSE112611_meta$characteristics_ch1.1)


###some NAs in the meth to deal with 

colMeans(is.na(GSE112611_meth))  # fraction of NAs per sample
rowMeans(is.na(GSE112611_meth))  # fraction of NAs per CpG

#Simple mean/mode imputation
library(impute)
GSE112611_meth <- apply(GSE112611_meth, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
GSE112611_meth <- t(GSE112611_meth)

GSE112611_meth <- as.data.frame(GSE112611_meth)

# construct new ExpressionSet
GSE112611_meta_fixed <- ExpressionSet(
  #assayData = as.matrix(GSE112611_meth),
  assayData = t(scale(t(as.matrix(GSE112611_meth)))), #scaling meth object as normalization. Rows=CpGs, cols=samples. Scale operates on columns so need to Transmute
  phenoData = phenoData(GSE112611_meta),
  experimentData = experimentData(GSE112611_meta),
  annotation = annotation(GSE112611_meta)
)

GSE112611_meta_fixed$`age:ch1` <- as.numeric(GSE112611_meta_fixed$`age:ch1`)



plot_data3 <- GSE112611_meta_fixed@phenoData@data


#### Loop applying all clocks ###########
for (clock_name in names(coefficient_data)) {
  coef_df <- coefficient_data[[clock_name]]
  
  #Intercept standard deviation, the first big value that all the other values are subracted from.
  sumest <- coef_df$coefficient[length(coef_df$coefficient)]
  
  print(clock_name)
  
  #filtering for only CpGs that are present in the dataset (order for CpGs will still be the same too for the pred)
  coef_df <- coef_df[coef_df$X %in% rownames(GSE112611_meta_fixed),]
  
  #list of coef without the intercept value, intercept is already filtered earlier step
  conew <- as.numeric(coef_df$coefficient)
  
  
  
  cpgdat <- GSE112611_meta_fixed[coef_df$X,]
  
  cpgdat <- cpgdat@assayData$exprs
  
  predout=0
  for(x in 1:ncol(cpgdat)){
    #print(x)
    predout[x] <- sumest + sum(conew*as.numeric(cpgdat[,x]))
  }
  
  #print(predout)
  
  # Add prediction to phenoData
  pred_varname <- paste0(clock_name, "_pred")
  plot_data3[[pred_varname]] <- predout
  
  
}




plot_data3$Status_Time <- plot_data3$`baseline vs follow-up:ch1`

plot_data3$Status_Time <- ifelse( plot_data3$`diagnosis:ch1` == "Crohn's disease", yes = paste("Crohns", plot_data3$Status_Time, sep = "_"),
                                  no = paste("Control", plot_data3$Status_Time, sep = "_"))

plot_data3$gender <- plot_data3$`gender:ch1`
plot_data3$age <- plot_data3$`age:ch1`



##### Scaling prediction values as Z Score ####


library(ggpubr)
library(dplyr)
library(broom)



senPreds <- c("p14Arf",    "p16Ink4a", "p21")


# Creating scaled (z score) data frame 
scaled_df <- scale(plot_data3[, paste(senPreds, "tcellrna_pred", sep = "_")])

colnames(scaled_df) <- gsub("_tcellrna_pred", "_Z_score", colnames(scaled_df) )


#Adding scaled predictions to plot_data

plot_data3 <- cbind(plot_data3, scaled_df)

#Creating aggregate senescence z-score
plot_data3$senescence_zmean2 <- rowMeans(scaled_df)


all_pvals <- data.frame(
  Clock = character(),
  #Comparison = character(),
  P = numeric(),
  stringsAsFactors = FALSE
)


# Loop getting all pvals
for (clock_name in senPreds) {
  
  pred_varname <- paste0(clock_name, "_Z_score")
  
  # Ensure factor levels are correct
  plot_data3$Status_Time <- factor(plot_data3$Status_Time, levels = c("Crohns_FU", "Crohns_BL", "Control_BL"))
  
  # Fit the linear model adjusting for sex
  model <- lm(as.formula(paste(pred_varname, "~ Status_Time + gender + age")), data = plot_data3)
  
  # Get tidy summary to extract coefficients
  tidy_model <- broom::tidy(model) ### go back here and double check how youd get the third comparison
  
  # Create a named vector of p-values for each contrast
  p_Crohns_FU_vs_Crohns_BL <- tidy_model %>% filter(term == "Status_TimeCrohns_BL") %>% pull(p.value)
  p_Control_BL_vs_Crohns_FU <- tidy_model %>% filter(term == "Status_TimeControl_BL") %>% pull(p.value)
  
  # Now relevel so we can get Crohns_BL vs Control_BL
  plot_data3$Status_Time <- relevel(plot_data3$Status_Time, ref = "Control_BL")
  model2 <- lm(as.formula(paste(pred_varname, "~ Status_Time + gender + age")), data = plot_data3)
  tidy_model2 <- broom::tidy(model2)
  p_Crohns_BL_vs_Control_BL <- tidy_model2 %>% filter(term == "Status_TimeCrohns_BL") %>% pull(p.value)
  
  
  # Prepare a data frame for stat_pvalue_manual
  pval_df <- data.frame(
    group1 = c("Control_BL", "Crohns_BL", "Control_BL"),
    group2 = c("Crohns_BL", "Crohns_FU", "Crohns_FU"),
    Clock = clock_name,
    p = c(p_Crohns_BL_vs_Control_BL, p_Crohns_FU_vs_Crohns_BL, p_Control_BL_vs_Crohns_FU),
    p.signif = c(format_p(p_Crohns_BL_vs_Control_BL), format_p(p_Crohns_FU_vs_Crohns_BL), format_p(p_Control_BL_vs_Crohns_FU)),
    y.position = max(plot_data3[[pred_varname]], na.rm = TRUE) * c(1.0, 1.3, 1.5)
  )
  
  all_pvals <- bind_rows(all_pvals, pval_df)
  
  
  # Restore default factor order for plotting
  plot_data3$Status_Time <- factor(plot_data3$Status_Time, levels = c("Crohns_FU",  "Crohns_BL", "Control_BL"))
  
}

# Apply multiple testing correction within each clock separately
all_pvals <- all_pvals %>%
  group_by(Clock) %>%
  mutate(P_BH = p.adjust(p, method = "BH"),
         P_Bonf = p.adjust(p, method = "bonferroni")) %>%
  ungroup()

all_pvals$P_BH_sig <- format_p(all_pvals$P_BH)
all_pvals$P_Bonf_sig <- format_p(all_pvals$P_Bonf)



plot_list <- list()

#Crohns_cols <- c(
#  "Crohns_FU"  = "#009E73", # Teal
#  "Crohns_BL"  = "#CC79A7", # Reddish-Purple
#  "Control_BL" = "#D55E00"  # Vermillion
#)
Crohns_cols <- c(
  "Crohns_FU"  = "#c7ceeaff", # Lavender
  "Crohns_BL"  = "#e1c1f6ff",  # Peach/Pastel Orange
  "Control_BL" = "#fff5baff" # Soft Teal"
)

new_labels <- c("Crohns_FU" = "Crohn's Follow Up", 
                "Crohns_BL" = "Crohn's Baseline", 
                "Control_BL" = "Control")


# Reverse the order to: Control, Baseline, Follow Up
plot_data3$Status_Time <- factor(plot_data3$Status_Time, 
                                 levels = c("Control_BL", "Crohns_BL", "Crohns_FU"))


for (clock_name in senPreds) {
  
  pred_varname <- paste0(clock_name, "_Z_score")
  
  # Plot
  p2 <- ggboxplot(plot_data3,
                  x = "Status_Time",
                  y = pred_varname,
                  fill = "Status_Time",
                  outlier.shape = NA) +
    scale_fill_manual(values = Crohns_cols, labels = new_labels) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
    
    # mult = c(bottom_expansion, top_expansion). 
    # 0.2 adds 20% extra space to the top, giving room for p-values.
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
    
    
    
    #labs(x = "Health Status", y = paste(clock_name, "\n Predicted Expression (Z-Score)")) +
    labs(x = "", y = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    stat_pvalue_manual(all_pvals[all_pvals$Clock == clock_name,], label = "P_Bonf_sig", tip.length = 0.01)
  
  # Store the plot in the list
  plot_list[[clock_name]] <- p2
}

# Arrange all plots into one figure
combined_plot3 <- ggarrange(plotlist = plot_list, 
                            ncol = 1, nrow = 3, 
                            common.legend = TRUE, legend = "bottom")





all_pvals <- data.frame(
  Clock = character(),
  #Comparison = character(),
  P = numeric(),
  stringsAsFactors = FALSE
)


# Ensure factor levels are correct
plot_data3$Status_Time <- factor(plot_data3$Status_Time, levels = c("Crohns_FU", "Crohns_BL", "Control_BL"))

# Fit the linear model adjusting for sex
model <- lm(as.formula(paste("senescence_zmean2 ~ Status_Time + gender + age")), data = plot_data3)

# Get tidy summary to extract coefficients
tidy_model <- broom::tidy(model) ### go back here and double check how youd get the third comparison

# Create a named vector of p-values for each contrast
p_Crohns_FU_vs_Crohns_BL <- tidy_model %>% filter(term == "Status_TimeCrohns_BL") %>% pull(p.value)
p_Control_BL_vs_Crohns_FU <- tidy_model %>% filter(term == "Status_TimeControl_BL") %>% pull(p.value)

# Now relevel so we can get Crohns_BL vs Control_BL
plot_data3$Status_Time <- relevel(plot_data3$Status_Time, ref = "Control_BL")
model2 <- lm(as.formula(paste("senescence_zmean2 ~ Status_Time + gender + age")), data = plot_data3)
tidy_model2 <- broom::tidy(model2)
p_Crohns_BL_vs_Control_BL <- tidy_model2 %>% filter(term == "Status_TimeCrohns_BL") %>% pull(p.value)

# Prepare a data frame for stat_pvalue_manual
pval_df <- data.frame(
  group1 = c("Control_BL", "Crohns_BL", "Control_BL"),
  group2 = c("Crohns_BL", "Crohns_FU", "Crohns_FU"),
  Clock = clock_name,
  p = c(p_Crohns_BL_vs_Control_BL, p_Crohns_FU_vs_Crohns_BL, p_Control_BL_vs_Crohns_FU),
  p.signif = c(format_p(p_Crohns_BL_vs_Control_BL), format_p(p_Crohns_FU_vs_Crohns_BL), format_p(p_Control_BL_vs_Crohns_FU)),
  y.position = max(plot_data3$senescence_zmean, na.rm = TRUE) * c(1.0, 1.3, 1.5)
)

all_pvals <- bind_rows(all_pvals, pval_df)


# Restore default factor order for plotting
plot_data3$Status_Time <- factor(plot_data3$Status_Time, levels = c("Crohns_FU",  "Crohns_BL", "Control_BL"))



# 3. Apply corrections across all tests (doing both just in case)
all_pvals$P_BH <- p.adjust(all_pvals$p, method = "BH")
all_pvals$P_Bonf <- p.adjust(all_pvals$p, method = "bonferroni")



all_pvals$P_BH_sig <- format_p(all_pvals$P_BH)
all_pvals$P_Bonf_sig <- format_p(all_pvals$P_Bonf)

#altering just the p21 y position
#all_pvals[10,"y.position"] <- -0.001
#all_pvals[11,"y.position"] <- -0.01

#reordering the 


new_labels <- c("Crohns_FU" = "Crohn's Follow Up", 
                "Crohns_BL" = "Crohn's Baseline", 
                "Control_BL" = "Control")

# Reverse the order to: Control, Baseline, Follow Up
plot_data3$Status_Time <- factor(plot_data3$Status_Time, 
                                 levels = c("Control_BL", "Crohns_BL", "Crohns_FU"))



# Plot
p2 <- ggboxplot(plot_data3,
                x = "Status_Time",
                y = "senescence_zmean2",
                fill = "Status_Time",
                outlier.shape = NA) +
  scale_fill_manual(values = Crohns_cols, labels = new_labels) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  #labs(x = "Health Status", y = "Senescence Score (Mean Z-Score)") +
  labs(x = "", y = "") +
  
  # mult = c(bottom_expansion, top_expansion). 
  # 0.2 adds 20% extra space to the top, giving room for p-values.
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  
  
  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_pvalue_manual(all_pvals, label = "P_Bonf_sig", tip.length = 0.01)


plot_list[["MeanZScore"]] <- p2

# Arrange all plots into one figure
combined_plot3 <- ggarrange(plotlist = plot_list, 
                            ncol = 1, nrow = 4, 
                            common.legend = TRUE, legend = "bottom")


# Save to PDF
ggsave("senclocksFebver_scaled_MultiAdjP_Zscore_Crohns_BF.pdf", combined_plot3,
       width = 4.8, height = 13)



### All three datasets together


#temp plot 
plotlist2 <- list(combined_plot, combined_plot2, combined_plot3)
combined_plotTem <- ggarrange(plotlist = plotlist2, 
                              ncol = 3, nrow = 1, 
                              common.legend = TRUE, legend = "bottom",label.x = "")

ggsave("senclocksFebver_MultiAdjP_Zscore_COVIDs_BF.pdf", combined_plotTem,
       width = 7, height = 13.1)


