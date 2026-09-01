###############################################################################
# kaplan_meier.r
#
# PURPOSE
#   Produces the descriptive survival figures for the senescence scores:
#   Kaplan-Meier curves of 10-year all-cause mortality by score quartile,
#   stratified within 10-year age bands.
#
# INPUTS   (paths relative to project root)
#   data/additional_genscot_files/2024-04-19_deaths.csv
#       One row per participant recorded as deceased (linkage extract).
#   data/additional_genscot_files/covar_temp2.csv
#       QC'd covariates and senescence scores. Requires: id, age, t_censor.
#
# OUTPUTS  (one pair per score x age band)
#   results/score_results/KM_<score>_age<lo>-<hi>.pdf   figure
#   results/score_results/KM_<score>_age<lo>-<hi>.csv   underlying survival
#                                                       estimates, for
#                                                       reproducing the figure
#
# DEPENDENCIES
#   R with survival, survminer, ggplot2 and dplyr.
###############################################################################

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

### read in death data and QCd clock data ###
mortality_data <- read.csv("data/additional_genscot_files/2024-04-19_deaths.csv")
d <- read.csv("data/additional_genscot_files/covar_temp2.csv")

### create 10-year death phenotype and tte ###
temp <- which(d$id %in% mortality_data$id)
d$dead <- 0
d$dead[temp] <- 1
temp2 <- which(d$dead == 1 & d$t_censor > 10)
d$dead[temp2] <- 0
d$t_censor[d$t_censor > 10] <- 10

# List of clocks to loop through
clocks <- c("p14", "p16", "p21", "mean_senesc")
clock_labels <- c("p14" = "p14 score", "p16" = "p16 score",
                  "p21" = "p21 score", "mean_senesc" = "Mean senescence score")

#create age windows based on ten-year intervals starting from the minimum age
age_min <- floor(min(d$age) / 10) * 10
age_max <- ceiling(max(d$age) / 10) * 10
age_windows <- seq(age_min, age_max, by = 10)

#which age window has the most events?
event_counts <- sapply(seq_along(age_windows[-length(age_windows)]), function(i) {
  age_lo <- age_windows[i]
  age_hi <- age_windows[i + 1]
  sum(d$dead[d$age >= age_lo & d$age < age_hi])
})
age_window_labels <- paste0(age_windows[-length(age_windows)], "-", age_windows[-1])
event_summary <- data.frame(age_window = age_window_labels, events = event_counts)
print(event_summary)

# Helper: assign quartile labels
assign_quartiles <- function(x) {
  cuts <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  if (length(unique(cuts)) < 5) return(NULL)
  cut(x, breaks = cuts, labels = paste0("Q", 1:4), include.lowest = TRUE)
}

# Loop over clocks and age windows
for (clock in clocks) {
  for (i in seq_along(age_windows[-length(age_windows)])) {

    age_lo <- age_windows[i]
    age_hi <- age_windows[i + 1]

    # Subset to this age window
    sub <- d[d$age >= age_lo & d$age < age_hi, ]
    if (nrow(sub) < 20 || sum(sub$dead) == 0) next

    # Assign quartiles
    sub$quartile <- assign_quartiles(sub[[clock]])
    if (is.null(sub$quartile)) next
    sub <- sub[!is.na(sub$quartile), ]
    if (nrow(sub) < 10) next

    # Fit KM
    fit <- survfit(Surv(t_censor, dead) ~ quartile, data = sub)

    # Compute dynamic y-axis lower bound
    km_summary <- summary(fit)

    #export km_summary as csv
    km_df <- survminer::surv_summary(fit, data = sub)
    km_df$clock <- clock; km_df$age_lo <- age_lo; km_df$age_hi <- age_hi
    write.csv(km_df, sprintf("results/score_results/KM_%s_age%d-%d.csv", clock, round(age_lo), round(age_hi)), row.names = FALSE)

    y_min <- min(km_summary$surv)

    # Build plot object before opening any device
    p <- suppressWarnings(ggsurvplot(
      fit,
      data = sub,
      palette = c("#2196F3", "#4CAF50", "#FF9800", "#F44336"),
      legend.labs = paste0("Q", 1:4),
      legend.title = clock_labels[[clock]],
      title = sprintf("Age %d-%d", round(age_lo), round(age_hi)),
      xlab = "Time (years)",
      ylab = "Survival probability",
      conf.int = FALSE,
      risk.table = FALSE,
      ylim = c(y_min, 1),
      ggtheme = theme_bw(base_size = 11) +
        theme(
          plot.title      = element_text(size = 11, face = "bold"),
          legend.position = "right",
          legend.text     = element_text(size = 9),
          legend.title    = element_text(size = 9)
        ),
      xlim = c(0, 10),
      break.time.by = 2
    ))

    # Save via ggsave to avoid blank first page
    outfile <- sprintf("results/score_results/KM_%s_age%d-%d.pdf", clock, round(age_lo), round(age_hi))
    ggsave(outfile, plot = p$plot, width = 6, height = 4.5)

    message("Saved: ", outfile)
  }
}