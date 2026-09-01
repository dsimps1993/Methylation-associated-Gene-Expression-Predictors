"""
correlation_plot.py
 
PURPOSE
    Descriptive figure showing how the methylation-based senescence scores
    (p14, p16, p21 and their mean) relate to each other and to the panel of
    published epigenetic clocks in Generation Scotland.
 
INPUTS   (paths relative to this script)
    ../../data/additional_genscot_files/covar_temp2.csv
        QC'd covariates and senescence scores. Keyed on DNAm_ID.
    ../../data/additional_genscot_files/GS_clocks.csv
        Published epigenetic clock estimates. Keyed on DNAm_ID.
 
OUTPUT
    ../results/correlation_matrix_clustered.pdf
        The clustered heatmap, at 600 dpi. The other three figures are
        exploratory and display to screen only.
 
NOTES
    Pearson correlations on complete pairs (pandas .corr() default). No age
    adjustment is applied, so the strong positive block among the age
    predictors partly reflects their shared dependence on chronological age.
 
DEPENDENCIES
    pandas, seaborn, matplotlib.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

d = pd.read_csv("../../data/additional_genscot_files/covar_temp2.csv")
additional_clocks = pd.read_csv("../../data/additional_genscot_files/GS_clocks.csv")

#merge clock data with covariate data
d = pd.merge(d, additional_clocks, on="DNAm_ID")


clocks = ["p14", "p16", "p21", "mean_senesc", "Horvathv1" , "Hannum", "Lin", 
          "PhenoAge", "YingCausAge", "YingDamAge", "YingAdaptAge", "Horvathv2", "Zhang_10", 
          "DunedinPoAm38", "DunedinPACE", "DNAmGrimAge", "DNAmGrimAge.1", "DNAmTL"]

# Create a correlation matrix
corr_matrix = d[clocks].corr()

# Plot the correlation matrix
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix, annot=False, cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Correlation Matrix of Clocks')
plt.show()

#as above, but clustered
sns.clustermap(corr_matrix, annot=False, cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Clustered Correlation Matrix of Clocks')
plt.show()

#Use the order from the clustered correlation matrix to plot a heatmap without the dendrogram
clocks_clustered = corr_matrix.columns[sns.clustermap(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1).dendrogram_col.reordered_ind]  
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix.loc[clocks_clustered, clocks_clustered], annot=False, cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Correlation Matrix of Clocks (Clustered)')
#export at 600 dpi pdf
plt.savefig('../results/correlation_matrix_clustered.pdf', dpi=600, bbox_inches='tight')
plt.show()


#as above, but without DNAmTL
clocks_no_tl = ["p14", "p16", "p21", "mean_senesc", "Horvathv1" , "Hannum", "Lin", 
          "PhenoAge", "YingCausAge", "YingDamAge", "YingAdaptAge", "Horvathv2", "Zhang_10", 
          "DunedinPoAm38", "DunedinPACE", "DNAmGrimAge", "DNAmGrimAge.1"]
corr_matrix_no_tl = d[clocks_no_tl].corr()
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix_no_tl, annot=False, cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Correlation Matrix of Clocks (without DNAmTL)')
plt.show()
