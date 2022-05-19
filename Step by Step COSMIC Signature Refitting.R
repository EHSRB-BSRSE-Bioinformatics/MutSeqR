# Script to analyze distances between two dose groups of a chemical exposure using trinucleotide mutational patterns
library(ggplot2)
library(tidyverse)
library(spgs)
library(MutationalPatterns)
#library(SomaticSignatures)
library(phyloseq)
library(vegan)
library(IRanges)
library(GenomicRanges)
library(fuzzyjoin)
set.seed(3135) #What is set.seed? What does the number mean?

#Load the count data
PRC_BM_Trinuc_R <- read.delim("PRC_BM_Trinuc_Count_R.txt")
mutation.df <- PRC_BM_Trinuc_R

# Clean column names
colnames(mutation.df) <- c("mutation_type", "0", "6.25", "12.5", "25")

# Remove final row showing totals
mutation.df <- mutation.df[1:96,]

# Make data into long format #column names refer to dose and values refer to freq
mutation.df.long <- tidyr::pivot_longer(mutation.df, cols=2:5, names_to="dose", values_to="frequency")

# Convert things to matrix format
#Dose columns (D1-D3) become values, Mutation Type become row names
mutation.matrix <- as.matrix(mutation.df[,2:5])


row.names(mutation.matrix) <- mutation.df[,1]

plot_96_profile(mutation.matrix, ymax=0.1)

###########Signature Refitting - MutationalPatterns #########################

#retrieve SNV signatures from COSMIC databases
signatures<-get_known_signatures(muttype = "snv", genome = "mm10")

#Create CONTROL signature using 0 dose - convert count data into mutation type proportions and add it into our signature object
sigs_df <- as.data.frame(signatures)
sigs_df$CONTROL <- mutation.df[,"0"]/sum(mutation.df[,"0"])
sigs_new <- as.matrix(sigs_df) #If you want to run this in the Fit, replace "signatures" with "sigs_new"

#Fit to signatures - normal fit
fit_res<-fit_to_signatures(mutation.matrix, signatures)
#with controll signatures
fit_res_w_control<-fit_to_signatures(mutation.matrix, sigs_new)


#visualize
#Show the relative contribution of each signature to the fitted model
plot_contribution(fit_res$contribution, coord_flip=FALSE, mode="absolute")
plot_contribution_heatmap(fit_res$contribution)

#Fit to signatures - Strict fit, iterative removal of signatures
strict_refit<-fit_to_signatures_strict(mutation.matrix, signatures, max_delta=0.004)


fit_res_strict <- strict_refit$fit_res
#visualise the strict fit contribution of signatures
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)
#visualize the order in which signatures were removed during strict fit
fig_list <- strict_refit$sim_decay_fig
fig_list[[1]]
fig_list[[2]]
fig_list[[3]]
fig_list[[4]]



#Bootstrap: indicates the relative stability of your refitted model
contri_boots <- fit_to_signatures_bootstrapped(mutation.matrix,
                                               signatures,
                                               n_boots = 1000,
                                               method = "strict")

plot_bootstrapped_contribution(contri_boots)

#dotplot - color of dots shows  % of iterations in which the sig is found (contribution>0) and the dot size reps the avg conitrbution of the sig
plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "dotplot")
plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "none")
plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "barplot")

#### Table of relative signature contributions

# Calculate relative contributions and generate Table
df <- fit_res$contribution 
df_relative <- apply(df, 2, function(x){x/sum(x)})
df_relative

library(writexl)
df_relative.export<-as.data.frame(df_relative)
write_xlsx(df_relative.export, "Relative sig contribution with control.xlsx")

# Cosine Similarity between reconstructed profile vs original: how well does the refitted model compare to the raw data?
cos_sim_v1<- cos_sim_matrix(fit_res$reconstructed, mutation.matrix)
plot_original_vs_reconstructed(mutation.matrix, fit_res$reconstructed, y_intercept = 0.8)
plot_cosine_heatmap(cos_sim_v1)    
cos_sim_v1

#Cosine similarity between reconstructed profile vs signatures
cos_sim_v2 <- cos_sim_matrix(fit_res$reconstructed, signatures)
plot_cosine_heatmap(cos_sim_v2)

#Cosine Similarities between raw data and Signatures
cos_sim_samples_signatures <- cos_sim_matrix(mutation.matrix, signatures)
#Visualise with a heatmap
plot_cosine_heatmap(cos_sim_samples_signatures, 
                    cluster_rows = TRUE, cluster_cols = TRUE)

#Look at the cosine similarities between samples 
cos_sim_samples <- cos_sim_matrix(mutation.matrix, mutation.matrix)
plot_cosine_heatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE)
plot_cosine_heatmap(cos_sim_samples, cluster_rows = FALSE, cluster_cols = FALSE)
library(writexl)
cosinesim.export<-as.data.frame(cos_sim_samples)
write_xlsx(cosinesim.export, "Cosine similarity mut spectra vs signatures with control.xlsx")


#per mouse
pheatmap::pheatmap(t(pairwise_cosine_similarity_sigs),
                   cluster_cols=F,
                   annotation=sampledata %>% dplyr::select(Dose))
