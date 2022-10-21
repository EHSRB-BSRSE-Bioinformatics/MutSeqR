#!/usr/bin/R
library(ggplot2)
theme_set(theme_bw())
library(tidyverse)

per.nucleotide.data <- import_ds_data(mut_file = "./data/Jonatan_Mutations_in_blood_and_sperm_samples_221021_MM.txt",
                                      rsids = T)

# All data
ggplot(per.nucleotide.data, aes(x=sample, y=-log(vaf))) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = is_known)) + facet_wrap(~tissue)

trinucleotide_frequencies_global <- per.nucleotide.data %>%
  dplyr::group_by(context, subtype) %>%
  dplyr::summarise(depth_per_trinucleotide = sum(total_depth))

trinucleotide_frequencies_dose_group <- per.nucleotide.data %>%
  group_by(context, subtype, Dose) %>%
  summarise(depth_per_trinucleotide = sum(total_depth))

trinucleotide_frequencies_per_mouse <- per.nucleotide.data %>%
  group_by(context, subtype, sample) %>%
  summarise(depth_per_trinucleotide = sum(total_depth))

trinucleotide_counts_per_dose <- per.nucleotide.data %>%
  group_by(context, subtype, Dose) %>%
  summarise(count = sum(mut_depth))
