library(tidyverse)
library(reticulate)
reticulate::use_python("/usr/bin/python3")
reticulate::py_config()
library("SigProfilerMatrixGeneratorR")
# SigProfilerMatrixGeneratorR::install('GRCh38') # Don't run this every time...
# Note, to use this you must first install python dependencies.
# e.g., on the linux command line:
# pip install SigProfilerAssignment
# pip install SigProfilerMatrixGenerator
# ...other python tools from Alexandrov lab

project_name <- "Axelson_2022"
project_genome <- "GRCh38"
output_path <- file.path(here::here(),"SigProfiler")
group_var <- "sample"

# Clean data into required format for Alexandrov Lab tools...
signature_data <- as.data.frame(mutation_data) %>%
  dplyr::select(!!ensym(group_var), id, variation_type, seqnames, start, end, ref, alt) %>%
  dplyr::rename("Sample" = group_var,
                "ID" = "id",
                "mut_type" = "variation_type",
                "chrom" = "seqnames",
                "pos_start" = "start",
                "pos_end" = "end") %>%
  dplyr::mutate(Sample = stringr::str_replace_all(Sample, " ", "_")) %>%
  dplyr::mutate(chrom = stringr::str_replace(chrom, "chr", "")) %>%
  dplyr::mutate(Project = project_name) %>%
  dplyr::mutate(Genome = project_genome) %>%
  dplyr::mutate(Type = "SOMATIC") %>%
  dplyr::relocate(Project) %>%
  dplyr::relocate(Genome, .after = ID) %>%
  dplyr::mutate(mut_type = "SNP")

if (!dir.exists(output_path)) { dir.create(output_path) }
write.table(signature_data, file = file.path(output_path,"mutation_data.txt"),
            sep = "\t", row.names = F, quote = F)
                 
signature_matrices <- SigProfilerMatrixGeneratorR(project = project_name,
                                                  genome = project_genome,
                                                  matrix_path = file.path(here::here(),"SigProfiler"),
                                                  plot = T, exome = F, bed_file = NULL,
                                                  chrom_based = F, tsb_stat = T, seqInfo = T,
                                                  cushion=100)

source_python(file.path(here::here(),"inst","signatures.py"))
cosmic_fitR(file.path(output_path,"/output/SBS",paste0(project_name,".SBS96.all")),
            output_path)

sigProfilerExtractorR(file.path(output_path,"/output/SBS",paste0(project_name,".SBS96.all")),
           file.path(output_path,"SigProfilerExtractor"),
           project_genome)
