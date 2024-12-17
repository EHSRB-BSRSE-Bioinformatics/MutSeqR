#' multithread vcf import
#' @description this function will run the import_genomic_vcf function for
#' multiple VCF files, with the ability for multi-threading. 
#' 
# library(future)
# library(future.apply)

# Function to process multiple VCF files
process_vcf_files <- function(vcf_files, use_multithreading = FALSE) {
    if (use_multithreading) {
        # Plan for parallel processing
        plan(multisession)  # Use 'multicore' on Linux or 'multisession' for Windows
    } else {
        plan(sequential)  # Use sequential processing
    }
    
    # Apply the import function to each VCF file
    vcf_data_list <- future_lapply(vcf_files, import_vcf)

    # Optionally combine the results into a single data frame
    combined_data <- do.call(rbind, vcf_data_list)  # Adjust as needed for your data structure
    
    return(combined_data)
}