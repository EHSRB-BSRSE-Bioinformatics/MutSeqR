# Small valid VCFs with supplied context, avoiding large genome downloads.
write_depth_vcf <- function(directory, filename, sample = "sample1",
                            format_fields = list(VD = c("1", "2"),
                                                 AD = c("99,1", "198,2")),
                            info_fields = list(), starts = c(10, 20),
                            refs = c("C", "T"), alts = c("A", "C")) {
  n <- length(starts)
  stopifnot(all(lengths(format_fields) == n), all(lengths(info_fields) == n))
  format_headers <- vapply(names(format_fields), function(field) {
    sprintf('##FORMAT=<ID=%s,Number=%s,Type=Integer,Description="Test depth">',
            field, if (field == "AD") "R" else "1")
  }, character(1))
  info_headers <- vapply(names(info_fields), function(field) {
    sprintf('##INFO=<ID=%s,Number=1,Type=%s,Description="Test annotation">',
            field, if (field %in% c("sample", "tag")) "String" else "Integer")
  }, character(1))
  records <- vapply(seq_len(n), function(i) {
    context <- if (refs[i] == "T") "ATA" else "ACA"
    info <- c(paste0("context=", context), vapply(names(info_fields), function(field) {
      paste0(field, "=", info_fields[[field]][i])
    }, character(1)))
    genotype <- paste(vapply(format_fields, `[`, character(1), i), collapse = ":")
    paste("chr1", starts[i], ".", refs[i], alts[i], ".", ".",
          paste(info, collapse = ";"), paste(names(format_fields), collapse = ":"),
          genotype, sep = "\t")
  }, character(1))
  file <- file.path(directory, filename)
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000>",
    '##INFO=<ID=context,Number=1,Type=String,Description="Reference context">',
    info_headers, format_headers,
    paste("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", sample, sep = "\t"),
    records), file)
  file
}

expect_depth_values <- function(dat, total, alt = c(1, 2)) {
  expect_equal(dat$total_depth, total)
  expect_equal(dat$alt_depth, alt)
  expect_equal(dat$ref_depth, total - alt)
  expect_equal(dat$vaf, alt / total)
}
