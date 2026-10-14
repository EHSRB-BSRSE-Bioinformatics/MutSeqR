#!/usr/bin/env Rscript

devtools::load_all(quiet = TRUE)
devtools::test(stop_on_failure = TRUE)
