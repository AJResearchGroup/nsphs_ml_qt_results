setwd("~/GitHubs/nsphs_ml_qt_results/issue_24/")
base_input_filename <- "issue_24_subset"
bed_filename <- "issue_24_subset.bed"
bim_filename <- "issue_24_subset.bim"
fam_filename <- "issue_24_subset.fam"
plink_bin_data <- plinkr::read_plink_bin_data(base_input_filename = base_input_filename)
