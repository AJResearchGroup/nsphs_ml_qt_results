folder_names <- "~/GitHubs/nsphs_ml_qt/scripts_rackham/logs/20220304_issue_122/sim_data_1"
bed_filename <- file.path(folder_names, "sim_data_1.bed")
bim_filename <- file.path(folder_names, "sim_data_1.bim")
fam_filename <- file.path(folder_names, "sim_data_1.fam")
# phe_filename <- file.path(folder_names, "sim_data_1.phe")
data <- plinkr::create_plink_bin_filenames(
  bed_filename = bed_filename,
  bim_filename = bim_filename,
  fam_filename = fam_filename
)
phenotype_data <- plinkr::create_phenotype_data_filename(
  phe_filename = phe_filename
)
phe_table <- plinkr::read_plink_phe_file(phe_filename)

plinkr::assoc_qt(
  assoc_qt_data = plinkr::create_assoc_qt_data(
    data = data,
    phenotype_data = phenotype_data
  ),
  assoc_qt_params = plinkr::create_assoc_qt_params(
    data = data,
    phe_table = phe_table
  ),
  verbose = TRUE
)


df <- plinkr::read_plink_qassoc_file("~/assoc_output.additive.qassoc")
hist(df$P, breaks = seq(0, 1, by = 0.01))
