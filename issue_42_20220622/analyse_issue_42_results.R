all_root_folders <- list.dirs(recursive = FALSE)
all_folders <- list.dirs(recursive = TRUE)

successful_root_folders <- sort(unique(dirname(dirname(stringr::str_subset(all_folders, "pheno_weights")))))
testthat::expect_equal(29, length(successful_root_folders))

matches <- stringr::str_match(successful_root_folders, "_(M.*)_(p.)_(1.*)_ae")
model_ids_that_work <- unique(sort(matches[, 2]))
testthat::expect_equal(4, length(model_ids_that_work))

pheno_model_ids_that_work <- unique(sort(matches[, 3]))
testthat::expect_equal(2, length(pheno_model_ids_that_work))

window_kbs_that_work <- unique(sort(matches[, 4]))
testthat::expect_equal(4, length(window_kbs_that_work))

#############################################################################
t_combinations <- tidyr::expand_grid(
  model_id = model_ids_that_work,
  pheno_model_id = pheno_model_ids_that_work,
  window_kb = window_kbs_that_work
)
t_combinations$folder_name <- paste0(
    "data_issue_42_", t_combinations$model_id,
    "_", t_combinations$pheno_model_id, "_", t_combinations$window_kb,
    "_ae"
)
t_combinations$log_filename_25 <- paste0(
  "25_run_issue_42_", t_combinations$model_id,
  "_", t_combinations$pheno_model_id, "_", t_combinations$window_kb,
  ".log"
)
t_combinations$folder_exists <- dir.exists(t_combinations$folder_name)
t_combinations$log_filename_25_exists <- file.exists(t_combinations$log_filename_25)
#############################################################################
# Explore when dir does not exist, but log file does
#############################################################################
if ("caused by" == "concurrency") {
  is_interesting <- t_combinations$log_filename_25_exists & !t_combinations$folder_exists
  interesting_log_files <- t_combinations$log_filename_25[is_interesting]
  readLines(interesting_log_files)
  interesting_log_files
  file.copy(from = interesting_log_files, to = "~/25_run_issue_42_M1_p1_100.log")
}
#############################################################################
# Genotype concordance in time
#############################################################################
genotype_concordance_in_time_list <- list()
for (i in seq_len(nrow(t_combinations))) {
  this_row <- t_combinations[i, ]
  filename <- file.path(
    this_row$folder_name,
    "genotype_concordances.csv"
  )
  if (!file.exists(filename)) {
    message(filename, " does not exist")
    next
  }
  t_genotype_concordances_in_time <- gcaer::read_genotype_concordances_file(
    filename
  )
  t_genotype_concordances_in_time$model_id <- this_row$model_id
  t_genotype_concordances_in_time$pheno_model_id <- this_row$pheno_model_id
  t_genotype_concordances_in_time$window_kb <- this_row$window_kb
  genotype_concordance_in_time_list[[i]] <- t_genotype_concordances_in_time
}
t_genotype_concordances_in_time <- dplyr::bind_rows(genotype_concordance_in_time_list)

window_kb_to_data <- function(window_kb) {
  if (window_kb == 1) return("1 kb, 2 SNPs")
  if (window_kb == 10) return("10 kb, 5 SNPs")
  if (window_kb == 100) return("100 kb, 5 SNPs")
  if (window_kb == 1000) return("1000 kb, 6 SNPs")
  stop("Unknown 'window_kb': ", window_kb)
}

t_genotype_concordances_in_time$genetic_data <- Vectorize(window_kb_to_data)(t_genotype_concordances_in_time$window_kb)

ggplot2::ggplot(
  t_genotype_concordances_in_time,
  ggplot2::aes(x = epoch, y = genotype_concordance, color = genetic_data)
) + ggplot2::geom_line(size = 1) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  ggplot2::facet_grid(model_id ~ pheno_model_id) +
  gcaer::get_gcaer_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
  )
ggplot2::ggsave(filename = "issue_42_genotype_concordances_in_time.png", width = 7, height = 7)


###
### OLD
###
t <- tidyr::expand_grid(
  model_id = model_ids_that_work,
  pheno_model_id = pheno_model_ids_that_work,
  window_kb = window_kbs_that_work,
  genotype_concordance = NA,
  nmse = NA,
  r_squared = NA
)
for (row_index in seq_len(nrow(t))) {
  message(row_index)
  model_id <- t$model_id[row_index]
  pheno_model_id <- t$pheno_model_id[row_index]
  window_kb <- t$window_kb[row_index]
  folder_name <- paste0(
    "data_issue_42_", model_id,
    "_", pheno_model_id, "_", window_kb,
    "_ae"
  )
  if (!dir.exists(folder_name)) next
  testthat::expect_true(dir.exists(folder_name))

  # GC
  filename <- file.path(folder_name, "genotype_concordances.csv")
  if (!file.exists(filename)) {
    message(filename)
    next
  }
  testthat::expect_true(file.exists(filename))
  t$genotype_concordance[row_index] <- tail(
    gcaer::read_genotype_concordances_file(filename), n = 1
  )$genotype_concordance

  # NMSE
  filename <- file.path(folder_name, "nmse_in_time.csv")
  if (!file.exists(filename)) {
    message(filename)
    next
  }
  testthat::expect_true(file.exists(filename))
  t$nmse[row_index] <- tail(
    gcaer::read_nmse_in_time_file(filename), n = 1
  )$nmse

  # r_squared
  filename <- file.path(folder_name, "r_squared_in_time.csv")
  if (!file.exists(filename)) {
    message(filename)
    next
  }
  testthat::expect_true(file.exists(filename))
  t$r_squared[row_index] <- tail(
    gcaer::read_r_squared_in_time_file(filename), n = 1
  )$r_squared
}

ggplot2::ggplot(
  t, ggplot2::aes(x = model_id, y = genotype_concordance, fill = pheno_model_id)
) + ggplot2::geom_col(color = "black", position = "dodge") +
  ggplot2::facet_grid(window_kb ~ .) +
  ggplot2::labs(
    title = "Genotype concordance per model per SNP selection window",
    caption = paste0(
      "Rows are the SNP selection window centered around rs4819959, ",
      "the main hit of IL-17RA (panel CVD3_105_IL-17RA)"
    )
  )
ggplot2::ggsave(filename = "genotype_concordance_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = model_id, y = nmse, fill = pheno_model_id)
) + ggplot2::scale_y_log10("log of NMSE (less is better)") +
  ggplot2::geom_col(color = "black", position = "dodge") +
  ggplot2::facet_grid(window_kb ~ .) +
  ggplot2::labs(
    title = "NMSE per model per SNP selection window",
    caption = paste0(
      "Rows are the SNP selection window centered around rs4819959, ",
      "the main hit of IL-17RA (panel CVD3_105_IL-17RA)"
    )
  )
ggplot2::ggsave(filename = "nmse_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = model_id, y = r_squared, fill = pheno_model_id)
) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::geom_col(color = "black", position = "dodge") +
  ggplot2::facet_grid(window_kb ~ .) +
  ggplot2::labs(
    title = "R squared per model per SNP selection window",
    caption = paste0(
      "Rows are the SNP selection window centered around rs4819959, ",
      "the main hit of IL-17RA (panel CVD3_105_IL-17RA)"
    )
  )
ggplot2::ggsave(filename = "r_squared_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = genotype_concordance, y = nmse, color = model_id, shape = pheno_model_id)
) + ggplot2::geom_point(size = 10) +
  ggplot2::scale_y_log10() +
  ggplot2::facet_grid(window_kb ~ .) +
  ggplot2::labs(
    title = "Genotype concordance versus NMSE per SNP selection window",
    caption = paste0(
      "Rows are the SNP selection window centered around rs4819959, ",
      "the main hit of IL-17RA (panel CVD3_105_IL-17RA)"
    )
  )
ggplot2::ggsave(filename = "genotype_concordance_to_nmse_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = genotype_concordance, y = r_squared, color = model_id, shape = pheno_model_id)
) + ggplot2::geom_point(size = 10) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::facet_grid(window_kb ~ .) +
  ggplot2::labs(
    title = "Genotype concordance versus R squared per SNP selection window",
    caption = paste0(
      "Rows are the SNP selection window centered around rs4819959, ",
      "the main hit of IL-17RA (panel CVD3_105_IL-17RA)"
    )
  )
ggplot2::ggsave(filename = "genotype_concordance_to_r_squared_per_model_combination.png", width = 7, height = 7)

knitr::kable(t[order(t$genotype_concordance, decreasing = TRUE), ])
knitr::kable(t[order(t$nmse), ])
knitr::kable(t[order(t$r_squared, decreasing = TRUE), ])

knitr::kable(t |> dplyr::filter(model_id == "M1" & pheno_model_id == "p1"))

