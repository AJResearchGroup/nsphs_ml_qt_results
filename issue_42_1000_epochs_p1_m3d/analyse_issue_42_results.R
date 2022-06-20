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
) + ggplot2::geom_boxplot()
ggplot2::ggsave(filename = "genotype_concordance_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = model_id, y = nmse, fill = pheno_model_id)
) + ggplot2::scale_y_log10() + ggplot2::geom_boxplot()
ggplot2::ggsave(filename = "nmse_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = model_id, y = r_squared, fill = pheno_model_id)
) + ggplot2::scale_y_continuous(limits = c(0, 1)) + ggplot2::geom_boxplot()
ggplot2::ggsave(filename = "r_squared_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = genotype_concordance, y = nmse, color = model_id, shape = pheno_model_id)
) + ggplot2::geom_point(size = 10) +
  ggplot2::scale_y_log10()
ggplot2::ggsave(filename = "genotype_concordance_to_nmse_per_model_combination.png", width = 7, height = 7)

ggplot2::ggplot(
  t, ggplot2::aes(x = genotype_concordance, y = r_squared, color = model_id, shape = pheno_model_id)
) + ggplot2::geom_point(size = 10) +
  ggplot2::scale_y_continuous(limits = c(0, 1))
ggplot2::ggsave(filename = "genotype_concordance_to_r_squared_per_model_combination.png", width = 7, height = 7)

knitr::kable(t[order(t$genotype_concordance, decreasing = TRUE), ])
knitr::kable(t[order(t$nmse, decreasing = FALSE), ])

knitr::kable(t |> dplyr::filter(model_id == "M1" & pheno_model_id == "p1"))

