# Analysis of:
#
# * (0. Functions needed for the analysis)
# * 1. NMSE
# * 2. Genotype concordance
# * 3. Runtime
#
################################################################################
# * 0. Functions needed for the analysis
################################################################################
issue_to_n_markers <- function(issues) {
  markers <- issues
  markers[issues == 28] <- "one"
  markers[issues == 29] <- "more"
  markers
}
testthat::expect_equal(issue_to_n_markers(28), "one")
testthat::expect_equal(issue_to_n_markers(29), "more")
testthat::expect_identical(issue_to_n_markers(c(28, 29)), c("one", "more"))

folder_names <- paste0("issue_", c(28, 29), "_1000_epochs_p1_m3d")
testthat::expect_true(all(dir.exists(folder_names)))
testthat::expect_identical(unique(folder_names), folder_names)

################################################################################
# * 1. NMSE
################################################################################
if (nchar("analysis_about_nsme")) {
  nmse_in_time_filenames_table <- tibble::tibble(
    filename = c(
      list.files(folder_names, pattern = "nmse_in_time.csv", recursive = TRUE, full.names = TRUE)
    )
  )
  nmse_in_times_list <- list()
  for (row_index in seq_len(nrow(nmse_in_time_filenames_table))) {
    nmse_in_time_filename <- nmse_in_time_filenames_table$filename[row_index]
    t <- gcaer::read_nmse_in_time_file(
      nmse_in_time_filename = nmse_in_time_filename
    )
    matches <- stringr::str_match(
      nmse_in_time_filename,
      "issue_([[:digit:]]{2})_1000_epochs_p([[:digit:]]{1})/data_issue_([[:digit:]]{2})_([[:digit:]]{1,4})_ae"
    )
    testthat::expect_true(!is.na(matches)[1][1])
    testthat::expect_equal(matches[1, 2], matches[1, 4])
    issue <- matches[1, 2]
    t$n_markers <- issue_to_n_markers(as.numeric(issue))
    t$phenotype_model <- as.numeric(matches[1, 3])
    t$window_kb <- as.numeric(matches[1, 5])
    nmse_in_times_list[[row_index]] <- t
  }
  nmse_in_times_table <- dplyr::bind_rows(nmse_in_times_list)
  nmse_in_times_table$n_markers <- as.factor(nmse_in_times_table$n_markers)
  nmse_in_times_table$window_kb <- as.factor(nmse_in_times_table$window_kb)

  for (phenotype_model_index in c(0, 1)) {
    p <- ggplot2::ggplot(
      data = dplyr::filter(nmse_in_times_table, phenotype_model == phenotype_model_index),
      ggplot2::aes(x = epoch, y = nmse, color = window_kb, lty = n_markers)
    ) + ggplot2::geom_line() + ggplot2::scale_y_log10()

    testthat::expect_equal(1, length(phenotype_model_index))
    png_filename <- paste0("nmse_28_and_29_1_plot_p", phenotype_model_index,".png")
    message(png_filename)
    ggplot2::ggsave(
      filename = png_filename,
      plot = p, width = 7, height = 7)
    p

    q <- p + ggplot2::facet_grid(n_markers ~ window_kb)
    ggplot2::ggsave(
      filename = paste0("nmse_28_and_29_facet_grid_p", phenotype_model_index,".png"),
      plot = q, width = 7, height = 7
    )
    q
  }

}
################################################################################
# * 2. Genotype concordance
################################################################################
if (nchar("analysis_about_genotype_concordances")) {
  all_genotype_concordances_filenames <- c(
    list.files(folder_names, pattern = "genotype_concordances.csv", recursive = TRUE, full.names = TRUE)
  )
  # GCAE also stores a file called genotype_concordances
  genotype_concordances_filenames <- stringr::str_subset(
    all_genotype_concordances_filenames,
    pattern = "ae/geno"
  )
  genotype_concordances_filenames_table <- tibble::tibble(
    filename = genotype_concordances_filenames
  )
  genotype_concordancess_list <- list()
  for (row_index in seq_len(nrow(genotype_concordances_filenames_table))) {
    genotype_concordances_filename <- genotype_concordances_filenames_table$filename[row_index]

    t <- gcaer::read_genotype_concordances_file(
      genotype_concordances_filename = genotype_concordances_filename
    )
    matches <- stringr::str_match(
      genotype_concordances_filename,
      "issue_([[:digit:]]{2})_1000_epochs_p([[:digit:]]{1})/data_issue_([[:digit:]]{2})_([[:digit:]]{1,4})_ae"
    )
    testthat::expect_true(!is.na(matches)[1][1])
    testthat::expect_equal(matches[1, 2], matches[1, 4])
    issue <- matches[1, 2]
    t$n_markers <- issue_to_n_markers(as.numeric(issue))
    t$phenotype_model <- as.numeric(matches[1, 3])
    t$window_kb <- as.numeric(matches[1, 5])
    genotype_concordancess_list[[row_index]] <- t
  }
  genotype_concordancess_table <- dplyr::bind_rows(genotype_concordancess_list)
  genotype_concordancess_table$n_markers <- as.factor(genotype_concordancess_table$n_markers)
  genotype_concordancess_table$window_kb <- as.factor(genotype_concordancess_table$window_kb)

  for (phenotype_model_index in c(0, 1)) {
    p <- ggplot2::ggplot(
      data = dplyr::filter(genotype_concordancess_table, phenotype_model == phenotype_model_index),
      ggplot2::aes(x = epoch, y = genotype_concordance, color = window_kb, lty = n_markers)
    ) + ggplot2::geom_line()
    p
    ggplot2::ggsave(
      filename = paste0("genotype_concordance_28_and_29_1_plot_p", phenotype_model_index,".png"),
      plot = p, width = 7, height = 7
    )
    q <- p + ggplot2::facet_grid(n_markers ~ window_kb)
    q
    ggplot2::ggsave(
      filename = paste0("genotype_concordance_28_and_29_facet_grid_p", phenotype_model_index, ".png"),
      plot = q, width = 7, height = 7
    )
  }
}
################################################################################
# * 3. Runtime
################################################################################
if (nchar("analysis_about_runtime")) {
  all_log_filenames <- c(
    list.files(folder_names, pattern = ".log", recursive = TRUE, full.names = TRUE)
  )
  t <- tibble::tibble(
    filename = stringr::str_subset(all_log_filenames, "25_")
  )
  matches <- stringr::str_match(
    t$filename,
    "issue_([[:digit:]]{2})_([[:digit:]]{1,4})\\.log"
  )
  t$n_matches <- issue_to_n_markers(as.numeric(matches[, 2]))
  t$window_kb <- as.numeric(matches[, 3])
  t$runtime_sec <- NA

  for (row_index in seq_len(nrow(t))) {
    text <- readr::read_lines(t$filename[row_index])
    line <- stringr::str_subset(text, "Duration")
    testthat::expect_equal(1, length(line))
    matches <- stringr::str_match(line, "Duration: (.*) seconds")
    t$runtime_sec[row_index] <- as.numeric(matches[1, 2])
  }
  t$runtime_min <- t$runtime_sec / 60
  t$runtime_hours <- t$runtime_min / 60
  p <- ggplot2::ggplot(
    t,
    ggplot2::aes(x = window_kb, y = runtime_hours, color = n_matches)
  ) + ggplot2::geom_point(size = 10) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::geom_smooth(method = "lm", se = FALSE, lty = "dashed") +
    ggplot2::geom_smooth(
      method = "lm", formula = y ~ x + I(x^2),  se = FALSE, lty = "dashed"
    )
  p
  ggplot2::ggsave(filename = "runtime_hours_28_and_29_1_plot.png", plot = p, width = 7, height = 7)
}

################################################################################
# * 4. R sqaured
################################################################################
if (nchar("analysis_about_nsme")) {
  r_squared_in_time_filenames_table <- tibble::tibble(
    filename = c(
      list.files(folder_names, pattern = "r_squared_in_time.csv", recursive = TRUE, full.names = TRUE)
    )
  )
  r_squared_in_times_list <- list()
  for (row_index in seq_len(nrow(r_squared_in_time_filenames_table))) {
    r_squared_in_time_filename <- r_squared_in_time_filenames_table$filename[row_index]
    t <- gcaer::read_r_squared_in_time_file(
      r_squared_in_time_filename = r_squared_in_time_filename
    )
    matches <- stringr::str_match(
      r_squared_in_time_filename,
      "issue_([[:digit:]]{2})_1000_epochs_p([[:digit:]]{1})_m3d/data_issue_([[:digit:]]{2})_([[:digit:]]{1,4})_ae"
    )
    testthat::expect_true(!is.na(matches)[1][1])
    testthat::expect_equal(matches[1, 2], matches[1, 4])
    issue <- matches[1, 2]
    t$n_markers <- issue_to_n_markers(as.numeric(issue))
    t$phenotype_model <- as.numeric(matches[1, 3])
    t$window_kb <- as.numeric(matches[1, 5])
    r_squared_in_times_list[[row_index]] <- t
  }
  r_squared_in_times_table <- dplyr::bind_rows(r_squared_in_times_list)
  r_squared_in_times_table$n_markers <- as.factor(r_squared_in_times_table$n_markers)
  r_squared_in_times_table$window_kb <- as.factor(r_squared_in_times_table$window_kb)

  p <- ggplot2::ggplot(
    data = r_squared_in_times_table,
    ggplot2::aes(x = epoch, y = r_squared, color = window_kb, lty = n_markers)
  ) + ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      name = "r_squared",
      limits = c(0, 1)
    )
  p
  testthat::expect_equal(1, length(phenotype_model_index))
  png_filename <- paste0("r_squared_28_and_29_1_plot_p1.png")
  message(png_filename)
  ggplot2::ggsave(
    filename = png_filename,
    plot = p, width = 7, height = 7)
  p

  q <- p + ggplot2::facet_grid(n_markers ~ window_kb)
  ggplot2::ggsave(
    filename = paste0("r_squared_28_and_29_facet_grid_p1.png"),
    plot = q, width = 7, height = 7
  )
  q
}
