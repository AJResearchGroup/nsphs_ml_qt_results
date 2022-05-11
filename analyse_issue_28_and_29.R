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

################################################################################
# * 1. NMSE
################################################################################
if (nchar("analysis_about_nsme")) {
  nmse_in_time_filenames_table <- tibble::tibble(
    filename = c(
      list.files("issue_28", pattern = "nmse_in_time.csv", recursive = TRUE, full.names = TRUE),
      list.files("issue_29", pattern = "nmse_in_time.csv", recursive = TRUE, full.names = TRUE)
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
      "data_issue_([[:digit:]]{2})_([[:digit:]]{1,4})_ae"
    )
    t$n_markers <- issue_to_n_markers(as.numeric(matches[1, 2]))
    t$window_kb <- as.numeric(matches[1, 3])
    nmse_in_times_list[[row_index]] <- t
  }
  nmse_in_times_table <- dplyr::bind_rows(nmse_in_times_list)
  nmse_in_times_table$n_markers <- as.factor(nmse_in_times_table$n_markers)
  nmse_in_times_table$window_kb <- as.factor(nmse_in_times_table$window_kb)

  p <- ggplot2::ggplot(
    data = nmse_in_times_table,
    ggplot2::aes(x = epoch, y = nmse, color = window_kb, lty = n_markers)
  ) + ggplot2::geom_line() + ggplot2::scale_y_log10()
  ggplot2::ggsave(filename = "nmse_28_and_29_1_plot.png", plot = p, width = 7, height = 7)

  q <- p + ggplot2::facet_grid(n_markers ~ window_kb)
  ggplot2::ggsave(filename = "nmse_28_and_29_facet_grid.png", plot = q, width = 7, height = 7)

}
################################################################################
# * 2. Genotype concordance
################################################################################
if (nchar("analysis_about_genotype_concordances")) {
  all_genotype_concordances_filenames <- c(
    list.files("issue_28", pattern = "genotype_concordances.csv", recursive = TRUE, full.names = TRUE),
    list.files("issue_29", pattern = "genotype_concordances.csv", recursive = TRUE, full.names = TRUE)
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
      "data_issue_([[:digit:]]{2})_([[:digit:]]{1,4})_ae"
    )
    t$n_markers <- issue_to_n_markers(as.numeric(matches[1, 2]))
    t$window_kb <- as.numeric(matches[1, 3])
    genotype_concordancess_list[[row_index]] <- t
  }
  genotype_concordancess_table <- dplyr::bind_rows(genotype_concordancess_list)
  genotype_concordancess_table$n_markers <- as.factor(genotype_concordancess_table$n_markers)
  genotype_concordancess_table$window_kb <- as.factor(genotype_concordancess_table$window_kb)

  p <- ggplot2::ggplot(
    data = genotype_concordancess_table,
    ggplot2::aes(x = epoch, y = genotype_concordance, color = window_kb, lty = n_markers)
  ) + ggplot2::geom_line()
  p
  ggplot2::ggsave(filename = "genotype_concordance_28_and_29_1_plot.png", plot = p, width = 7, height = 7)
  q
  q <- p + ggplot2::facet_grid(n_markers ~ window_kb)
  ggplot2::ggsave(filename = "genotype_concordance_28_and_29_facet_grid.png", plot = q, width = 7, height = 7)

}
################################################################################
# * 3. Runtime
################################################################################
if (nchar("analysis_about_runtime")) {
  all_log_filenames <- c(
    list.files("issue_28", pattern = ".log", recursive = TRUE, full.names = TRUE),
    list.files("issue_29", pattern = ".log", recursive = TRUE, full.names = TRUE)
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
    ggplot2::geom_smooth(method = "lm", formula = y ~ x + I(x^2),  se = FALSE, lty = "dashed")
    p
  ggplot2::ggsave(filename = "runtime_hours_28_and_29_1_plot.png", plot = p, width = 7, height = 7)
}
