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
    t$issue <- as.numeric(matches[1, 2])
    t$window_kb <- as.numeric(matches[1, 3])
    nmse_in_times_list[[row_index]] <- t
  }
  nmse_in_times_table <- dplyr::bind_rows(nmse_in_times_list)
  nmse_in_times_table$issue <- as.factor(nmse_in_times_table$issue)
  nmse_in_times_table$window_kb <- as.factor(nmse_in_times_table$window_kb)

  p <- ggplot2::ggplot(
    data = nmse_in_times_table,
    ggplot2::aes(x = epoch, y = nmse, color = window_kb, lty = issue)
  ) + ggplot2::geom_line() + ggplot2::scale_y_log10()
  ggplot2::ggsave(filename = "nmse_28_and_29_1_plot.png", plot = p, width = 7, height = 7)

  q <- p + ggplot2::facet_grid(issue ~ window_kb)
  ggplot2::ggsave(filename = "nmse_28_and_29_facet_grid.png", plot = q, width = 7, height = 7)

}

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
    t$issue <- as.numeric(matches[1, 2])
    t$window_kb <- as.numeric(matches[1, 3])
    genotype_concordancess_list[[row_index]] <- t
  }
  genotype_concordancess_table <- dplyr::bind_rows(genotype_concordancess_list)
  genotype_concordancess_table$issue <- as.factor(genotype_concordancess_table$issue)
  genotype_concordancess_table$window_kb <- as.factor(genotype_concordancess_table$window_kb)

  p <- ggplot2::ggplot(
    data = genotype_concordancess_table,
    ggplot2::aes(x = epoch, y = genotype_concordance, color = window_kb, lty = issue)
  ) + ggplot2::geom_line()
  p
  ggplot2::ggsave(filename = "genotype_concordance_28_and_29_1_plot.png", plot = p, width = 7, height = 7)
  q
  q <- p + ggplot2::facet_grid(issue ~ window_kb)
  ggplot2::ggsave(filename = "genotype_concordance_28_and_29_facet_grid.png", plot = q, width = 7, height = 7)

}

