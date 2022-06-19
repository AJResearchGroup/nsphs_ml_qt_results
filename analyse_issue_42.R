# Analysis of R squareds in time
################################################################################
# * 0. Functions needed for the analysis
################################################################################
folder_name <- "issue_42_1000_epochs_p1_m3d"
testthat::expect_true(dir.exists(folder_name))
r_squared_filenames <- list.files(path = folder_name, pattern = "r_squared_in_time.csv", recursive = TRUE, full.names = TRUE)
testthat::expect_true(length(r_squared_filenames) > 0)


t_matches <- stringr::str_match(
  r_squared_filenames,
  ".*data_issue_42_(M.*)_(p[[:digit:]])_([[:digit:]]{1,4})_ae.*"
)
t_matches <- tibble::as_tibble(t_matches)
names(t_matches) <- c("filename", "model_id", "pheno_model_id", "window_kb")
t_matches
t_matches$last_r_squared <- NA

for (i in seq_len(nrow(t_matches))) {
  t_matches$last_r_squared[i] <- as.numeric(
    utils::tail(
      gcaer::read_r_squared_in_time_file(
        t_matches$filename[i]
      )$r_squared,
      n = 1
    )
  )
}


p <- ggplot2::ggplot(
  data = t_matches,
  ggplot2::aes(x = window_kb, y = last_r_squared)
) + ggplot2::geom_col() +
  ggplot2::scale_y_continuous(
    limits = c(0, 1)
  ) + ggplot2::facet_grid(model_id ~ pheno_model_id)
p
ggplot2::ggsave(filename = "issue_42.png", plot = p, width = 7, height = 7)
