# Run training of GCAE example on Rackham #78
gcae_options <- gcaer::create_gcae_options(
  gcae_folder = tempfile()
)
datadir <- file.path(
  gcaer::get_gcae_subfolder(gcae_options = gcae_options),
  "example_tiny/"
)
data <- "HumanOrigins249_tiny"
gcae_setup <- gcaer::create_gcae_setup()
epochs <- 1000
save_interval <- 1

testthat::expect_false(gcaer::is_gcae_installed(gcae_options = gcae_options))
gcaer::install_gcae(gcae_options = gcae_options)
testthat::expect_true(gcaer::is_gcae_installed(gcae_options = gcae_options))

train_filenames <- gcaer::gcae_train(
  gcae_options = gcae_options,
  datadir = datadir,
  data = data,
  gcae_setup = gcae_setup,
  epochs = epochs,
  save_interval = save_interval
)

train_results <- gcaer::parse_train_filenames(
  train_filenames
)

readr::write_csv(
  train_results$train_times_table,
  "issue_78_train_times.csv"
)
ggplot2::ggplot(
  train_results$train_times_table,
  ggplot2::aes(x = epoch, y = train_times_sec)
) +
  ggplot2::geom_line() +
  ggplot2::scale_y_continuous(limits = c(0.0, NA))
ggplot2::ggsave("issue_78_train_times.png")

readr::write_csv(
  train_results$losses_from_train_t_table,
  "issue_78_losses_from_train_t.csv"
)
ggplot2::ggplot(
  train_results$losses_from_train_t_table,
  ggplot2::aes(x = epoch, y = losses_from_train_t)
) + ggplot2::geom_line() +
  ggplot2::scale_y_continuous(limits = c(0.0, 1.0))
ggplot2::ggsave("issue_78_losses_from_train_t.png")

readr::write_csv(
  train_results$losses_from_train_v_table,
  "issue_78_losses_from_train_v.csv"
)
ggplot2::ggplot(
  train_results$losses_from_train_v_table,
  ggplot2::aes(x = epoch, y = losses_from_train_v)
) + ggplot2::geom_line() +
  ggplot2::scale_y_continuous(limits = c(0.0, 1.0))
ggplot2::ggsave("issue_78_losses_from_train_v.png")

unlink(gcae_options$gcae_folder, recursive = TRUE)


