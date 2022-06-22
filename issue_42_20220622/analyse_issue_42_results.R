#############################################################################
# General setup
#############################################################################
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
# Helper functions
#############################################################################
window_kb_to_data <- function(window_kb) {
  if (window_kb == 1) return("1 kb, 2 SNPs")
  if (window_kb == 10) return("10 kb, 5 SNPs")
  if (window_kb == 100) return("100 kb, 5 SNPs")
  if (window_kb == 1000) return("1000 kb, 6 SNPs")
  stop("Unknown 'window_kb': ", window_kb)
}

testthat::expect_true(stringr::str_detect(window_kb_to_data(1), "2 SNPs"))
testthat::expect_true(stringr::str_detect(window_kb_to_data(10), "5 SNPs"))
testthat::expect_true(stringr::str_detect(window_kb_to_data(100), "5 SNPs"))
testthat::expect_true(stringr::str_detect(window_kb_to_data(1000), "6 SNPs"))

#############################################################################
# Get all combinations
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

#############################################################################
# NMSE in time
#############################################################################
nmse_in_time_list <- list()
for (i in seq_len(nrow(t_combinations))) {
  this_row <- t_combinations[i, ]
  filename <- file.path(
    this_row$folder_name,
    "nmse_in_time.csv"
  )
  if (!file.exists(filename)) {
    message(filename, " does not exist")
    next
  }
  t_nmses_in_time <- gcaer::read_nmse_in_time_file(
    filename
  )
  t_nmses_in_time$model_id <- this_row$model_id
  t_nmses_in_time$pheno_model_id <- this_row$pheno_model_id
  t_nmses_in_time$window_kb <- this_row$window_kb
  nmse_in_time_list[[i]] <- t_nmses_in_time
}
t_nmses_in_time <- dplyr::bind_rows(nmse_in_time_list)

t_nmses_in_time$genetic_data <- Vectorize(window_kb_to_data)(t_nmses_in_time$window_kb)

ggplot2::ggplot(
  t_nmses_in_time,
  ggplot2::aes(x = epoch, y = nmse, color = genetic_data)
) + ggplot2::geom_line(size = 1) +
  ggplot2::scale_y_continuous(limits = c(0, 2), oob = scales::squish) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  ggplot2::facet_grid(model_id ~ pheno_model_id) +
  gcaer::get_gcaer_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
  )
ggplot2::ggsave(filename = "issue_42_nmses_in_time.png", width = 7, height = 7)


#############################################################################
# r_squared in time
#############################################################################
r_squared_in_time_list <- list()
for (i in seq_len(nrow(t_combinations))) {
  this_row <- t_combinations[i, ]
  filename <- file.path(
    this_row$folder_name,
    "r_squared_in_time.csv"
  )
  if (!file.exists(filename)) {
    message(filename, " does not exist")
    next
  }
  t_r_squareds_in_time <- gcaer::read_r_squared_in_time_file(
    filename
  )
  t_r_squareds_in_time$model_id <- this_row$model_id
  t_r_squareds_in_time$pheno_model_id <- this_row$pheno_model_id
  t_r_squareds_in_time$window_kb <- this_row$window_kb
  r_squared_in_time_list[[i]] <- t_r_squareds_in_time
}
t_r_squareds_in_time <- dplyr::bind_rows(r_squared_in_time_list)

t_r_squareds_in_time$genetic_data <- Vectorize(window_kb_to_data)(t_r_squareds_in_time$window_kb)

ggplot2::ggplot(
  t_r_squareds_in_time,
  ggplot2::aes(x = epoch, y = r_squared, color = genetic_data)
) + ggplot2::geom_line(size = 1) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  ggplot2::facet_grid(model_id ~ pheno_model_id) +
  gcaer::get_gcaer_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
  )
ggplot2::ggsave(filename = "issue_42_r_squareds_in_time.png", width = 7, height = 7)



#############################################################################
# last genotype_concordance in time
#############################################################################
all_epochs <- sort(unique(t_genotype_concordances_in_time$epoch))
top_epochs <- all_epochs[all_epochs > max(all_epochs) * 0.89]

t_last_genotype_concordances_in_time <- dplyr::filter(
  t_genotype_concordances_in_time,
  epoch %in% top_epochs
)
t_last_genotype_concordances_in_time$genetic_data <- as.factor(t_last_genotype_concordances_in_time$genetic_data)
ggplot2::ggplot(
  t_last_genotype_concordances_in_time,
  ggplot2::aes(x = genetic_data, y = genotype_concordance, fill = genetic_data, color = genetic_data)
) + ggplot2::geom_boxplot() +
  ggplot2::scale_x_discrete(name = "", drop = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::facet_grid(model_id ~ pheno_model_id) +
  gcaer::get_gcaer_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
  )

ggplot2::ggsave(filename = "issue_42_genotype_concordances_at_end.png", width = 7, height = 7)

#############################################################################
# last nmse in time
#############################################################################
all_epochs <- sort(unique(t_nmses_in_time$epoch))
top_epochs <- all_epochs[all_epochs > max(all_epochs) * 0.89]

t_last_nmses_in_time <- dplyr::filter(
  t_nmses_in_time,
  epoch %in% top_epochs
)
t_last_nmses_in_time$genetic_data <- as.factor(t_last_nmses_in_time$genetic_data)
ggplot2::ggplot(
  t_last_nmses_in_time,
  ggplot2::aes(x = 1, y = nmse, fill = genetic_data)
) + ggplot2::geom_boxplot() +
  ggplot2::scale_x_discrete(name = "", drop = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 2), oob = scales::squish) +
  ggplot2::facet_grid(model_id ~ pheno_model_id) +
  gcaer::get_gcaer_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
  )

ggplot2::ggsave(filename = "issue_42_nmses_at_end.png", width = 7, height = 7)

#############################################################################
# last r_squared in time
#############################################################################
all_epochs <- sort(unique(t_r_squareds_in_time$epoch))
top_epochs <- all_epochs[all_epochs > max(all_epochs) * 0.89]

t_last_r_squareds_in_time <- dplyr::filter(
  t_r_squareds_in_time,
  epoch %in% top_epochs
)

if ("calculate averages" == "myself") {
  t_average_last_r_squareds_in_time <- dplyr::summarise(
    dplyr::group_by(
      t_last_r_squareds_in_time,
      model_id, pheno_model_id, genetic_data
    ),
    average_last_r_squared = mean(r_squared),
    .groups = "drop"
  )
}

t_last_r_squareds_in_time$genetic_data <- as.factor(t_last_r_squareds_in_time$genetic_data)
ggplot2::ggplot(
  t_last_r_squareds_in_time,
  ggplot2::aes(x = 1, y = r_squared, fill = genetic_data)
) + ggplot2::geom_boxplot() +
  ggplot2::scale_x_discrete(name = "", drop = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::facet_grid(model_id ~ pheno_model_id) +
  gcaer::get_gcaer_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
  )

ggplot2::ggsave(filename = "issue_42_r_squareds_at_end.png", width = 7, height = 7)
