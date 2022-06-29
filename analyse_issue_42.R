#############################################################################
# General setup
#############################################################################
path <- "~/GitHubs/nsphs_ml_qt_results/issue_42_20220622/"

all_root_folders <- list.dirs(
  path = path,
  recursive = FALSE
)
all_folders <- list.dirs(
  path = path,
  recursive = TRUE
)

successful_root_folders <- sort(unique(dirname(dirname(stringr::str_subset(all_folders, "pheno_weights")))))
testthat::expect_equal(32, length(successful_root_folders))

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
  n_associations <- nsphsmlqt::window_kb_to_n_associations(
    protein_name = "IL-17RA",
    window_kb = window_kb
  )
  paste0(window_kb, " kb, ", n_associations, " associations")
}

testthat::expect_true(stringr::str_detect(window_kb_to_data(1), "2 associations")) #nolint indeed a long line
testthat::expect_true(stringr::str_detect(window_kb_to_data(10), "4 associations")) #nolint indeed a long line
testthat::expect_true(stringr::str_detect(window_kb_to_data(100), "5 associations")) #nolint indeed a long line
testthat::expect_true(stringr::str_detect(window_kb_to_data(1000), "6 associations")) #nolint indeed a long line

window_kbs_to_data <- function(window_kbs) {
  # LUT = Look Up Table
  t_lut <- tibble::tibble(
    window_kb = unique(window_kbs)
  )
  t_lut$n_associations <- Vectorize(window_kb_to_data)(t_lut$window_kb)
  t <- tibble::tibble(window_kb = window_kbs)
  t_merged <- dplyr::left_join(t, t_lut, by = "window_kb")
  testthat::expect_true(identical(t$window_kb, t_merged$window_kb))
  testthat::expect_equal(length(window_kbs), length(t_merged$n_associations))
  t_merged$n_associations
}

one_kb <- window_kb_to_data(1)
ten_kb <- window_kb_to_data(10)
hundred_kb <- window_kb_to_data(100)
thousand_kb <- window_kb_to_data(1000)
window_kbs <- rep(c(1, 10, 100, 1000), each = 3, times = 2)
data <- rep(c(one_kb, ten_kb, hundred_kb, thousand_kb), each = 3, times = 2)
testthat::expect_equal(
  window_kbs_to_data(window_kbs),
  data
)

#############################################################################
# Get all combinations
#############################################################################
t_combinations <- tidyr::expand_grid(
  model_id = model_ids_that_work,
  pheno_model_id = pheno_model_ids_that_work,
  window_kb = window_kbs_that_work
)
t_combinations$window_kb <- as.numeric(t_combinations$window_kb)
t_combinations$folder_name <- file.path(
  path,
  paste0(
    "data_issue_42_", t_combinations$model_id,
    "_", t_combinations$pheno_model_id, "_", t_combinations$window_kb,
    "_ae"
  )
)

t_combinations$log_filename_25 <- file.path(
  path,
  paste0(
    "25_run_issue_42_", t_combinations$model_id,
    "_", t_combinations$pheno_model_id, "_", t_combinations$window_kb,
    ".log"
  )
)
dir.exists("~/GitHubs/nsphs_ml_qt_results/issue_42_20220622//data_issue_42_M1_p1_100_ae")
t_combinations$folder_exists <- dir.exists(t_combinations$folder_name)
t_combinations$log_filename_25_exists <- file.exists(t_combinations$log_filename_25)

# To re-run:
t_combinations[!t_combinations$folder_exists | !t_combinations$log_filename_25_exists, ]

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
  testthat::expect_true(is.numeric(this_row$window_kb))
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

t_genotype_concordances_in_time$genetic_data <- window_kbs_to_data(t_genotype_concordances_in_time$window_kb)

dplyr::distinct(dplyr::select(t_genotype_concordances_in_time, window_kb, genetic_data))
#dplyr::select(t_genotype_concordances_in_time, model_id, pheno_model_id, )

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

t_nmses_in_time$genetic_data <- window_kbs_to_data(t_nmses_in_time$window_kb)

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

t_r_squareds_in_time$genetic_data <- window_kbs_to_data(t_r_squareds_in_time$window_kb)

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
  ggplot2::aes(x = genetic_data, y = nmse, fill = genetic_data)
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

t_last_r_squareds_in_time$genetic_data <- as.factor(t_last_r_squareds_in_time$genetic_data)
ggplot2::ggplot(
  t_last_r_squareds_in_time,
  ggplot2::aes(x = genetic_data, y = r_squared, fill = genetic_data)
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

#############################################################################
# get the top averages for each
#############################################################################
t_average_last_genotype_concordance_in_time <- dplyr::summarise(
  dplyr::group_by(
    t_last_genotype_concordances_in_time,
    model_id, pheno_model_id, genetic_data
  ),
  average_last_genotype_concordance = mean(genotype_concordance),
  .groups = "drop"
)
top_10_gc <- head(
  t_average_last_genotype_concordance_in_time[order(t_average_last_genotype_concordance_in_time$average_last_genotype_concordance, decreasing = TRUE),  ],
  n = 10
)

top_10_gc_tex <- print(
  xtable::xtable(
    top_10_gc,
    label = "tab:top_10_gc",
    caption = "Top 10 models with the highest average genotype concordance"
  ),
  include.colnames = TRUE
)
writeLines(top_10_gc_tex, "~/issue_42_top_10_gc.tex")

t_average_last_nmse_in_time <- dplyr::summarise(
  dplyr::group_by(
    t_last_nmses_in_time,
    model_id, pheno_model_id, genetic_data
  ),
  average_last_nmse = mean(nmse),
  .groups = "drop"
)

top_10_nmse <- head(
  t_average_last_nmse_in_time[order(t_average_last_nmse_in_time$average_last_nmse, decreasing = FALSE),  ],
  n = 10
)

top_10_nmse_tex <- print(
  xtable::xtable(
    top_10_nmse,
    label = "tab:top_10_nmse",
    caption = "Top 10 models with the lowest normalized mean squared error"
  ),
  include.colnames = TRUE
)
writeLines(top_10_nmse_tex, "~/issue_42_top_10_nmse.tex")


t_average_last_r_squareds_in_time <- dplyr::summarise(
  dplyr::group_by(
    t_last_r_squareds_in_time,
    model_id, pheno_model_id, genetic_data
  ),
  average_last_r_squared = mean(r_squared),
  .groups = "drop"
)

top_10_rs <- head(
  t_average_last_r_squareds_in_time[order(t_average_last_r_squareds_in_time$average_last_r_squared, decreasing = TRUE),  ],
  n = 10
)

top_10_rs_tex <- print(
  xtable::xtable(
    top_10_rs,
    label = "tab:top_10_rs",
    caption = "Top 10 models with the heighest r squared"
  ),
  include.colnames = TRUE
)
writeLines(top_10_rs_tex, "~/issue_42_top_10_rs.tex")
