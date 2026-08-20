#' Create Simulated Reference and Focal Data
#'
#' Generates simulated combined reference/focal datasets for predictive bias
#' research. Varying simulation parameters produce separate CSV files; fixed
#' parameters must each be supplied as a single value.
#'
#' @param output_location Character string path to the parent directory where
#'   dated output folders are created.
#' @param reference_m_x Reference group predictor mean. Single value only.
#' @param reference_sd_x Reference group predictor SD. Single value only.
#' @param reference_sd_y Reference group criterion SD. Single value only.
#' @param focal_m_x Focal group predictor mean. Single value only.
#' @param focal_sd_x Focal group predictor SD. Single value only.
#' @param iteration Number of iterations per parameter combination. Single
#'   value only.
#' @param x_var_name Name for the predictor column in output files. Single
#'   value only.
#' @param y_var_name Name for the criterion column in output files. Single
#'   value only.
#' @param reference_validity Reference group validity (correlation). One or
#'   more values.
#' @param focal_validity Focal group validity (correlation). One or more values.
#' @param focal_sd_y Focal group criterion SD. One or more values.
#' @param intercept_diffs Focal group intercept difference. One or more values.
#' @param D Predictor mean difference between groups. One or more values.
#' @param n Total sample size. One or more values.
#' @param population_percentage Proportion of the sample in the focal group.
#'   One or more values.
#' @param seed Optional integer seed for reproducible simulation.
#'
#' @return A list with elements `output_folder` (path to the created folder),
#'   `key_file` (path to the parameter key text file), `files` (vector of CSV
#'   paths), and `conditions` (tibble of parameter combinations).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tmp_dir <- tempdir()
#' results <- create_sim_dat(
#'   output_location = tmp_dir,
#'   reference_validity = c(0.1, 0.2),
#'   focal_validity = 0.1,
#'   focal_sd_y = 1,
#'   intercept_diffs = 0.2,
#'   D = 0.7,
#'   n = 500,
#'   population_percentage = 0.5,
#'   iteration = 1
#' )
#' results$output_folder
#' }
create_sim_dat <- function(
    output_location,
    reference_m_x = 0,
    reference_sd_x = 1,
    reference_sd_y = 1,
    focal_m_x = 0,
    focal_sd_x = 1,
    iteration = 1,
    x_var_name = "x",
    y_var_name = "y",
    reference_validity = 0.1,
    focal_validity = 0.1,
    focal_sd_y = 1,
    intercept_diffs = 0,
    D = 0.7,
    n = 5000,
    population_percentage = 0.5,
    seed = NULL) {
  single_value_params <- list(
    output_location = output_location,
    reference_m_x = reference_m_x,
    reference_sd_x = reference_sd_x,
    reference_sd_y = reference_sd_y,
    focal_m_x = focal_m_x,
    focal_sd_x = focal_sd_x,
    iteration = iteration,
    x_var_name = x_var_name,
    y_var_name = y_var_name
  )

  for (param_name in names(single_value_params)) {
    param_value <- single_value_params[[param_name]]
    if (length(param_value) != 1) {
      stop(
        param_name,
        " must be a single value, but ",
        length(param_value),
        " values were supplied.",
        call. = FALSE
      )
    }
  }

  if (!is.character(output_location) || is.na(output_location) || output_location == "") {
    stop("output_location must be a non-empty character string.", call. = FALSE)
  }

  if (!is.character(x_var_name) || !is.character(y_var_name)) {
    stop("x_var_name and y_var_name must be character strings.", call. = FALSE)
  }

  if (x_var_name == y_var_name) {
    stop("x_var_name and y_var_name must be different.", call. = FALSE)
  }

  multi_value_params <- list(
    reference_validity = reference_validity,
    focal_validity = focal_validity,
    focal_sd_y = focal_sd_y,
    intercept_diffs = intercept_diffs,
    D = D,
    n = n,
    population_percentage = population_percentage
  )

  for (param_name in names(multi_value_params)) {
    param_value <- multi_value_params[[param_name]]
    if (length(param_value) < 1) {
      stop(param_name, " must contain at least one value.", call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  conditions <- tidyr::expand_grid(
    reference_m_x = reference_m_x,
    reference_sd_x = reference_sd_x,
    reference_sd_y = reference_sd_y,
    focal_m_x = focal_m_x,
    focal_sd_x = focal_sd_x,
    x_var_name = x_var_name,
    y_var_name = y_var_name,
    reference_validity = reference_validity,
    focal_validity = focal_validity,
    focal_sd_y = focal_sd_y,
    intercept_diffs = intercept_diffs,
    D = D,
    n = n,
    population_percentage = population_percentage,
    iteration = seq_len(iteration)
  ) %>%
    dplyr::mutate(condition = dplyr::row_number())

  output_folder <- .create_dated_output_folder(output_location)
  file_paths <- character(nrow(conditions))
  key_lines <- c(
    "Simulation Data Parameter Key",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Output folder: ", basename(output_folder)),
    ""
  )

  for (i in seq_len(nrow(conditions))) {
    condition_row <- conditions[i, ]
    focal_n <- round(condition_row$population_percentage * condition_row$n, digits = 0)
    reference_n <- condition_row$n - focal_n

    combined_dat <- .generate_combined_sim_dat(
      reference_n = reference_n,
      focal_n = focal_n,
      reference_m_x = condition_row$reference_m_x,
      reference_sd_x = condition_row$reference_sd_x,
      reference_sd_y = condition_row$reference_sd_y,
      reference_validity = condition_row$reference_validity,
      focal_m_x = condition_row$focal_m_x,
      focal_sd_x = condition_row$focal_sd_x,
      focal_sd_y = condition_row$focal_sd_y,
      focal_validity = condition_row$focal_validity,
      intercept_diffs = condition_row$intercept_diffs,
      D = condition_row$D,
      condition = condition_row$condition,
      iteration = condition_row$iteration,
      x_var_name = condition_row$x_var_name,
      y_var_name = condition_row$y_var_name
    )

    file_name <- paste0("simdat", i, ".csv")
    file_path <- file.path(output_folder, file_name)
    readr::write_csv(combined_dat, file_path)
    file_paths[i] <- file_path

    key_lines <- c(
      key_lines,
      file_name,
      .format_condition_key(condition_row),
      ""
    )
  }

  key_file <- file.path(output_folder, "parameter_key.txt")
  writeLines(key_lines, con = key_file)

  list(
    output_folder = output_folder,
    key_file = key_file,
    files = file_paths,
    conditions = conditions
  )
}

.create_dated_output_folder <- function(output_location) {
  if (!dir.exists(output_location)) {
    dir.create(output_location, recursive = TRUE)
  }

  date_cleaned <- gsub("-", "_", as.character(Sys.Date()))
  existing_dirs <- list.dirs(output_location, full.names = FALSE, recursive = FALSE)
  pattern <- paste0("^", date_cleaned, "_[0-9]+$")
  matching_dirs <- existing_dirs[grepl(pattern, existing_dirs)]

  if (length(matching_dirs) == 0) {
    folder_name <- paste0(date_cleaned, "_1")
  } else {
    suffixes <- as.integer(sub(paste0("^", date_cleaned, "_"), "", matching_dirs))
    folder_name <- paste0(date_cleaned, "_", max(suffixes, na.rm = TRUE) + 1)
  }

  folder_path <- file.path(output_location, folder_name)
  dir.create(folder_path, recursive = TRUE)
  folder_path
}

.generate_combined_sim_dat <- function(
    reference_n,
    focal_n,
    reference_m_x,
    reference_sd_x,
    reference_sd_y,
    reference_validity,
    focal_m_x,
    focal_sd_x,
    focal_sd_y,
    focal_validity,
    intercept_diffs,
    D,
    condition,
    iteration,
    x_var_name,
    y_var_name) {
  combined_dat <-
    tibble::tibble(
      ref_foc = "a_reference",
      e = stats::rnorm(reference_n, 0, 1),
      x = stats::rnorm(reference_n, reference_m_x, reference_sd_x)
    ) %>%
    dplyr::mutate(
      y = reference_sd_y * (0 + reference_validity * x + sqrt(1 - reference_validity^2) * e)
    ) %>%
    dplyr::bind_rows(
      tibble::tibble(
        ref_foc = "b_focal",
        e = stats::rnorm(focal_n, 0, 1),
        x = stats::rnorm(focal_n, focal_m_x, focal_sd_x)
      ) %>%
        dplyr::mutate(
          x = x - D,
          y = focal_sd_y * (intercept_diffs + focal_validity * x + sqrt(1 - focal_validity^2) * e)
        )
    ) %>%
    dplyr::mutate(
      condition = condition,
      iteration = iteration,
      dummy_column = 1
    )

  combined_dat %>%
    dplyr::select(ref_foc, x, y, condition, iteration) %>%
    dplyr::rename(
      !!x_var_name := x,
      !!y_var_name := y
    )
}

.format_condition_key <- function(condition_row) {
  param_names <- c(
    "reference_m_x",
    "reference_sd_x",
    "reference_sd_y",
    "focal_m_x",
    "focal_sd_x",
    "x_var_name",
    "y_var_name",
    "reference_validity",
    "focal_validity",
    "focal_sd_y",
    "intercept_diffs",
    "D",
    "n",
    "population_percentage",
    "iteration",
    "condition"
  )

  paste0(
    "  ",
    param_names,
    ": ",
    vapply(
      param_names,
      function(param_name) {
        as.character(condition_row[[param_name]])
      },
      character(1)
    ),
    collapse = "\n"
  )
}
