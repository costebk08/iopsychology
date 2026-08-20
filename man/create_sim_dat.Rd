% Generated manually. Edit R/create_sim_dat.R and run devtools::document() to update.
\name{create_sim_dat}
\alias{create_sim_dat}
\title{Create Simulated Reference and Focal Data}
\usage{
create_sim_dat(
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
  seed = NULL
)
}
\arguments{
\item{output_location}{Character string path to the parent directory where
dated output folders are created.}

\item{reference_m_x}{Reference group predictor mean. Single value only.}

\item{reference_sd_x}{Reference group predictor SD. Single value only.}

\item{reference_sd_y}{Reference group criterion SD. Single value only.}

\item{focal_m_x}{Focal group predictor mean. Single value only.}

\item{focal_sd_x}{Focal group predictor SD. Single value only.}

\item{iteration}{Number of iterations per parameter combination. Single
value only.}

\item{x_var_name}{Name for the predictor column in output files. Single
value only.}

\item{y_var_name}{Name for the criterion column in output files. Single
value only.}

\item{reference_validity}{Reference group validity (correlation). One or
more values.}

\item{focal_validity}{Focal group validity (correlation). One or more values.}

\item{focal_sd_y}{Focal group criterion SD. One or more values.}

\item{intercept_diffs}{Focal group intercept difference. One or more values.}

\item{D}{Predictor mean difference between groups. One or more values.}

\item{n}{Total sample size. One or more values.}

\item{population_percentage}{Proportion of the sample in the focal group.
One or more values.}

\item{seed}{Optional integer seed for reproducible simulation.}
}
\value{
A list with elements \code{output_folder} (path to the created folder),
\code{key_file} (path to the parameter key text file), \code{files} (vector of CSV
paths), and \code{conditions} (tibble of parameter combinations).
}
\description{
Generates simulated combined reference/focal datasets for predictive bias
research. Varying simulation parameters produce separate CSV files; fixed
parameters must each be supplied as a single value.
}
