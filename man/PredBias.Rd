% Generated manually. Edit R/PredBias.R and run devtools::document() to update.
\name{PredBias}
\alias{PredBias}
\title{Predictive Bias Analysis}
\usage{
PredBias(
  Combined,
  x,
  y,
  output_file_name,
  scatterplot_file_name,
  predictor_histogram_file_name,
  criterion_histogram_file_name
)
}
\arguments{
\item{Combined}{A data frame containing predictor, criterion, and group
columns. The group column must be named \code{ref_foc} with values
\code{"a_reference"} and \code{"b_focal"}.}

\item{x}{Character string naming the predictor column in \code{Combined}.}

\item{y}{Character string naming the criterion column in \code{Combined}.}

\item{output_file_name}{Path for the CSV file of analysis results.}

\item{scatterplot_file_name}{Path for the scatterplot PNG file.}

\item{predictor_histogram_file_name}{Path for the predictor histogram PNG.}

\item{criterion_histogram_file_name}{Path for the criterion histogram PNG.}
}
\value{
A tibble containing all computed metrics (also written to
\code{output_file_name}).
}
\description{
Runs a predictive bias analysis on a combined reference and focal group
dataset. Computes Gulliksen-Wilks tests, Cleary regression tests,
separate-group regressions, Fisher's z comparisons, residual metrics,
Levene tests, homoscedasticity diagnostics, and descriptive statistics.
Writes a CSV of results and optional PNG plots.
}
\examples{
\dontrun{
data_path <- system.file("extdata", "SampleData.csv", package = "iopsychology")
sample_data <- readr::read_csv(data_path, show_col_types = FALSE)

PredBias(
  Combined = sample_data,
  x = "test",
  y = "criterion",
  output_file_name = "SampleOutput.csv",
  scatterplot_file_name = "SampleScatterPlot.png",
  predictor_histogram_file_name = "SamplePredictorHistogram.png",
  criterion_histogram_file_name = "SampleCriterionHistogram.png"
)
}
}
\seealso{
The original script is preserved in the repository history as
\code{PredBias Function.R}.
}
