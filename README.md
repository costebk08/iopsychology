# iopsychology

Industrial-organizational psychology tools for R.

## Installation

Install from GitHub with [devtools](https://devtools.r-lib.org/):

```r
# install.packages("devtools")
devtools::install_github("costebk08/iopsychology")
```

## Quick start

The package includes sample data in `inst/extdata/SampleData.csv`. Your input
data must include a `ref_foc` column with values `"a_reference"` and
`"b_focal"`.

```r
library(iopsychology)
library(readr)

data_path <- system.file("extdata", "SampleData.csv", package = "iopsychology")
sample_data <- read_csv(data_path, show_col_types = FALSE)

results <- PredBias(
  Combined = sample_data,
  x = "test",
  y = "criterion",
  output_file_name = "SampleOutput.csv",
  scatterplot_file_name = "SampleScatterPlot.png",
  predictor_histogram_file_name = "SamplePredictorHistogram.png",
  criterion_histogram_file_name = "SampleCriterionHistogram.png"
)
```

## PredBias

`PredBias()` runs predictive bias analyses including:

- Gulliksen-Wilks tests (slopes, intercepts, residual variance)
- Cleary regression tests (intercept, slope, and interaction)
- Separate-group regressions and Fisher's z correlation comparisons
- Residual metrics, Levene tests, and Breusch-Pagan homoscedasticity checks
- Group descriptives and effect sizes

Output is returned as a tibble and written to the CSV path you provide. Three
PNG plots (scatterplot and histograms) are also written to the paths you
specify.

## create_sim_dat

`create_sim_dat()` generates simulated reference/focal datasets for Monte Carlo
work. Fixed parameters must be single values; varying parameters can take
multiple values and produce separate CSV files for each combination.

```r
sim_results <- create_sim_dat(
  output_location = "C:/sim_output",
  reference_validity = c(0.1, 0.2),
  focal_validity = c(0.1, 0.2),
  focal_sd_y = seq(0.9, 1.1, by = 0.1),
  intercept_diffs = c(0.2, 0.3),
  D = 0.7,
  n = 5000,
  population_percentage = 0.5,
  iteration = 1,
  x_var_name = "test",
  y_var_name = "criterion"
)
```

Each run creates a dated folder such as `2026_08_20_1` containing
`simdat1.csv`, `simdat2.csv`, and a `parameter_key.txt` file mapping each
dataset to the parameters used.

## Dependencies

The following CRAN packages are required:

`broom`, `car`, `DescTools`, `dplyr`, `e1071`, `effsize`, `ggplot2`,
`lm.beta`, `lmtest`, `readr`, `tibble`, and `tidyr`.

## Authors

- Brian Costello
- Jeff Cucina (modifications for local empirical data files)
- Kimberly Wilson
- Phil Walmsley

## License

MIT
