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
