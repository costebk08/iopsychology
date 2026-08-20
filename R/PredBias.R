#' Predictive Bias Analysis
#'
#' Runs a predictive bias analysis on a combined reference and focal group
#' dataset. Computes Gulliksen-Wilks tests, Cleary regression tests,
#' separate-group regressions, Fisher's z comparisons, residual metrics,
#' Levene tests, homoscedasticity diagnostics, and descriptive statistics.
#' Writes a CSV of results and optional PNG plots.
#'
#' @param Combined A data frame containing predictor, criterion, and group
#'   columns. The group column must be named `ref_foc` with values
#'   `"a_reference"` and `"b_focal"`.
#' @param x Character string naming the predictor column in `Combined`.
#' @param y Character string naming the criterion column in `Combined`.
#' @param output_file_name Path for the CSV file of analysis results.
#' @param scatterplot_file_name Path for the scatterplot PNG file.
#' @param predictor_histogram_file_name Path for the predictor histogram PNG.
#' @param criterion_histogram_file_name Path for the criterion histogram PNG.
#'
#' @return A tibble containing all computed metrics (also written to
#'   `output_file_name`).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data_path <- system.file("extdata", "SampleData.csv", package = "iopsychology")
#' sample_data <- readr::read_csv(data_path, show_col_types = FALSE)
#'
#' PredBias(
#'   Combined = sample_data,
#'   x = "test",
#'   y = "criterion",
#'   output_file_name = "SampleOutput.csv",
#'   scatterplot_file_name = "SampleScatterPlot.png",
#'   predictor_histogram_file_name = "SamplePredictorHistogram.png",
#'   criterion_histogram_file_name = "SampleCriterionHistogram.png"
#' )
#' }
PredBias <- function(Combined, x, y, output_file_name,
                     scatterplot_file_name, predictor_histogram_file_name,
                     criterion_histogram_file_name) {
  Combined$x <- unlist(Combined[[x]])
  Combined$y <- unlist(Combined[[y]])

  Combined <- dplyr::mutate(Combined, case = 1)

  reference_dat <- Combined[Combined$ref_foc == "a_reference", ]
  focal_dat <- Combined[Combined$ref_foc == "b_focal", ]

  GMX <- mean(Combined$x)
  GMY <- mean(Combined$y)
  n_total <- nrow(Combined)

  gw_test <-
    Combined %>%
    dplyr::group_by(ref_foc) %>%
    dplyr::summarise(
      mean_x = mean(x),
      mean_y = mean(y),
      n = sum(case),
      ex = sum((x - mean_x)^2),
      ey = sum((y - mean_y)^2),
      exy = sum((x - mean_x) * (y - mean_y)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      C7 = ey - ((exy^2) / ex)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      C8 = sum(C7) / sum(n),
      C9 = sum(n) * log(C8),
      C10 = C7 / n,
      C11 = n * log(C10),
      GWA = C9 - sum(C11),
      E = sum(exy),
      AB = mean_y - GMY,
      CD = mean_x - GMX,
      SDOT = sum(ey) - ((E^2) / sum(ex)),
      GWB = sum(n) * (log(SDOT / sum(n)) - log(sum(C7) / sum(n))),
      SA = sum(ey) + sum(n * (AB^2)),
      SB = sum(exy) + sum(n * AB * CD),
      SC = sum(ex) + sum(n * (CD^2)),
      SDASH = SA - ((SB^2) / SC),
      GWC = n_total * (log(SDASH / n_total) - log(SDOT / n_total))
    ) %>%
    dplyr::summarise(
      GWA = mean(GWA),
      GWA_p = stats::pchisq(GWA, 1, lower.tail = FALSE),
      GWB = mean(GWB),
      GWB_p = stats::pchisq(GWB, 1, lower.tail = FALSE),
      GWC = mean(GWC),
      GWC_p = stats::pchisq(GWC, 1, lower.tail = FALSE),
      .groups = "drop"
    )

  reg1 <- stats::lm(y ~ x, data = Combined)
  reg1_summary <- summary(reg1)
  step1_r2 <- reg1_summary$r.squared
  step1_beta <- lm.beta::lm.beta(reg1)
  step1_r <- sqrt(step1_r2)

  reg2 <- stats::lm(y ~ x + ref_foc, data = Combined)
  reg2_summary <- summary(reg2)
  step2_r2 <- reg2_summary$r.squared
  step2_r <- sqrt(step2_r2)
  step2_beta <- lm.beta::lm.beta(reg2)
  step2_step1_comp <- stats::anova(reg1, reg2)["2", ]
  step2_delta_r <- step2_r - step1_r
  step2_delta_r2 <- step2_r2 - step1_r2
  step2_b_intercept <- step2_beta$coefficients["(Intercept)"]
  step2_b_main <- step2_beta$coefficients["ref_focb_focal"]
  step2_beta_main <- step2_beta$standardized.coefficients["ref_focb_focal"]
  step2_delta_r_f <- step2_step1_comp$`F`
  step2_delta_r_p <- step2_step1_comp$`Pr(>F)`

  reg3 <- stats::lm(y ~ x * ref_foc, data = Combined)
  reg3_summary <- summary(reg3)
  step3_r2 <- reg3_summary$r.squared
  step3_r <- sqrt(step3_r2)
  step3_beta <- lm.beta::lm.beta(reg3)
  step3_step2_comp <- stats::anova(reg2, reg3)["2", ]
  step3_delta_r <- step3_r - step2_r
  step3_delta_r2 <- step3_r2 - step2_r2
  step3_b_interaction <- step3_beta$coefficients["x:ref_focb_focal"]
  step3_beta_interaction <- step3_beta$standardized.coefficients["x:ref_focb_focal"]
  step3_delta_r_f <- step3_step2_comp$`F`
  step3_delta_r_p <- step3_step2_comp$`Pr(>F)`

  cohen_f2 <- step3_delta_r2 / (1 - (step1_r2 + step2_delta_r2 + step3_delta_r2))
  liu_yuan <- step3_delta_r2 / (step1_r2 + step3_delta_r2)

  cleary_test <-
    tibble::tibble(
      step1_r = step1_r,
      step2_delta_r = step2_delta_r,
      step2_delta_r2 = step2_delta_r2,
      step2_b_intercept = step2_b_intercept,
      step2_b_main = step2_b_main,
      step2_beta_main = step2_beta_main,
      step2_delta_r_f = step2_delta_r_f,
      step2_delta_r_p = step2_delta_r_p,
      step3_delta_r = step3_delta_r,
      step3_delta_r2 = step3_delta_r2,
      step3_b_interaction = step3_b_interaction,
      step3_beta_interaction = step3_beta_interaction,
      step3_delta_r_f = step3_delta_r_f,
      step3_delta_r_p = step3_delta_r_p,
      cohen_f2 = cohen_f2,
      liu_yuan = liu_yuan
    )

  ref_reg <- stats::lm(y ~ x, data = reference_dat)
  ref_reg_beta <- lm.beta::lm.beta(ref_reg)
  ref_constant <- ref_reg_beta$coefficients["(Intercept)"]
  ref_b <- ref_reg_beta$coefficients["x"]
  ref_beta <- ref_reg_beta$standardized.coefficients["x"]
  ref_p <- broom::tidy(ref_reg) %>%
    dplyr::filter(term == "x") %>%
    dplyr::pull(p.value)

  focal_reg <- stats::lm(y ~ x, data = focal_dat)
  focal_reg_beta <- lm.beta::lm.beta(focal_reg)
  focal_constant <- focal_reg_beta$coefficients["(Intercept)"]
  focal_b <- focal_reg_beta$coefficients["x"]
  focal_beta <- focal_reg_beta$standardized.coefficients["x"]
  focal_p <- broom::tidy(focal_reg) %>%
    dplyr::filter(term == "x") %>%
    dplyr::pull(p.value)

  sep_reg_tests <-
    tibble::tibble(
      ref_constant = ref_constant,
      ref_b = ref_b,
      ref_beta = ref_beta,
      ref_p = ref_p,
      focal_constant = focal_constant,
      focal_b = focal_b,
      focal_beta = focal_beta,
      focal_p = focal_p
    )

  levene_test_mean <-
    broom::tidy(
      car::leveneTest(y ~ as.factor(ref_foc), Combined, center = mean)
    ) %>%
    dplyr::filter(!is.na(statistic)) %>%
    dplyr::transmute(
      levene_mean_f = statistic,
      levene_mean_p = p.value
    )

  levene_test_median <-
    broom::tidy(
      car::leveneTest(y ~ as.factor(ref_foc), Combined)
    ) %>%
    dplyr::filter(!is.na(statistic)) %>%
    dplyr::transmute(
      levene_median_f = statistic,
      levene_median_p = p.value
    )

  n <- sum(Combined$case)

  residual_metrics <-
    tibble::tibble(
      ref_foc = Combined$ref_foc,
      residuals = step1_beta$residuals,
      residuals_z = scale(step1_beta$residuals)[1:n]
    ) %>%
    dplyr::mutate(
      ref_foc = ifelse(ref_foc == "a_reference", "reference", "focal")
    ) %>%
    dplyr::group_by(ref_foc) %>%
    dplyr::summarise(
      residuals_m = mean(residuals),
      residuals_sd = stats::sd(residuals),
      residuals_z_m = mean(residuals_z),
      residuals_z_sd = stats::sd(residuals_z),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = -ref_foc,
      names_to = "key",
      values_to = "value"
    ) %>%
    dplyr::transmute(
      metric = paste(ref_foc, key, sep = "_"),
      value = value
    ) %>%
    tidyr::pivot_wider(
      names_from = metric,
      values_from = value
    )

  n_reference <- nrow(reference_dat)
  n_focal <- nrow(focal_dat)
  se_z_corr <- sqrt((1 / (n_focal - 3)) + (1 / (n_reference - 3)))
  se_z_corr_results <- tibble::tibble(se_z_corr = se_z_corr)

  fisher_z_reference <- DescTools::FisherZ(ref_beta)
  fisher_z_focal <- DescTools::FisherZ(focal_beta)
  fisher_z_corr_compare_z <- (fisher_z_reference - fisher_z_focal) / se_z_corr
  fisher_z_corr_compare_p <- 2 * (1 - stats::pnorm(abs(fisher_z_corr_compare_z)))
  fisher_z_corr_compare_q <- fisher_z_reference - fisher_z_focal

  fisher_z_test <-
    tibble::tibble(
      fisher_z_reference = fisher_z_reference,
      fisher_z_focal = fisher_z_focal,
      fisher_z_corr_compare_z = fisher_z_corr_compare_z,
      fisher_z_corr_compare_p = fisher_z_corr_compare_p,
      fisher_z_corr_compare_q = fisher_z_corr_compare_q
    )

  MeanXreference <- mean(reference_dat$x)
  MeanYreference <- mean(reference_dat$y)
  MeanXfocal <- mean(focal_dat$x)
  MeanYfocal <- mean(focal_dat$y)

  MinXreference <- min(reference_dat$x)
  MinYreference <- min(reference_dat$y)
  MinXfocal <- min(focal_dat$x)
  MinYfocal <- min(focal_dat$y)

  MaxXreference <- max(reference_dat$x)
  MaxYreference <- max(reference_dat$y)
  MaxXfocal <- max(focal_dat$x)
  MaxYfocal <- max(focal_dat$y)

  SDXreference <- stats::sd(reference_dat$x)
  SDYreference <- stats::sd(reference_dat$y)
  SDXfocal <- stats::sd(focal_dat$x)
  SDYfocal <- stats::sd(focal_dat$y)

  skewnessXreference <- e1071::skewness(reference_dat$x)
  skewnessYreference <- e1071::skewness(reference_dat$y)
  skewnessXfocal <- e1071::skewness(focal_dat$x)
  skewnessYfocal <- e1071::skewness(focal_dat$y)

  kurtosisXreference <- e1071::kurtosis(reference_dat$x)
  kurtosisYreference <- e1071::kurtosis(reference_dat$y)
  kurtosisXfocal <- e1071::kurtosis(focal_dat$x)
  kurtosisYfocal <- e1071::kurtosis(focal_dat$y)

  t_test_results <-
    broom::tidy(stats::t.test(x ~ ref_foc, data = Combined)) %>%
    dplyr::transmute(
      t_test_x = statistic,
      t_test_x_p = p.value
    ) %>%
    dplyr::bind_cols(
      broom::tidy(stats::t.test(y ~ ref_foc, data = Combined)) %>%
        dplyr::transmute(
          t_test_y = statistic,
          t_test_y_p = p.value
        )
    )

  cohen_d_x <- effsize::cohen.d(reference_dat$x, focal_dat$x)
  hedges_g_x <- effsize::cohen.d(reference_dat$x, focal_dat$x, hedges.correction = TRUE)
  cohen_d_y <- effsize::cohen.d(reference_dat$y, focal_dat$y)
  hedges_g_y <- effsize::cohen.d(reference_dat$y, focal_dat$y, hedges.correction = TRUE)

  cohen_d_x_results <-
    tibble::tibble(
      cohen_d_x_estimate = cohen_d_x$estimate,
      cohen_d_x_within_group_sd = cohen_d_x$sd
    )

  cohen_d_y_results <-
    tibble::tibble(
      cohen_d_y_estimate = cohen_d_y$estimate,
      cohen_d_y_within_group_sd = cohen_d_y$sd
    )

  hedges_g_x_results <-
    tibble::tibble(
      hedges_g_x_estimate = hedges_g_x$estimate
    )

  hedges_g_y_results <-
    tibble::tibble(
      hedges_g_y_estimate = hedges_g_y$estimate
    )

  homoscedasticity_tests <-
    broom::tidy(lmtest::bptest(reg1)) %>%
    dplyr::transmute(
      breusch_pagan_lm1 = statistic,
      breusch_pagan_p_lm1 = p.value
    ) %>%
    dplyr::bind_cols(
      broom::tidy(lmtest::bptest(reg2)) %>%
        dplyr::transmute(
          breusch_pagan_lm2 = statistic,
          breusch_pagan_p_lm2 = p.value
        )
    ) %>%
    dplyr::bind_cols(
      broom::tidy(lmtest::bptest(reg3)) %>%
        dplyr::transmute(
          breusch_pagan_lm3 = statistic,
          breusch_pagan_p_lm3 = p.value
        )
    )

  referenceValidity <- ref_beta
  focalValidity <- focal_beta

  sep_descriptives <-
    tibble::tibble(
      n_reference = n_reference,
      MeanXreference = MeanXreference,
      SDXreference = SDXreference,
      MinXreference = MinXreference,
      MaxXreference = MaxXreference,
      skewnessXreference = skewnessXreference,
      kurtosisXreference = kurtosisXreference,
      MeanYreference = MeanYreference,
      SDYreference = SDYreference,
      MinYreference = MinYreference,
      MaxYreference = MaxYreference,
      skewnessYreference = skewnessYreference,
      kurtosisYreference = kurtosisYreference,
      referenceValidity = referenceValidity,
      n_focal = n_focal,
      MeanXfocal = MeanXfocal,
      SDXfocal = SDXfocal,
      MinXfocal = MinXfocal,
      MaxXfocal = MaxXfocal,
      skewnessXfocal = skewnessXfocal,
      kurtosisXfocal = kurtosisXfocal,
      MeanYfocal = MeanYfocal,
      SDYfocal = SDYfocal,
      MinYfocal = MinYfocal,
      MaxYfocal = MaxYfocal,
      skewnessYfocal = skewnessYfocal,
      kurtosisYfocal = kurtosisYfocal,
      focalValidity = focalValidity
    )

  final_dat <-
    dplyr::bind_cols(
      x = x,
      y = y,
      sep_descriptives,
      t_test_results,
      cohen_d_x_results,
      hedges_g_x_results,
      cohen_d_y_results,
      hedges_g_y_results,
      gw_test,
      cleary_test,
      sep_reg_tests,
      fisher_z_test,
      se_z_corr_results,
      residual_metrics,
      levene_test_mean,
      levene_test_median,
      homoscedasticity_tests
    )

  readr::write_csv(final_dat, output_file_name)

  ref_r <- round(ref_reg_beta$standardized.coefficients["x"], digits = 2)
  ref_x <- round(ref_reg$coefficients["x"], digits = 2)
  ref_intercept <- round(ref_reg$coefficients["(Intercept)"], digits = 2)

  focal_r <- round(focal_reg_beta$standardized.coefficients["x"], digits = 2)
  focal_x <- round(focal_reg$coefficients["x"], digits = 2)
  focal_intercept <- round(focal_reg$coefficients["(Intercept)"], digits = 2)

  grDevices::png(
    filename = scatterplot_file_name,
    width = 480,
    height = 480,
    pointsize = 12,
    bg = "white",
    res = NA
  )

  scatterplot <-
    Combined %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, color = ref_foc)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
    ggplot2::labs(
      caption = paste0(
        "Ref: R = ", ref_r, "; y = ", ref_x, "*x + ", ref_intercept,
        "; ", "Focal: R = ", focal_r, "; y = ", focal_x, "*x + ", focal_intercept
      )
    )

  print(scatterplot + ggplot2::theme_bw())
  grDevices::dev.off()

  combined_dat1 <- Combined %>%
    dplyr::mutate(
      ref_foc = dplyr::case_when(
        ref_foc == "a_reference" ~ "Reference",
        ref_foc == "b_focal" ~ "Focal"
      )
    ) %>%
    dplyr::rename(`Focal/Reference` = ref_foc)

  grDevices::png(
    filename = predictor_histogram_file_name,
    width = 480,
    height = 480,
    pointsize = 12,
    bg = "white",
    res = NA
  )

  hist_x <- combined_dat1 %>%
    ggplot2::ggplot(ggplot2::aes(x = x, fill = `Focal/Reference`, color = `Focal/Reference`)) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.5) +
    ggplot2::labs(
      title = "Predictor Histogram",
      x = "Predictor Values",
      y = "Count"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(hist_x)
  grDevices::dev.off()

  grDevices::png(
    filename = criterion_histogram_file_name,
    width = 480,
    height = 480,
    pointsize = 12,
    bg = "white",
    res = NA
  )

  hist_y <- combined_dat1 %>%
    ggplot2::ggplot(ggplot2::aes(x = y, fill = `Focal/Reference`, color = `Focal/Reference`)) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.5) +
    ggplot2::labs(
      title = "Criterion Histogram",
      x = "Criterion Values",
      y = "Count"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(hist_y)
  grDevices::dev.off()

  final_dat
}
