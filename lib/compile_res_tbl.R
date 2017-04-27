compile_res_tbl <- function(x, prefix = NULL, suffix = NULL) {
  res_tbl <- cbind(
    tibble(coef = rownames(summary(x)$coefficients)),
    summary(x)$conf.int[, c("exp(coef)", "lower .95", "upper .95"), drop = FALSE],
    summary(x)$coefficients[, "Pr(>|z|)", drop = FALSE]) %>%
    as_tibble() %>%
    (function(x) {
      names(x) <- c("coef", "HR", "CI_lo", "CI_up", "P")
      x
    }) %>%
    mutate(P = cut(
      P,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1),
      labels = c("***", "**", "*", ".", " "))
    ) %>%
    dplyr::select(coef, HR, P)

  if (!is.null(prefix)) {
    res_tbl <- rename(res_tbl, c(
      "coef" = "coef",
      "HR" = paste(prefix, "HR"),
      "P" = paste(prefix, "P")))
  }

  if (!is.null(suffix)) {
    res_tbl <- rename(res_tbl, c(
      "coef" = "coef",
      "HR" = paste("HR", suffix),
      "P" = paste("P", suffix)))
  }

  res_tbl
}
