FIG_DIR <- here::here("figures")

analysis_setup <- function(
  load_plot_themes = TRUE,
  load_sessioninfo = FALSE
) {
  packages <- c("tidyverse", "here", "patchwork")
  if (load_sessioninfo) {
    packages <- c(packages, "sessioninfo")
  }

  for (pkg in unique(packages)) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }

  source(here::here("R/data_loading.R"))
  if (load_plot_themes) {
    source(here::here("R/plot_themes.R"))
  }

  invisible(TRUE)
}
