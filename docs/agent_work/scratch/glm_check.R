## Standalone GLM verification script for external_evaluation.qmd
##
## Purpose: replicate analysis/external_evaluation.qmd's data loading and
## "Binomial Model" section EXACTLY (same code, adapted only to run outside
## Quarto), then print the values needed to resolve two manuscript
## placeholders and confirm two manuscript claims:
##   1. small_var_glm_lcb_pct / sv_glm_lcb_pct  -- minimum one-sided 95% LCB
##      across strata for the conservative (unsure = incorrect) primary GLM,
##      per bench type (small variant vs. structural variant).
##   2. GLM N + callset counts per model (manuscript claims: smvar N=350/5
##      callsets; stvar N=280/7 callsets).
##   3. Firth/brglm2 sensitivity failure count (manuscript claims: "5
##      failures across 350 observations" for the small-variant sensitivity
##      model).
##
## Run from repo root with renv active:
##   Rscript docs/agent_work/scratch/glm_check.R
##
## This script intentionally mirrors analysis/external_evaluation.qmd
## chunk-for-chunk (Data Loading, Strata Label Mapping, ci_data_prep,
## ci_models, ci_emmeans) so any mismatch reflects a real discrepancy in the
## qmd's own pipeline rather than a re-implementation difference.

set.seed(42) ## no randomization occurs in this pipeline; set defensively

source(here::here("analysis/_notebook_setup.R"))
analysis_setup(load_plot_themes = TRUE)

library(emmeans)
library(brglm2)

## ---------------------------------------------------------------------
## Data Loading (verbatim from analysis/external_evaluation.qmd)
## ---------------------------------------------------------------------

## Loading collaborator curations
eval_csvs <- list.files(
  path = here("data/external-evaluations"),
  pattern = "Miqa.csv",
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = TRUE
) %>%
  {
    set_names(., str_extract(., "(?<=evaluations/).*(?= On)"))
  }

## Subsetting small and structural variant tables based on number of rows
## - small variants expect 71, structural variants 41
## a few files have less as less than 5 variants per strata to subset
nlines <- eval_csvs %>% map_int(~ length(read_lines(.x)))

cols_to_keep <- c("callset", "correct.in.callset", "correct.in.benchmark", "chrom", "chromStart", "var_type", "strata")
smvar_evals <- eval_csvs[nlines > 50] %>%
  map_dfr(read_csv, name_repair = "universal", show_col_types = FALSE, .id = "callset") %>%
  select(all_of(cols_to_keep))

stvar_evals <- eval_csvs[nlines < 50] %>%
  map_dfr(read_csv,
    name_repair = "universal",
    show_col_types = FALSE,
    .id = "callset"
  )

## Load pre-combined evaluation data with JZ curations from local TSV files
smvar_evals_jz <- read_tsv(
  here("data/external-evaluations/Q100-ext-evals-2025-04-11-smvars.tsv"),
  show_col_types = FALSE,
  name_repair = "universal"
)

stvar_evals_jz <- read_tsv(
  here("data/external-evaluations/Q100-ext-evals-2025-04-11-stvars.tsv"),
  show_col_types = FALSE,
  name_repair = "universal"
)

## Manually Revising callset names for figures
stvar_callsets <- c(
  "Baylor-Hifi" = "sniffles2\nHiFi",
  "Baylor-Ont" = "sniffles2\nONT",
  "Comenius-Svdss2" = "svdss2",
  "Illumina-Dragen" = "DRAGEN\nIll",
  "Iter-Dragen" = "DRAGEN\nITER",
  "Pacbio-Sawfish" = "sawfish",
  "Ucla-Delly-Newvcf" = "Delly",
  "Ucla-Manta-Newvcf" = "Manta",
  "Ucla-Ensamble" = "Ensamble"
)

smvar_callsets <- c(
  "Illumina-Dragen" = "DRAGEN",
  "Jhu-Imputefirst" = "Impute",
  "Pacbio-Deepvariant" = "PB-DV",
  "Roche-Neusomatic" = "roche",
  "Ultima-Ug" = "ultima"
)

## Cleaned up tables with updated curations
## Note: TSV files already contain JZ curations, no join needed
smvar_curations <- full_join(smvar_evals, smvar_evals_jz) %>%
  rename(unique_to = label) %>%
  select(
    callset,
    GRCh38_chr, GRCh38_start, GRCh38_end,
    correct.in.benchmark, JZcorrect.in.benchmark,
    correct.in.callset, JZcorrect.in.callset,
    var_type, strata, unique_to
  ) %>%
  mutate(
    curation_benchmark = if_else(!is.na(JZcorrect.in.benchmark),
      JZcorrect.in.benchmark,
      correct.in.benchmark
    ),
    curation_benchmark = factor(curation_benchmark,
      levels = c("no", "unsure", "yes")
    ),
    curation_callset = if_else(!is.na(JZcorrect.in.callset),
      JZcorrect.in.callset,
      correct.in.callset
    ),
    curation_callset = factor(curation_callset,
      levels = c("no", "unsure", "yes")
    )
  ) %>%
  mutate(callset = smvar_callsets[callset])

stvar_curations <- full_join(stvar_evals, stvar_evals_jz) %>%
  rename(var_type = SVTYPE) %>%
  mutate(
    unique_to = if_else(unique_to == "benchmark", "bench", "query"),
    strata = str_glue("{TR_anno}\n{sv_cat}\n{unique_to}")
  ) %>%
  select(
    callset,
    GRCh38_chr, GRCh38_start, GRCh38_end,
    correct.in.benchmark, JZcorrect.in.benchmark,
    correct.in.callset, JZcorrect.in.callset,
    var_type, strata, unique_to, TR_anno, sv_cat, unique_to
  ) %>%
  mutate(
    curation_benchmark = if_else(!is.na(JZcorrect.in.benchmark),
      JZcorrect.in.benchmark,
      correct.in.benchmark
    ),
    curation_benchmark = factor(curation_benchmark,
      levels = c("no", "unsure", "yes")
    ),
    curation_callset = if_else(!is.na(JZcorrect.in.callset),
      JZcorrect.in.callset,
      correct.in.callset
    ),
    curation_callset = factor(curation_callset,
      levels = c("no", "unsure", "yes")
    )
  ) %>%
  mutate(callset = stvar_callsets[callset]) %>%
  filter(callset != "Manta", callset != "Delly")

## ---------------------------------------------------------------------
## Strata Label Mapping (verbatim)
## ---------------------------------------------------------------------

smvar_strata_labels <- c(
  "S01" = "SNV: TR",
  "S02" = "SNV: Low Map",
  "S03" = "SNV: XY nonPAR",
  "S04" = "SNV: Autosome",
  "S05" = "INDEL: HP >6bp",
  "S06" = "INDEL: TR",
  "S07" = "INDEL: Low Map",
  "S08" = "INDEL: XY nonPAR",
  "S09" = "INDEL: Autosome",
  "S10" = "SNV: TR",
  "S11" = "SNV: Low Map",
  "S12" = "SNV: XY nonPAR",
  "S13" = "SNV: Autosome",
  "S14" = "INDEL: HP >6bp",
  "S15" = "INDEL: TR",
  "S16" = "INDEL: Low Map",
  "S17" = "INDEL: XY nonPAR",
  "S18" = "INDEL: Autosome"
)

smvar_curations <- smvar_curations %>%
  mutate(
    unique_to = if_else(strata %in% paste0("S0", 1:9), "Query", "Benchmark"),
    strata_label = smvar_strata_labels[strata]
  )

## ---------------------------------------------------------------------
## Binomial Model: Probability Discrepancy is Error in Comparison Callset
## (verbatim ci_data_prep / ci_models / ci_emmeans chunks)
## ---------------------------------------------------------------------

prepare_smvar_features <- function(df) {
  df |>
    mutate(
      genomic_context = str_extract(strata_label, "(?<=: ).*"),
      genomic_context = factor(genomic_context,
        levels = c(
          "Autosome", "TR", "HP >6bp",
          "Low Map", "XY nonPAR"
        )
      ),
      genomic_context = droplevels(genomic_context),
      var_type = factor(var_type, levels = c("SNP", "INDEL")),
      unique_to = factor(unique_to, levels = c("Query", "Benchmark"))
    )
}

prepare_stvar_features <- function(df) {
  df |>
    mutate(
      TR_anno = factor(TR_anno, levels = c("nonTR", "TR")),
      sv_cat = factor(sv_cat, levels = c("DEL", "INS")),
      unique_to = factor(unique_to, levels = c("query", "bench"))
    )
}

## Primary analysis: unsure = incorrect (conservative)
smvar_model_df <- smvar_curations |>
  mutate(correct_binary = as.integer(curation_benchmark == "yes")) |>
  prepare_smvar_features()

stvar_model_df <- stvar_curations |>
  mutate(correct_binary = as.integer(curation_benchmark == "yes")) |>
  prepare_stvar_features()

## Sensitivity analysis: unsure = correct
smvar_sens_df <- smvar_curations |>
  mutate(correct_binary = as.integer(curation_benchmark %in% c("yes", "unsure"))) |>
  prepare_smvar_features()

stvar_sens_df <- stvar_curations |>
  mutate(correct_binary = as.integer(curation_benchmark %in% c("yes", "unsure"))) |>
  prepare_stvar_features()

## Primary models (standard GLM)
smvar_glm <- glm(correct_binary ~ var_type + genomic_context + unique_to,
  family = binomial, data = smvar_model_df
)
stvar_glm <- glm(correct_binary ~ TR_anno + sv_cat + unique_to,
  family = binomial, data = stvar_model_df
)

## Sensitivity models
## smvar sensitivity has quasi-complete separation (only 5 "no" curations
## with 0 in the Low Map context), so bias-reduced estimation is used.
smvar_sens_glm <- glm(correct_binary ~ var_type + genomic_context + unique_to,
  family = binomial, data = smvar_sens_df,
  method = "brglmFit"
)
stvar_sens_glm <- glm(correct_binary ~ TR_anno + sv_cat + unique_to,
  family = binomial, data = stvar_sens_df
)

## Helper to extract emmeans with one-sided 95% lower bound
extract_emmeans <- function(em, benchmarkset, feature_name) {
  as_tibble(summary(em, side = ">")) |>
    mutate(
      benchmarkset = benchmarkset,
      feature = feature_name,
      level = if ("1" %in% names(pick(everything()))) {
        "Overall"
      } else {
        as.character(pick(1)[[1]])
      }
    ) |>
    select(benchmarkset, feature, level, prob, SE, asymp.LCL)
}

compute_all_emmeans <- function(smvar_fit, stvar_fit) {
  bind_rows(
    extract_emmeans(
      emmeans(smvar_fit, ~1, type = "response", weights = "proportional"),
      "Small Variant", "Overall"
    ),
    extract_emmeans(
      emmeans(smvar_fit, ~var_type, type = "response", weights = "proportional"),
      "Small Variant", "Variant Type"
    ),
    extract_emmeans(
      emmeans(smvar_fit, ~genomic_context, type = "response", weights = "proportional"),
      "Small Variant", "Genomic Context"
    ),
    extract_emmeans(
      emmeans(smvar_fit, ~unique_to, type = "response", weights = "proportional"),
      "Small Variant", "Unique To"
    ),
    extract_emmeans(
      emmeans(stvar_fit, ~1, type = "response", weights = "proportional"),
      "Structural Variant", "Overall"
    ),
    extract_emmeans(
      emmeans(stvar_fit, ~TR_anno, type = "response", weights = "proportional"),
      "Structural Variant", "TR Annotation"
    ),
    extract_emmeans(
      emmeans(stvar_fit, ~sv_cat, type = "response", weights = "proportional"),
      "Structural Variant", "SV Category"
    ),
    extract_emmeans(
      emmeans(stvar_fit, ~unique_to, type = "response", weights = "proportional"),
      "Structural Variant", "Unique To"
    )
  )
}

primary_results <- compute_all_emmeans(smvar_glm, stvar_glm) |>
  mutate(analysis = "Conservative (unsure = incorrect)")
sens_results <- compute_all_emmeans(smvar_sens_glm, stvar_sens_glm) |>
  mutate(analysis = "Sensitive (unsure = correct)")

ci_results <- bind_rows(primary_results, sens_results)

## =======================================================================
## REPORTING (new code -- not in the qmd; prints the values requested by
## the manuscript-value-verification task)
## =======================================================================

cat("\n================ GLM Observation / Callset Counts ================\n")
cat(sprintf(
  "Small variant model (smvar_model_df): N = %d observations, %d callsets (%s)\n",
  nrow(smvar_model_df),
  length(unique(smvar_model_df$callset)),
  paste(sort(unique(smvar_model_df$callset)), collapse = ", ")
))
cat(sprintf(
  "Structural variant model (stvar_model_df): N = %d observations, %d callsets (%s)\n",
  nrow(stvar_model_df),
  length(unique(stvar_model_df$callset)),
  paste(sort(unique(stvar_model_df$callset)), collapse = ", ")
))
cat("Manuscript claims: smvar N=350 (5 callsets); stvar N=280 (7 callsets)\n")

cat("\n================ Per-Stratum One-Sided 95% LCBs (Primary/Conservative Model) ================\n")
primary_print <- primary_results |>
  mutate(across(c(prob, asymp.LCL), ~ round(.x * 100, 2))) |>
  select(benchmarkset, feature, level, prob_pct = prob, lcb_pct = asymp.LCL)
print(as.data.frame(primary_print), row.names = FALSE)

cat("\n================ Minimum LCB per Bench Type (resolves manuscript placeholders) ================\n")
min_lcb_all <- primary_results |>
  group_by(benchmarkset) |>
  summarise(min_lcb_pct = round(min(asymp.LCL) * 100, 2), .groups = "drop")
cat("-- Including 'Overall' row --\n")
print(as.data.frame(min_lcb_all), row.names = FALSE)

min_lcb_no_overall <- primary_results |>
  filter(feature != "Overall") |>
  group_by(benchmarkset) |>
  summarise(min_lcb_pct = round(min(asymp.LCL) * 100, 2), .groups = "drop")
cat("-- Excluding 'Overall' row (per-stratum only) --\n")
print(as.data.frame(min_lcb_no_overall), row.names = FALSE)

cat("\nsmall_var_glm_lcb_pct candidate(s) above (Small Variant row(s))\n")
cat("sv_glm_lcb_pct candidate(s) above (Structural Variant row(s))\n")

cat("\n================ Firth/brglm2 Sensitivity Failure Count ================\n")
smvar_sens_failures <- sum(smvar_sens_df$correct_binary == 0)
smvar_sens_n <- nrow(smvar_sens_df)
cat(sprintf(
  "Small variant sensitivity model (unsure = correct): %d failures across %d observations\n",
  smvar_sens_failures, smvar_sens_n
))
cat("Manuscript claims: \"5 failures across 350 observations\"\n")

stvar_sens_failures <- sum(stvar_sens_df$correct_binary == 0)
stvar_sens_n <- nrow(stvar_sens_df)
cat(sprintf(
  "(For reference) Structural variant sensitivity model (unsure = correct): %d failures across %d observations\n",
  stvar_sens_failures, stvar_sens_n
))

cat("\n================ Done ================\n")
