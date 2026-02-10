# Source schema and cache infrastructure
source(here::here("R/schemas.R"))
source(here::here("R/cache.R"))

# Ensure magrittr pipe works when sourcing this file in a clean R session
if (!exists("%>%", mode = "function")) {
  `%>%` <- magrittr::`%>%`
}

#' Parse Benchmark Identifier from File Path
#'
#' Extracts benchmark metadata (version, reference, variant type) from a file path
#' or benchmark ID string using the standard naming pattern:
#' `{benchmark_version}_{ref}_{bench_type}`.
#'
#' @param file_path Character string containing file path or benchmark ID
#'
#' @return Named list with elements:
#'   - bench_version: Benchmark version (e.g., "v5.0q")
#'   - ref: Reference genome (e.g., "GRCh38")
#'   - bench_type: benchmark set type small or structural variants (factored: smvar, stvar)
#'
#' @examples
#' \dontrun{
#' parse_benchmark_id("results/var_counts/v5.0q_GRCh38_smvar/genomic_context_combined_metrics.csv")
#' # Returns list(bench_version = "v5.0q", ref = "GRCh38", bench_type = "smvar")
#' }
#'
#' @export
parse_benchmark_id <- function(file_path) {
  if (!is.character(file_path) || length(file_path) != 1L || is.na(file_path)) {
    stop(
      "file_path must be a single non-missing character string.",
      call. = FALSE
    )
  }

  # Extract all benchmark-like tokens and use the last one (works for full paths)
  benchmark_tokens <- stringr::str_extract_all(
    file_path,
    "[A-Za-z0-9.-]+_[A-Za-z0-9.-]+_(smvar|stvar)"
  )[[1]]

  if (length(benchmark_tokens) == 0) {
    stop(
      glue::glue(
        "Could not parse benchmark ID from: {file_path}\n",
        "Expected pattern: [benchmark_version]_[reference]_[smvar|stvar]"
      ),
      call. = FALSE
    )
  }

  benchmark_id <- benchmark_tokens[length(benchmark_tokens)]
  components <- stringr::str_match(
    benchmark_id,
    "^([A-Za-z0-9.-]+)_([A-Za-z0-9.-]+)_(smvar|stvar)$"
  )[1, ]

  list(
    bench_version = components[2L],
    ref = components[3L],
    bench_type = components[4L]
  )
}


#' Parse Snakemake Configuration for Benchmarks
#'
#' Reads the Snakemake pipeline config.yaml and extracts benchmark metadata.
#' Used for dynamic validation without hardcoding expected values.
#'
#' @param config_path Path to config.yaml file. Default: `here::here("config/config.yaml")`
#'
#' @return Named list with:
#'   - benchmarks: Vector of benchmark IDs
#'   - num_benchmarks: Number of benchmarks
#'   - references: Vector of unique reference genomes
#'   - bench_type: benchmark set type small or structural variants (factored: smvar, stvar)
#'
#' @examples
#' \dontrun{
#' config <- parse_pipeline_config()
#' config$num_benchmarks # Get number of benchmarks
#' config$benchmarks # Get all benchmark IDs
#' }
#'
#' @export
parse_pipeline_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    config_path <- here::here("config/config.yaml")
  }

  if (!fs::file_exists(config_path)) {
    stop(glue::glue("Config file not found: {config_path}"), call. = FALSE)
  }

  # Parse YAML
  config <- yaml::read_yaml(config_path)

  # Extract benchmarks section
  if (is.null(config$benchmarksets)) {
    stop("No 'benchmarksets' section found in config.yaml", call. = FALSE)
  }

  benchmarks <- names(config$benchmarksets)

  # Parse benchmark metadata
  benchmark_meta <- benchmarks %>%
    purrr::map_df(
      ~ {
        meta <- parse_benchmark_id(.x)
        tibble::tibble(
          benchmark_id = .x,
          bench_version = meta$bench_version,
          ref = meta$ref,
          bench_type = meta$bench_type
        )
      }
    )

  list(
    benchmarks = benchmarks,
    num_benchmarks = length(benchmarks),
    benchmark_meta = benchmark_meta,
    references = unique(benchmark_meta$ref),
    bench_types = unique(benchmark_meta$bench_type),
    num_contexts = 6 # HP, MAP, SD, SD10kb, TR, TR10kb
  )
}

## Helper functions to standardize variables and factor levels
refs <- c("GRCh37", "GRCh38", "CHM13v2.0")
bench_version <- c("v0.6", "v4.2.1", "v5.0q")
bench_type <- c("smvar", "stvar")
chromosomes <- stringr::str_c("chr", c(1:22, "X", "Y"))
context_names <- c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
var_types <- c("SNV", "INDEL", "DEL", "INS")

std_ <- function(x, levels) {
  factor(x, levels = levels, labels = levels, ordered = TRUE)
}
std_references <- function(x) std_(x, levels = refs)
std_bench_versions <- function(x) std_(x, levels = bench_version)
std_bench_types <- function(x) std_(x, levels = bench_type)
std_context_name <- function(x) std_(x, levels = context_names)
std_chrom <- function(x, chromx = chromosomes) {
  dplyr::if_else(
    stringr::str_detect(x, "^chr"),
    x,
    stringr::str_c("chr", x)
  ) %>%
    std_(levels = chromx)
}
std_var_type <- function(x) std_(x, levels = var_types)

.write_cache_safely <- function(
  data,
  dataset_name,
  source_files,
  cache_params = list()
) {
  tryCatch(
    write_cache(data, dataset_name, source_files, cache_params),
    error = function(e) {
      warning(
        glue::glue("Failed to write cache for {dataset_name}: {e$message}"),
        call. = FALSE
      )
    }
  )
}

.verify_region_dataframe <- function(
  regions_df,
  bench_version_levels,
  ref_levels,
  bench_type_levels,
  source_label = "region BED files"
) {
  required_cols <- c(
    "bench_version",
    "ref",
    "bench_type",
    "chrom",
    "start",
    "end",
    "interval_size"
  )

  missing_cols <- setdiff(required_cols, names(regions_df))
  if (length(missing_cols) > 0) {
    stop(
      glue::glue(
        "Missing required columns while loading {source_label}: ",
        "{paste(missing_cols, collapse = ', ')}"
      ),
      call. = FALSE
    )
  }

  if (nrow(regions_df) == 0) {
    stop(glue::glue("No regions were loaded from {source_label}."), call. = FALSE)
  }

  if (
    anyNA(regions_df$bench_version) ||
      anyNA(regions_df$ref) ||
      anyNA(regions_df$bench_type)
  ) {
    stop(
      glue::glue(
        "Metadata parsing failed for at least one entry in {source_label}. ",
        "Expected bench_version in [{paste(bench_version_levels, collapse = ', ')}], ",
        "ref in [{paste(ref_levels, collapse = ', ')}], and ",
        "bench_type in [{paste(bench_type_levels, collapse = ', ')}]."
      ),
      call. = FALSE
    )
  }

  if (anyNA(regions_df$chrom) || any(!stringr::str_detect(regions_df$chrom, "^chr"))) {
    stop(
      glue::glue(
        "Chromosome normalization failed for {source_label}. ",
        "All chromosome values must start with 'chr'."
      ),
      call. = FALSE
    )
  }

  if (anyNA(regions_df$start) || any(regions_df$start < 0)) {
    stop(
      glue::glue(
        "Invalid start coordinates detected in {source_label}. ",
        "Start positions must be non-negative."
      ),
      call. = FALSE
    )
  }

  if (anyNA(regions_df$end) || any(regions_df$end <= regions_df$start)) {
    stop(
      glue::glue(
        "Invalid end coordinates detected in {source_label}. ",
        "End must be greater than start for all intervals."
      ),
      call. = FALSE
    )
  }

  if (anyNA(regions_df$interval_size) || any(regions_df$interval_size <= 0)) {
    stop(
      glue::glue(
        "Invalid interval sizes detected in {source_label}. ",
        "All interval sizes must be > 0."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.read_region_bed_files <- function(
  region_files,
  bench_version_levels,
  ref_levels,
  bench_type_levels = c("smvar", "stvar")
) {
  if (length(region_files) == 0) {
    stop("No region files provided.", call. = FALSE)
  }

  if (is.null(names(region_files)) || any(names(region_files) == "")) {
    stop("Region files must be a named character vector.", call. = FALSE)
  }

  regions_df <- region_files %>%
    purrr::map_dfr(
      ~ vroom::vroom(
        .x,
        col_names = c("chrom", "start", "end"),
        col_types = "cii",
        show_col_types = FALSE
      ),
      .id = "benchmark_id"
    ) %>%
    dplyr::mutate(
      bench_meta = purrr::map(benchmark_id, parse_benchmark_id),
      bench_version = purrr::map_chr(bench_meta, "bench_version"),
      ref = purrr::map_chr(bench_meta, "ref"),
      bench_type = purrr::map_chr(bench_meta, "bench_type")
    ) %>%
    dplyr::mutate(
      chrom = dplyr::if_else(
        stringr::str_detect(chrom, "^chr"),
        chrom,
        stringr::str_c("chr", chrom)
      ),
      interval_size = end - start,
      bench_version = factor(bench_version, levels = bench_version_levels),
      ref = factor(ref, levels = ref_levels),
      bench_type = factor(bench_type, levels = bench_type_levels)
    ) %>%
    dplyr::select(
      bench_version,
      ref,
      bench_type,
      chrom,
      start,
      end,
      interval_size
      )

  .verify_region_dataframe(
    regions_df = regions_df,
    bench_version_levels = bench_version_levels,
    ref_levels = ref_levels,
    bench_type_levels = bench_type_levels
  )

  regions_df
}

#' Load Genomic Context Metrics
#'
#' Loads primary analysis data files containing per-genomic-context metrics and variant counts.
#' These are the smallest, fastest-loading files and should be used for most analyses.
#'
#' @param results_dir Path to results directory. Default: `here::here("results")`
#' @param benchmark_filter Optional character vector of benchmark IDs to filter results
#'
#' @return Tibble with columns:
#'   - bench_version: Benchmark version (factored: v0.6, v4.2.1, v5.0q)
#'   - ref: Reference genome (factored: GRCh37, GRCh38, CHM13v2.0)
#'   - bench_type: benchmark set type small or structural variants (factored: smvar, stvar)
#'   - context_name: Genomic context region name (factored: HP, MAP, SD, SD10kb, TR, TR10kb)
#'   - context_bp: Total bases in genomic context
#'   - intersect_bp: Bases overlapping with benchmark
#'   - pct_of_context: Percentage of genomic context overlapping benchmark
#'   - pct_of_bench: Percentage of benchmark overlapping genomic context
#'   - total_variants: Total variant count
#'   - snv_count: Single Nucleotide Variant count (SNV)
#'   - indel_count: INDEL count
#'   - del_count: Deletion count (structural variants only)
#'   - ins_count: Insertion count (structural variants only)
#'   - complex_count: Complex variant count
#'   - other_count: Other variant count
#'   - variant_density_per_mb: Variants per megabase
#'
#' @examples
#' \dontrun{
#' # Load all metrics
#' metrics <- load_genomic_context_metrics()
#'
#' # Load specific benchmarks
#' metrics <- load_genomic_context_metrics(
#'   benchmark_filter = c("v5.0q_GRCh38_smvar", "v5.0q_GRCh37_stvar")
#' )
#' }
#'
#' @export
load_genomic_context_metrics <- function(
  results_dir = NULL,
  benchmark_filter = NULL
) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Find all genomic context metrics files
  metrics_files <- fs::dir_ls(
    results_dir,
    recurse = TRUE,
    glob = "**/genomic_context_combined_metrics.csv"
  )

  if (length(metrics_files) == 0) {
    stop(
      glue::glue(
        "No genomic_context_combined_metrics.csv files found in {results_dir}"
      ),
      call. = FALSE
    )
  }

  # Read all files and add benchmark metadata
  metrics_df <- metrics_files %>%
    purrr::map_dfr(function(file) {
      # Get benchmark ID from parent directory
      benchmark_dir <- fs::path_dir(file)
      benchmark_id <- fs::path_file(benchmark_dir)

      # Parse benchmark metadata
      meta <- parse_benchmark_id(benchmark_id)

      # Read file with explicit column types
      file %>%
        vroom::vroom(
          col_types = vroom::cols(
            context_name = "c",
            context_bp = "d",
            intersect_bp = "d",
            pct_of_context = "d",
            pct_of_bench = "d",
            total_variants = "i",
            snp_count = "i",
            indel_count = "i",
            del_count = "i",
            ins_count = "i",
            complex_count = "i",
            other_count = "i",
            variant_density_per_mb = "d",
            .default = "c"
          ),
          show_col_types = FALSE
        ) %>%
        tibble::add_column(
          bench_version = meta$bench_version,
          ref = meta$ref,
          bench_type = meta$bench_type,
          .before = 1
        )
    })

  # Apply benchmark filter if provided
  if (!is.null(benchmark_filter)) {
    metrics_df <- metrics_df %>%
      dplyr::filter(
        paste(bench_version, ref, bench_type, sep = "_") %in% benchmark_filter
      )

    if (nrow(metrics_df) == 0) {
      warning(
        glue::glue(
          "No rows remaining after applying benchmark filter. ",
          "Check that benchmark IDs match available data."
        )
      )
    }
  }

  # Rename snp_count to snv_count (SNV = Single Nucleotide Variant)
  # and factor categorical variables with standard levels
  metrics_df <- metrics_df %>%
    dplyr::rename(snv_count = snp_count) %>%
    dplyr::mutate(
      bench_version = std_bench_versions(bench_version),
      ref = std_references(ref),
      bench_type = std_bench_types(bench_type),
      context_name = std_context_name(context_name)
    )
  return(metrics_df)
}


#' Load Exclusion Metrics
#'
#' Loads exclusion overlap tables (v5.0q benchmarks only).
#' For other benchmark versions, this function will warn and return empty tibble.
#'
#' @param results_dir Path to results directory. Default: `here::here("results")`
#'
#' @return Tibble with columns:
#'   - bench_version: Benchmark version
#'   - ref: Reference genome
#'   - bench_type: benchmark set type small or structural variants (factored: smvar, stvar)
#'   - exclusions: Exclusion region name
#'   - exclusion_bp: Total bases in exclusion
#'   - intersect_bp: Bases overlapping with benchmark
#'   - pct_of_exclusion: Percentage of exclusion overlapping benchmark
#'   - pct_of_dip: Percentage of diploid genome overlapping benchmark
#'
#' @examples
#' \dontrun{
#' exclusions <- load_exclusion_metrics()
#' }
#'
#' @export
load_exclusion_metrics <- function(results_dir = NULL) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Find all exclusion metrics files
  exclusion_files <- fs::dir_ls(
    results_dir,
    recurse = TRUE,
    glob = "**/exclusions_intersection_table.csv",
    fail = FALSE
  )

  if (length(exclusion_files) == 0) {
    warning(
      "No exclusion files found. These are only available for v5.0q benchmarks.",
      call. = FALSE
    )
    return(tibble::tibble())
  }

  # Read all files and add benchmark metadata
  exclusion_df <- exclusion_files %>%
    purrr::map_dfr(function(file) {
      # Get benchmark ID from parent directory
      benchmark_dir <- fs::path_dir(file)
      benchmark_id <- fs::path_file(benchmark_dir)

      # Parse benchmark metadata
      meta <- parse_benchmark_id(benchmark_id)

      # Read file
      file %>%
        vroom::vroom(
          col_types = vroom::cols(
            exclusions = "c",
            exclusion_bp = "d",
            intersect_bp = "d",
            pct_of_exclusion = "d",
            pct_of_dip = "d",
            .default = "c"
          ),
          show_col_types = FALSE
        ) %>%
        tibble::add_column(
          bench_version = meta$bench_version,
          ref = meta$ref,
          bench_type = meta$bench_type,
          .before = 1
        )
    })

  return(exclusion_df)
}


#' Load Reference Genome Sizes
#'
#' Loads reference genome size files and calculates assembled bases
#' (excluding Ns).
#'
#' @param results_dir Path to results directory. Default: `here::here("results")`
#'
#' @return Tibble with columns:
#'   - ref: Reference genome name
#'   - chrom: Chromosome name (with "chr" prefix)
#'   - length: Total chromosome length
#'   - ns: Number of N bases
#'   - asm_bp: Assembled bases (length - ns)
#'
#' @examples
#' \dontrun{
#' ref_sizes <- load_reference_sizes()
#' }
#'
#' @export
load_reference_sizes <- function(results_dir = NULL) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Find all size files
  size_files <- fs::dir_ls(
    fs::path(results_dir, "ref_genome_sizes"),
    glob = "*_size.tsv",
    fail = FALSE
  )

  if (length(size_files) == 0) {
    stop(
      glue::glue(
        "No reference size files found in {results_dir}/ref_genome_sizes"
      ),
      call. = FALSE
    )
  }

  # Read all files
  size_df <- size_files %>%
    purrr::map_dfr(function(file) {
      # Extract reference name from filename
      ref_name <- fs::path_file(file) %>%
        stringr::str_remove("_size.tsv")

      file %>%
        vroom::vroom(
          col_names = c("chrom", "length", "atgcs", "ns"),
          col_types = "ciii",
          show_col_types = FALSE,
          comment = "#"
        ) %>%
        tibble::add_column(ref = ref_name, .before = 1) %>%
        dplyr::mutate(
          # Standardize chromosome naming
          chrom = dplyr::if_else(
            stringr::str_detect(chrom, "^chr"),
            chrom,
            paste0("chr", chrom)
          ),
          # Calculate assembled bases
          asm_bp = length - ns
        )
    })

  return(size_df)
}

## Reference Genome Sizes
## Get hg002q100v1.1 size
#' Get HG002 Q100 assembly size
#'
#' Retrieve the genome size associated with the HG002 Q100 maternal assembly used by the q100 variant benchmark.
#'
#' @details
#' Returns a single numeric value representing the total bp for the mat assembly.
#'
#' @return
#' genome in in base pairs (numeric)
#'
#' @examples
#' \dontrun{
#' size <- load_hg002q100_size()
#' print(size)
#' }
#'
load_hg002q100_size <- function(
  asm_version = "v1.1",
  asm_verison = asm_version,
  fai_url = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.mat_Y_EBV_MT.fasta.gz.fai",
  fai_md5 = "c84f852b1cd4a00b12e6439ae7a2dd87"
) {
  # Preserve backward compatibility for misspelled argument name.
  if (!missing(asm_verison)) {
    asm_version <- asm_verison
  }

  if (
    !is.character(asm_version) ||
      length(asm_version) != 1L ||
      is.na(asm_version)
  ) {
    stop(
      "asm_version must be a single non-missing character string.",
      call. = FALSE
    )
  }

  fai_path <- tempfile()
  download.file(url = fai_url, destfile = fai_path, quiet = TRUE)
  if (!identical(unname(tools::md5sum(fai_path)), fai_md5)) {
    stop(glue::glue(
      "Downloaded hg002q100 {asm_version} fai file has incorrect md5sum"
    ))
  }

  fai_df <- readr::read_tsv(
    fai_path,
    col_names = c("chrom", "length", "offset", "line_bases", "line_width"),
    col_types = "cicii"
  ) |>
    dplyr::mutate(chrom = stringr::str_remove(chrom, "_.ATERNAL")) |>
    dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y"))) |>
    dplyr::select(chrom, length)
  file.remove(fai_path)

  asm_size <- sum(fai_df$length)

  if (is.na(asm_size) || asm_size <= 0) {
    stop("Error calculating hg002q100 size from fai file")
  }

  ## Return size in base pairs
  return(asm_size)
}

#' Load Full Variant Table For a Single Benchmarkset
#'
#' Loads full variant-level data from variants.tsv files. These are large
#' files (~GB per benchmark) and should be used only when variant-level
#' detail is required. Consider using `load_stratification_metrics()` for
#' aggregated summaries instead.
#'
#' @param benchmark_id Single benchmark identifier (e.g., "v5.0q_GRCh38_smvar")
#' @param results_dir Path to results directory. Default: `here::here("results")`
#' @param filters Optional list with filters:
#'   - `variant_types`: Character vector of variant types to include
#'   - `chromosomes`: Character vector of chromosomes to include
#'   - `in_benchmark_only`: Logical, keep only variants in benchmark regions (default: TRUE)
#'
#' @return Tibble with variant-level data:
#'   - chrom, pos, end: Genomic coordinates
#'   - gt: Genotype
#'   - vkx: Variant class
#'   - var_type: Variant type (SNP, INDEL, DEL, INS, COMPLEX, OTHER)
#'   - len_ref, len_alt: Reference and alternate allele lengths
#'   - var_size: Size of variant (len_alt - len_ref for small variants)
#'   - region_ids: Region classification
#'   - Additional columns from original file
#'
#' @examples
#' \dontrun{
#' # Load with default settings
#' variants <- load_variant_table("v5.0q_GRCh38_smvar")
#'
#' # Load with filters
#' variants <- load_variant_table(
#'   "v5.0q_GRCh38_smvar",
#'   filters = list(
#'     chromosomes = c("chr1", "chr2"),
#'     variant_types = c("SNP"),
#'     in_benchmark_only = TRUE
#'   )
#' )
#' }
#'
#' @param use_cache Logical; if TRUE (default), use Parquet caching
#' @param force_refresh Logical; if TRUE, ignore existing cache and reload
#'
#' @export
load_variant_table <- function(
  benchmark_id,
  results_dir = NULL,
  filters = NULL
) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Parse benchmark ID to construct path
  meta <- parse_benchmark_id(benchmark_id)

  variant_file <- fs::path(
    results_dir,
    "variant_tables",
    glue::glue("{meta$bench_version}_{meta$ref}_{meta$bench_type}"),
    "variants.tsv"
  )

  if (!fs::file_exists(variant_file)) {
    stop(
      glue::glue(
        "Variant file not found at: {variant_file}\n",
        "Check that benchmark ID is correct."
      ),
      call. = FALSE
    )
  }

  message(
    glue::glue(
      "Loading variant table for {benchmark_id}.\n",
      "This may take several minutes as the file is large (~GB)."
    )
  )

  # Read variant table
  variants_df <- read_variant_table(variant_file)

  # Apply filters if provided
  if (!is.null(filters)) {
    if (!is.null(filters$in_benchmark_only) && filters$in_benchmark_only) {
      variants_df <- variants_df %>%
        dplyr::filter(stringr::str_detect(region_ids, "^BMKREGIONS"))
    }

    if (!is.null(filters$variant_types)) {
      variants_df <- variants_df %>%
        dplyr::filter(var_type %in% filters$variant_types)
    }

    if (!is.null(filters$chromosomes)) {
      variants_df <- variants_df %>%
        dplyr::filter(chrom %in% filters$chromosomes)
    }
  }
  return(variants_df)
}

get_bench_var_cols <- function(benchmarkset) {
  ## Common column names and types
  cnames <- c(
    "chrom",
    "pos",
    "end",
    "gt",
    "vkx",
    "var_type",
    "len_ref",
    "len_alt",
    "context_ids",
    "region_ids"
  )
  ctypes <- "ciiccciicc"

  if (str_detect(benchmarkset, "v0.6")) {
    cnames <- c(
      cnames,
      "END",
      "SVTYPE",
      "SVLEN",
      "ClusterIDs",
      "NumClusterSVs",
      "ExactMatchIDs",
      "NumExactMatchSVs",
      "ClusterMaxShiftDist",
      "ClusterMaxSizeDiff",
      "ClusterMaxEditDist",
      "PBcalls",
      "Illcalls",
      "TenXcalls",
      "CGcalls",
      "PBexactcalls",
      "Illexactcalls",
      "TenXexactcalls",
      "CGexactcalls",
      "HG2count",
      "HG3count",
      "HG4count",
      "NumTechs",
      "NumTechsExact",
      "DistBack",
      "DistForward",
      "DistMin",
      "DistMinlt1000",
      "MultiTech",
      "MultiTechExact",
      "sizecat",
      "DistPASSHG2gt49Minlt1000",
      "DistPASSMinlt1000",
      "MendelianError",
      "HG003_GT",
      "HG004_GT",
      "BREAKSIMLENGTH",
      "REFWIDENED",
      "REPTYPE",
      "TRall",
      "TRgt100",
      "TRgt10k",
      "segdup"
    )
    ctypes <- paste0(ctypes, "ccicicccccccccccccccccccccccccccccc")
    return(list(col_names = cnames, col_types = ctypes))
  } else if (str_detect(benchmarkset, "v4.2.1")) {
    cnames <- c(
      cnames,
      "DPSum",
      "platforms",
      "platformnames",
      "platformbias",
      "datasets",
      "datasetnames",
      "datasetsmissingcall",
      "callsets",
      "callsetnames",
      "varType",
      "filt",
      "callable",
      "difficultregion",
      "arbitrated",
      "callsetwiththisuniqgenopassing",
      "callsetwithotheruniqgenopassing"
    )
    ctypes <- paste0(ctypes, "icccccccccccccicccccc")
    return(list(col_names = cnames, col_types = ctypes))
  } else if (str_detect(benchmarkset, "v5.0q")) {
    cnames <- c(
      cnames,
      "TRF",
      "TRFdiff",
      "TRFrepeat",
      "TRFovl",
      "TRFstart",
      "TRFend",
      "TRFperiod",
      "TRFcopies",
      "TRFscore",
      "TRFentropy",
      "TRFsim"
    )
    ctypes <- paste0(ctypes, "cicdiiiinnn")
    if (str_detect(benchmarkset, "stvar")) {
      cnames <- c(
        cnames,
        "SVTYPE",
        "SVLEN",
        "RM_score",
        "RM_repeat",
        "RM_clsfam",
        "LCR",
        "REMAP"
      )
      ctypes <- paste0(ctypes, "cincccc")
    }
    return(list(col_names = cnames, col_types = ctypes))
  } else {
    stop(paste("Unknown benchmark set :", benchmarkset))
  }

  return(list(col_names = cnames, col_types = ctypes))
}

tidy_smvar <- function(var_df) {
  ## excluding variants less than 50bp
  lt_50bp <- var_df$len_ref < 50 & var_df$len_alt < 50
  var_df <- var_df[lt_50bp, ]

  ## changing var_type for OTHER and OVERLAP
  var_df <- var_df %>%
    mutate(
      var_type = case_when(
        ## When REF is different from the first base of ALT
        var_type == "OTHER" &
          (len_ref > 1 & len_alt == 1) |
          (len_ref == 1 & len_alt > 1) ~ "INDEL",
        ## Assigning variant types for overlapping (atmoic) variants
        var_type == "OVERLAP" & len_ref == 1 & len_alt == 1 ~ "SNP",
        var_type == "OVERLAP" &
          (len_ref > 1 & len_alt == 1) |
          (len_ref == 1 & len_alt > 1) ~ "INDEL",
        var_type == "OVERLAP" & len_ref > 1 & len_alt > 1 ~ "COMPLEX",
        TRUE ~ var_type
      )
    )

  ## Annotating Variant Length
  var_df$var_size <- var_df$len_alt - var_df$len_ref

  return(var_df)
}

tidy_stvar <- function(var_df) {
  ## excluding variants less than 50bp
  gt_50bp <- var_df$len_ref > 49 | var_df$len_alt > 49
  var_df <- var_df[gt_50bp, ]

  ## changing var_type for OTHER and OVERLAP
  var_df$var_type <- var_df$SVTYPE

  ## Annotating Variant Length
  var_df <- var_df %>%
    dplyr::mutate(
      var_size = case_when(
        SVTYPE == "INS" ~ SVLEN,
        SVTYPE == "DEL" & SVLEN < 0 ~ SVLEN,
        SVTYPE == "DEL" & SVLEN > 0 ~ -SVLEN,
        TRUE ~ 0
      )
    )

  return(var_df)
}


#' Read Variant Table
read_variant_table <- function(table_path) {
  ## Extracting benchmark set id from file path and metadata from benchmark_id
  benchmark_id <- str_extract(
    table_path,
    "(?<=variant_tables/).*(?=/variants.tsv)"
  )
  meta <- parse_benchmark_id(benchmark_id)

  ## Getting column info based on benchmarkset
  col_info <- get_bench_var_cols(benchmark_id)

  ## reading variant table with vroom and selecting relevant columns
  var_df <- vroom::vroom(
    table_path,
    delim = "\t",
    na = ".",
    skip = 1,
    col_names = col_info$col_names,
    col_types = col_info$col_types,
    progress = TRUE
  ) %>%
    dplyr::mutate(
      chrom = std_chrom(chrom)
    )

  ## remove variants outside benchmark regions
  in_bmk <- grepl(x = var_df$region_ids, pattern = "^BMKREGIONS")
  var_df <- var_df[in_bmk, ]

  ## Cleaning up variant tables
  if (meta$bench_type == "smvar") {
    print("Tidying small variant table")
    var_df <- tidy_smvar(var_df)
  } else if (meta$bench_type == "stvar") {
    print("Tidying SV table")
    var_df <- tidy_stvar(var_df)
  } else {
    print("Dirty Table!! Fix code")
  }

  ## Adding benchmark metadata columns
  if (!is.null(benchmark_id)) {
    var_df <- var_df %>%
      tibble::add_column(
        bench_version = meta$bench_version,
        ref = meta$ref,
        bench_type = meta$bench_type,
        .before = 1
      )
    var_df
  }

  ## Reducing number of columns
  var_df %>%
    dplyr::select(
      bench_version,
      ref,
      bench_type,
      chrom,
      pos,
      end,
      gt,
      vkx,
      var_type,
      var_size,
      len_ref,
      len_alt,
      region_ids,
      context_ids
    )
}

#' Load All Variant Tables
#'
#' Loads and combines all variant tables from the results directory into a single tibble.
#' This is a large operation and should be used with caution.
#' Consider using smaller input data with pre-aggregated metrics for most analyses.
#'
load_all_variant_tables <- function(
  results_dir = NULL,
  use_cache = TRUE,
  force_refresh = FALSE
) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Find all variant tables
  variant_files <- fs::dir_ls(
    fs::path(results_dir, "variant_tables"),
    recurse = TRUE,
    glob = "**/variants.tsv",
    fail = FALSE
  )

  if (length(variant_files) == 0) {
    stop(
      glue::glue(
        "No variant tables found in {results_dir}/variant_tables"
      ),
      call. = FALSE
    )
  }

  # Try cache first
  source_files <- as.character(variant_files)

  if (use_cache && !force_refresh) {
    cached <- read_cache("variant_table", source_files)
    if (!is.null(cached)) {
      return(cached)
    }
  }

  require(furrr)
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  ## Loading and combining large variant tables
  future::plan(future::multisession, workers = parallel::detectCores() - 1)

  variants_df <- furrr::future_map_dfr(
    variant_files,
    read_variant_table,
    .progress = TRUE
  ) %>%
    dplyr::mutate(
      bench_version = std_bench_versions(bench_version),
      ref = std_references(ref),
      bench_type = std_bench_types(bench_type),
      .before = 1
    )

  # Write to cache
  if (use_cache) {
    tryCatch(
      write_cache(variants_df, "variant_table", source_files),
      error = function(e) {
        warning(
          glue::glue("Failed to write cache for variant_table: {e$message}"),
          call. = FALSE
        )
      }
    )
  }

  return(variants_df)
}

#' Load Genomic Context Coverage
#'
#' Loads base-level coverage data for difficult genomic context generated by
#' bedtools coverage. These files are in BED format with overlap counts.
#'
#' @param benchmark_id Single benchmark identifier (e.g., "v5.0q_GRCh38_smvar")
#' @param results_dir Path to results directory. Default: `here::here("results")`
#' @param context_filter Optional character vector of genomic context names to include
#'   (e.g., c("HP", "TR")). Valid values: HP, MAP, SD, SD10kb, TR, TR10kb
#'
#' @return Tibble with columns:
#'   - context_name: Genomic context region name
#'   - chrom: Chromosome
#'   - start: Start position (0-based)
#'   - end: End position
#'   - n_overlap: Number of overlapping intervals (from bedtools)
#'   - bases_cov: Number of bases covered
#'   - ivl_len: Total interval length
#'   - frac_cov: Fraction of interval covered
#'
#' @examples
#' \dontrun{
#' # Load all genomic contexts
#' coverage <- load_genomic_context_coverage("v5.0q_GRCh38_smvar")
#'
#' # Load specific genomic contexts
#' coverage <- load_diff_coverage(
#'   "v5.0q_GRCh38_smvar",
#'   context_filter = c("HP", "TR")
#' )
#' }
#'
#' @param use_cache Logical; if TRUE (default), use Parquet caching
#' @param force_refresh Logical; if TRUE, ignore existing cache and reload
#'
#' @export
load_genomic_context_coverage <- function(
  benchmark_id,
  results_dir = NULL,
  context_filter = NULL,
  use_cache = TRUE,
  force_refresh = FALSE
) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Parse benchmark ID to construct path
  meta <- parse_benchmark_id(benchmark_id)

  coverage_dir <- fs::path(
    results_dir,
    "diff_region_coverage",
    glue::glue("{meta$bench_version}_{meta$ref}_{meta$bench_type}")
  )

  if (!fs::dir_exists(coverage_dir)) {
    stop(
      glue::glue(
        "Coverage directory not found at: {coverage_dir}\n",
        "Check that benchmark ID is correct."
      ),
      call. = FALSE
    )
  }

  # Find all coverage files
  coverage_files <- fs::dir_ls(
    coverage_dir,
    glob = "*_cov.bed"
  )

  if (length(coverage_files) == 0) {
    stop(
      glue::glue(
        "No coverage files found in {coverage_dir}"
      ),
      call. = FALSE
    )
  }

  # Try cache first
  source_files <- as.character(coverage_files)
  cache_params <- list(
    benchmark_id = benchmark_id,
    context_filter = context_filter
  )

  if (use_cache && !force_refresh) {
    cached <- read_cache("diff_coverage", source_files, cache_params)
    if (!is.null(cached)) {
      return(cached)
    }
  }

  # Read all coverage files
  coverage_df <- coverage_files %>%
    purrr::map_dfr(function(file) {
      # Extract genomic context name from filename
      context <- fs::path_file(file) %>%
        stringr::str_remove("_cov.bed")

      file %>%
        vroom::vroom(
          col_names = c(
            "chrom",
            "start",
            "end",
            "n_overlap",
            "bases_cov",
            "ivl_len"
          ),
          col_types = "ciiiii",
          show_col_types = FALSE
        ) %>%
        tibble::add_column(context_name = context, .before = 1) %>%
        dplyr::mutate(
          frac_cov = bases_cov / ivl_len
        )
    })

  # Apply genomic context filter if provided
  if (!is.null(context_filter)) {
    coverage_df <- coverage_df %>%
      dplyr::filter(context_name %in% context_filter)

    if (nrow(coverage_df) == 0) {
      warning(
        glue::glue(
          "No rows remaining after applying genomic context filter. ",
          "Check that genomic context names are valid."
        )
      )
    }
  }

  # Write to cache
  if (use_cache) {
    tryCatch(
      write_cache(coverage_df, "diff_coverage", source_files, cache_params),
      error = function(e) {
        warning(
          glue::glue("Failed to write cache for diff_coverage: {e$message}"),
          call. = FALSE
        )
      }
    )
  }

  return(coverage_df)
}


#' Load Benchmark Region Files
#'
#' Loads all benchmark region BED files from the resources directory and combines them
#' into a single tibble with standardized column names and factored variables.
#'
#' @param resources_dir Path to resources directory. Default: `here::here("resources/benchmarksets")`
#'
#' @return Tibble with columns:
#'   - bench_version: Benchmark version (factored)
#'   - ref: Reference genome (factored)
#'   - bench_type: Variant type classification (factored: smvar, stvar)
#'   - chrom: Chromosome name with "chr" prefix (factored)
#'   - start: Start position (0-based)
#'   - end: End position (1-based)
#'   - interval_size: Size of region in bases (end - start)
#'
#' @examples
#' \dontrun{
#' bench_regions <- load_benchmark_regions()
#'
#' # Get regions for specific benchmark
#' bench_regions %>%
#'   filter(bench_version == "v5.0q", ref == "GRCh38")
#' }
#'
#' @param use_cache Logical; if TRUE (default), use Parquet caching
#' @param force_refresh Logical; if TRUE, ignore existing cache and reload
#'
#' @export
load_benchmark_regions <- function(
  resources_dir = NULL,
  use_cache = TRUE,
  force_refresh = FALSE
) {
  if (is.null(resources_dir)) {
    resources_dir <- here::here("resources/benchmarksets")
  }

  if (!fs::dir_exists(resources_dir)) {
    stop(
      glue::glue(
        "Resources directory not found: {resources_dir}"
      ),
      call. = FALSE
    )
  }

  # Find all benchmark BED files
  bench_paths <- fs::dir_ls(
    resources_dir,
    glob = "*_benchmark.bed",
    fail = FALSE
  )

  if (length(bench_paths) == 0) {
    stop(
      glue::glue(
        "No benchmark BED files found in {resources_dir}"
      ),
      call. = FALSE
    )
  }

  bench_ids <- stringr::str_remove(
    fs::path_file(bench_paths),
    "_benchmark\\.bed$"
  )
  bench_files <- stats::setNames(bench_paths, bench_ids)

  # Try cache first
  source_files <- as.character(bench_files)
  cache_params <- list(resources_dir = resources_dir)

  if (use_cache && !force_refresh) {
    cached <- read_cache(
      "benchmark_regions",
      source_files,
      cache_params,
      validate = TRUE
    )
    if (!is.null(cached)) {
      return(cached)
    }
  }

  bench_regions_df <- .read_region_bed_files(
    region_files = bench_files,
    bench_version_levels = bench_version,
    ref_levels = refs,
    bench_type_levels = bench_type
  )

  # Write to cache
  if (use_cache) {
    .write_cache_safely(
      data = bench_regions_df,
      dataset_name = "benchmark_regions",
      source_files = source_files,
      cache_params = cache_params
    )
  }

  return(bench_regions_df)
}

#' Load Platinum Pedigree Region Files
#'
#' Downloads Platinum Pedigree v1.2 benchmark BED files from public S3 and
#' loads them into a standardized tibble.
#'
#' @param download_dir Path where the Platinum Pedigree files should be stored.
#'   Default: `here::here("resources/platinum-pedigree-data/truthset_v1.2")`
#' @param s3_uri Public S3 URI to sync files from.
#'   Default: `"s3://platinum-pedigree-data/truthset_v1.2/"`
#' @param use_cache Logical; if TRUE (default), use Parquet caching
#' @param force_refresh Logical; if TRUE, ignore existing cache and reload
#' @param sync_s3 Logical; if TRUE (default), run `aws s3 sync --no-sign-request`
#'   when files are missing locally (or when `force_refresh = TRUE`)
#'
#' @return Tibble with columns:
#'   - bench_version: Benchmark version (`PP`)
#'   - ref: Reference genome (`GRCh38`)
#'   - bench_type: Variant type classification (`smvar`, `stvar`)
#'   - chrom: Chromosome name with "chr" prefix
#'   - start: Start position (0-based)
#'   - end: End position (1-based)
#'   - interval_size: Size of region in bases (end - start)
#'
#' @examples
#' \dontrun{
#' pp_regions <- load_platinum_pedigree_regions()
#' pp_regions %>% dplyr::count(bench_type)
#' }
#'
#' @export
load_platinum_pedigree_regions <- function(
  download_dir = NULL,
  s3_uri = "s3://platinum-pedigree-data/truthset_v1.2/",
  use_cache = TRUE,
  force_refresh = FALSE,
  sync_s3 = TRUE
) {
  if (is.null(download_dir)) {
    download_dir <- here::here("resources/platinum-pedigree-data/truthset_v1.2")
  }

  fs::dir_create(download_dir, recurse = TRUE)

  pp_region_list <- c(
    PP_GRCh38_smvar = fs::path(download_dir, "NA12878_hq_v1.2.smallvar.bed.gz"),
    PP_GRCh38_stvar = fs::path(download_dir, "NA12878_hq_v1.2.svs.bed.gz")
  )

  source_files <- as.character(pp_region_list)
  cache_params <- list(download_dir = download_dir, s3_uri = s3_uri)

  if (use_cache && !force_refresh) {
    cached <- read_cache(
      "platinum_pedigree_regions",
      source_files,
      cache_params,
      validate = TRUE
    )
    if (!is.null(cached)) {
      return(cached)
    }
  }

  missing_files <- pp_region_list[!fs::file_exists(pp_region_list)]
  should_sync <- length(missing_files) > 0 || (force_refresh && sync_s3)

  if (should_sync) {
    if (!sync_s3) {
      stop(
        glue::glue(
          "Missing Platinum Pedigree BED files in {download_dir}: ",
          "{paste(names(missing_files), collapse = ', ')}. ",
          "Set sync_s3 = TRUE to download from S3."
        ),
        call. = FALSE
      )
    }

    aws_bin <- Sys.which("aws")
    if (identical(aws_bin, "")) {
      stop(
        glue::glue(
          "AWS CLI not found in PATH. Install AWS CLI or place files in {download_dir}."
        ),
        call. = FALSE
      )
    }

    sync_args <- c("s3", "sync", "--no-sign-request", s3_uri, download_dir)
    sync_output <- system2(
      aws_bin,
      args = sync_args,
      stdout = TRUE,
      stderr = TRUE
    )
    sync_status <- attr(sync_output, "status")

    if (!is.null(sync_status) && sync_status != 0) {
      stop(
        glue::glue(
          "Failed to sync Platinum Pedigree files from S3 (exit {sync_status}).\n",
          "{paste(sync_output, collapse = '\n')}"
        ),
        call. = FALSE
      )
    }
  }

  missing_files <- pp_region_list[!fs::file_exists(pp_region_list)]
  if (length(missing_files) > 0) {
    stop(
      glue::glue(
        "Missing Platinum Pedigree BED files in {download_dir}: ",
        "{paste(names(missing_files), collapse = ', ')}"
      ),
      call. = FALSE
    )
  }

  pp_regions_df <- .read_region_bed_files(
    region_files = pp_region_list,
    bench_version_levels = c("PP"),
    ref_levels = c("GRCh38"),
    bench_type_levels = bench_type
  )

  if (use_cache) {
    .write_cache_safely(
      data = pp_regions_df,
      dataset_name = "platinum_pedigree_regions",
      source_files = source_files,
      cache_params = cache_params
    )
  }

  return(pp_regions_df)
}

#' Backward Compatibility Aliases
#'
#' These functions maintain backward compatibility with code using old terminology
#' ("stratification" instead of "genomic context").
#'
#' @keywords internal
#' @name stratification-aliases

#' @rdname stratification-aliases
#' @export
load_stratification_metrics <- function(
  results_dir = NULL,
  benchmark_filter = NULL
) {
  load_genomic_context_metrics(
    results_dir = results_dir,
    benchmark_filter = benchmark_filter
  )
}
