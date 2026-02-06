#' Parse Benchmark Identifier from File Path
#'
#' Extracts benchmark metadata (version, reference, variant type) from a file path
#' or benchmark ID string using the standard naming pattern: v{version}_{ref}_{var_type}
#'
#' @param file_path Character string containing file path or benchmark ID
#'
#' @return Named list with elements:
#'   - bench_version: Benchmark version (e.g., "v5.0q")
#'   - ref: Reference genome (e.g., "GRCh38")
#'   - var_type: Variant type ("smvar" for small variants or "stvar" for structural variants)
#'
#' @examples
#' \dontrun{
#' parse_benchmark_id("results/var_counts/v5.0q_GRCh38_smvar/stratification_combined_metrics.csv")
#' # Returns list(bench_version = "v5.0q", ref = "GRCh38", var_type = "smvar")
#' }
#'
#' @export
parse_benchmark_id <- function(file_path) {
  # Extract benchmark ID pattern: v[number]q?_[alphanumeric]_[smvar|stvar]
  pattern <- "v([0-9.]+q?)_([A-Za-z0-9.]+)_(smvar|stvar)"

  matches <- stringr::str_extract(file_path, pattern)

  if (is.na(matches)) {
    stop(
      glue::glue(
        "Could not parse benchmark ID from: {file_path}\n",
        "Expected pattern: v[version]_[reference]_[smvar|stvar]"
      ),
      call. = FALSE
    )
  }

  # Extract components
  components <- stringr::str_match(matches, pattern)[1, ]

  list(
    bench_version = components[2],
    ref = components[3],
    var_type = components[4]
  )
}


#' Load Stratification Metrics
#'
#' Loads primary analysis data files containing per-stratification metrics and variant counts.
#' These are the smallest, fastest-loading files and should be used for most analyses.
#'
#' @param results_dir Path to results directory. Default: `here::here("results")`
#' @param benchmark_filter Optional character vector of benchmark IDs to filter results
#'
#' @return Tibble with columns:
#'   - bench_version: Benchmark version
#'   - ref: Reference genome
#'   - var_type: Variant type (smvar or stvar)
#'   - strat_name: Stratification region name (HP, MAP, SD, SD10kb, TR, TR10kb)
#'   - strat_bp: Total bases in stratification
#'   - intersect_bp: Bases overlapping with benchmark
#'   - pct_of_strat: Percentage of stratification overlapping benchmark
#'   - pct_of_bench: Percentage of benchmark overlapping stratification
#'   - total_variants: Total variant count
#'   - snp_count: SNP count
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
#' metrics <- load_stratification_metrics()
#'
#' # Load specific benchmarks
#' metrics <- load_stratification_metrics(
#'   benchmark_filter = c("v5.0q_GRCh38_smvar", "v5.0q_GRCh37_stvar")
#' )
#' }
#'
#' @export
load_stratification_metrics <- function(results_dir = NULL, benchmark_filter = NULL) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Find all stratification metrics files
  metrics_files <- fs::dir_ls(
    results_dir,
    recurse = TRUE,
    glob = "**/stratification_combined_metrics.csv"
  )

  if (length(metrics_files) == 0) {
    stop(
      glue::glue(
        "No stratification_combined_metrics.csv files found in {results_dir}"
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
            strat_name = "c",
            strat_bp = "d",
            intersect_bp = "d",
            pct_of_strat = "d",
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
          var_type = meta$var_type,
          .before = 1
        )
    })

  # Apply benchmark filter if provided
  if (!is.null(benchmark_filter)) {
    metrics_df <- metrics_df %>%
      dplyr::filter(
        paste(bench_version, ref, var_type, sep = "_") %in% benchmark_filter
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
#'   - var_type: Variant type
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
          var_type = meta$var_type,
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
          col_names = c("chrom", "length", "ns"),
          col_types = "cii",
          show_col_types = FALSE
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


#' Load Full Variant Table
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
#' @export
load_variant_table <- function(benchmark_id, results_dir = NULL, filters = NULL) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Parse benchmark ID to construct path
  meta <- parse_benchmark_id(benchmark_id)

  variant_file <- fs::path(
    results_dir,
    "variant_tables",
    glue::glue("{meta$bench_version}_{meta$ref}_{meta$var_type}"),
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
  variants_df <- vroom::vroom(
    variant_file,
    show_col_types = FALSE
  )

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


#' Load Difficult Region Coverage
#'
#' Loads base-level coverage data for difficult regions generated by
#' bedtools coverage. These files are in BED format with overlap counts.
#'
#' @param benchmark_id Single benchmark identifier (e.g., "v5.0q_GRCh38_smvar")
#' @param results_dir Path to results directory. Default: `here::here("results")`
#' @param strat_filter Optional character vector of stratification names to include
#'   (e.g., c("HP", "TR")). Valid values: HP, MAP, SD, SD10kb, TR, TR10kb
#'
#' @return Tibble with columns:
#'   - strat_name: Stratification region name
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
#' # Load all stratifications
#' coverage <- load_diff_coverage("v5.0q_GRCh38_smvar")
#'
#' # Load specific stratifications
#' coverage <- load_diff_coverage(
#'   "v5.0q_GRCh38_smvar",
#'   strat_filter = c("HP", "TR")
#' )
#' }
#'
#' @export
load_diff_coverage <- function(benchmark_id, results_dir = NULL, strat_filter = NULL) {
  if (is.null(results_dir)) {
    results_dir <- here::here("results")
  }

  # Parse benchmark ID to construct path
  meta <- parse_benchmark_id(benchmark_id)

  coverage_dir <- fs::path(
    results_dir,
    "diff_region_coverage",
    glue::glue("{meta$bench_version}_{meta$ref}_{meta$var_type}")
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

  # Read all coverage files
  coverage_df <- coverage_files %>%
    purrr::map_dfr(function(file) {
      # Extract stratification name from filename
      strat <- fs::path_file(file) %>%
        stringr::str_remove("_cov.bed")

      file %>%
        vroom::vroom(
          col_names = c("chrom", "start", "end", "n_overlap", "bases_cov", "ivl_len"),
          col_types = "ciiiii",
          show_col_types = FALSE
        ) %>%
        tibble::add_column(strat_name = strat, .before = 1) %>%
        dplyr::mutate(
          frac_cov = bases_cov / ivl_len
        )
    })

  # Apply stratification filter if provided
  if (!is.null(strat_filter)) {
    coverage_df <- coverage_df %>%
      dplyr::filter(strat_name %in% strat_filter)

    if (nrow(coverage_df) == 0) {
      warning(
        glue::glue(
          "No rows remaining after applying stratification filter. ",
          "Check that stratification names are valid."
        )
      )
    }
  }

  return(coverage_df)
}
