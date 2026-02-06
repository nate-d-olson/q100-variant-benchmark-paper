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
#' parse_benchmark_id("results/var_counts/v5.0q_GRCh38_smvar/genomic_context_combined_metrics.csv")
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
#'   - var_types: Vector of variant types (smvar, stvar)
#'
#' @examples
#' \dontrun{
#' config <- parse_pipeline_config()
#' config$num_benchmarks  # Get number of benchmarks
#' config$benchmarks      # Get all benchmark IDs
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
    purrr::map_df(~ {
      meta <- parse_benchmark_id(.x)
      tibble::tibble(
        benchmark_id = .x,
        bench_version = meta$bench_version,
        ref = meta$ref,
        var_type = meta$var_type
      )
    })

  list(
    benchmarks = benchmarks,
    num_benchmarks = length(benchmarks),
    benchmark_meta = benchmark_meta,
    references = unique(benchmark_meta$ref),
    var_types = unique(benchmark_meta$var_type),
    num_contexts = 6  # HP, MAP, SD, SD10kb, TR, TR10kb
  )
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
#'   - var_type: Variant type (factored: smvar, stvar)
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
load_genomic_context_metrics <- function(results_dir = NULL, benchmark_filter = NULL) {
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

  # Rename snp_count to snv_count (SNV = Single Nucleotide Variant)
  # and factor categorical variables with standard levels
  metrics_df <- metrics_df %>%
    dplyr::rename(snv_count = snp_count) %>%
    dplyr::mutate(
      bench_version = factor(
        bench_version,
        levels = c("v0.6", "v4.2.1", "v5.0q")
      ),
      ref = factor(
        ref,
        levels = c("GRCh37", "GRCh38", "CHM13v2.0")
      ),
      var_type = factor(
        var_type,
        levels = c("smvar", "stvar")
      ),
      context_name = factor(
        context_name,
        levels = c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
      )
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
#' coverage <- load_diff_coverage("v5.0q_GRCh38_smvar")
#'
#' # Load specific genomic contexts
#' coverage <- load_diff_coverage(
#'   "v5.0q_GRCh38_smvar",
#'   context_filter = c("HP", "TR")
#' )
#' }
#'
#' @export
load_diff_coverage <- function(benchmark_id, results_dir = NULL, context_filter = NULL) {
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
      # Extract genomic context name from filename
      context <- fs::path_file(file) %>%
        stringr::str_remove("_cov.bed")

      file %>%
        vroom::vroom(
          col_names = c("chrom", "start", "end", "n_overlap", "bases_cov", "ivl_len"),
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
#' @export
load_benchmark_regions <- function(resources_dir = NULL) {
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
  bench_files <- fs::dir_ls(
    resources_dir,
    glob = "*_benchmark.bed",
    fail = FALSE
  )

  if (length(bench_files) == 0) {
    stop(
      glue::glue(
        "No benchmark BED files found in {resources_dir}"
      ),
      call. = FALSE
    )
  }

  # Read all benchmark region files
  bench_regions_df <- bench_files %>%
    purrr::map_dfr(
      ~ vroom::vroom(
        .x,
        col_names = c("chrom", "start", "end"),
        col_types = "cii",
        show_col_types = FALSE
      ) %>%
        tibble::add_column(
          benchmark_file = fs::path_file(.x),
          .before = 1
        ),
      .id = "file_id"
    ) %>%
    # Extract benchmark identifier from filename
    dplyr::mutate(
      benchmark_id = stringr::str_remove(benchmark_file, "_benchmark.bed"),
      bench_meta = purrr::map(benchmark_id, parse_benchmark_id),
      bench_version = purrr::map_chr(bench_meta, "bench_version"),
      ref = purrr::map_chr(bench_meta, "ref"),
      bench_type = purrr::map_chr(bench_meta, "var_type")
    ) %>%
    # Standardize chromosome names and calculate interval size
    dplyr::mutate(
      chrom = dplyr::if_else(
        stringr::str_detect(chrom, "^chr"),
        chrom,
        stringr::str_c("chr", chrom)
      ),
      interval_size = end - start,
      # Factor variables with standard levels
      bench_version = factor(
        bench_version,
        levels = c("v0.6", "v4.2.1", "v5.0q")
      ),
      ref = factor(
        ref,
        levels = c("GRCh37", "GRCh38", "CHM13v2.0")
      ),
      bench_type = factor(
        bench_type,
        levels = c("smvar", "stvar")
      ),
      chrom = factor(
        chrom,
        levels = stringr::str_c("chr", c(1:22, "X", "Y"))
      )
    ) %>%
    dplyr::select(
      -file_id,
      -benchmark_file,
      -benchmark_id,
      -bench_meta
    )

  return(bench_regions_df)
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
load_stratification_metrics <- function(results_dir = NULL, benchmark_filter = NULL) {
  load_genomic_context_metrics(results_dir = results_dir, benchmark_filter = benchmark_filter)
}
