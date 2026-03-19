# BED interval arithmetic and variant subsetting helpers
# Used by benchmark_unique_regions.qmd for version comparison analysis
#
# AI Disclosure: Developed with assistance from Claude (Anthropic).

library(data.table)

empty_bed <- function() {
  data.table(
    chrom = character(),
    start = integer(),
    end = integer()
  )
}

normalize_chrom <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

resolve_data_root <- function() {
  current_root <- here::here()
  if (
    dir.exists(file.path(current_root, "results")) &&
      dir.exists(file.path(current_root, "resources"))
  ) {
    return(current_root)
  }

  git_common <- tryCatch(
    system2(
      "git",
      c("rev-parse", "--git-common-dir"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  if (length(git_common) > 0) {
    common_path <- git_common[[1]]
    common_abs <- if (grepl("^/", common_path)) {
      common_path
    } else {
      normalizePath(
        file.path(current_root, common_path),
        mustWork = FALSE
      )
    }
    candidate_root <- dirname(common_abs)
    if (
      dir.exists(file.path(candidate_root, "results")) &&
        dir.exists(file.path(candidate_root, "resources"))
    ) {
      return(candidate_root)
    }
  }

  stop(
    "Could not locate a data root containing both 'results/' and 'resources/'.",
    call. = FALSE
  )
}

merge_intervals <- function(dt) {
  if (nrow(dt) == 0) {
    return(empty_bed())
  }
  x <- copy(dt)[, .(chrom, start, end)]
  x <- x[end > start]
  if (nrow(x) == 0) {
    return(empty_bed())
  }
  setorder(x, chrom, start, end)
  x[, cum_end := cummax(end), by = chrom]
  x[,
    grp_break := ifelse(
      is.na(shift(cum_end)),
      1L,
      as.integer(start > shift(cum_end))
    ),
    by = chrom
  ]
  x[, grp := cumsum(grp_break), by = chrom]
  merged <- x[,
    .(
      start = min(start),
      end = max(end)
    ),
    by = .(chrom, grp)
  ][, grp := NULL]
  setorder(merged, chrom, start, end)
  merged[]
}

read_bed <- function(path, merge = TRUE) {
  if (!file.exists(path)) {
    stop(glue::glue("BED file not found: {path}"), call. = FALSE)
  }

  dt <- if (grepl("\\.gz$", path)) {
    data.table::fread(
      cmd = paste("gzip -dc", shQuote(path)),
      select = 1:3,
      col.names = c("chrom", "start", "end"),
      showProgress = FALSE
    )
  } else {
    data.table::fread(
      path,
      select = 1:3,
      col.names = c("chrom", "start", "end"),
      showProgress = FALSE
    )
  }

  dt[, chrom := normalize_chrom(chrom)]
  dt[, `:=`(
    start = as.integer(start),
    end = as.integer(end)
  )]
  dt <- dt[end > start]

  if (merge) merge_intervals(dt) else dt
}

interval_bases <- function(dt) {
  if (nrow(dt) == 0) {
    return(0)
  }
  sum(as.numeric(dt$end - dt$start))
}

intersect_intervals <- function(a, b) {
  if (nrow(a) == 0 || nrow(b) == 0) {
    return(empty_bed())
  }

  a <- merge_intervals(a)
  b <- merge_intervals(b)
  shared_chroms <- intersect(unique(a$chrom), unique(b$chrom))
  if (length(shared_chroms) == 0) {
    return(empty_bed())
  }

  out <- vector("list", length(shared_chroms))

  for (idx in seq_along(shared_chroms)) {
    chrom_id <- shared_chroms[[idx]]
    ad <- a[chrom == chrom_id]
    bd <- b[chrom == chrom_id]
    i <- 1L
    j <- 1L
    starts <- integer()
    ends <- integer()

    while (i <= nrow(ad) && j <= nrow(bd)) {
      s <- max(ad$start[[i]], bd$start[[j]])
      e <- min(ad$end[[i]], bd$end[[j]])
      if (s < e) {
        starts <- c(starts, s)
        ends <- c(ends, e)
      }
      if (ad$end[[i]] <= bd$end[[j]]) {
        i <- i + 1L
      } else {
        j <- j + 1L
      }
    }

    out[[idx]] <- if (length(starts) == 0) {
      NULL
    } else {
      data.table(chrom = chrom_id, start = starts, end = ends)
    }
  }

  merged <- rbindlist(out, use.names = TRUE, fill = TRUE)
  if (nrow(merged) == 0) {
    return(empty_bed())
  }
  merge_intervals(merged)
}

subtract_intervals <- function(a, b) {
  if (nrow(a) == 0) {
    return(empty_bed())
  }
  if (nrow(b) == 0) {
    return(merge_intervals(a))
  }

  a <- merge_intervals(a)
  b <- merge_intervals(b)
  chroms <- unique(a$chrom)
  out <- vector("list", length(chroms))

  for (idx in seq_along(chroms)) {
    chrom_id <- chroms[[idx]]
    ad <- a[chrom == chrom_id]
    bd <- b[chrom == chrom_id]

    if (nrow(bd) == 0) {
      out[[idx]] <- ad
      next
    }

    starts <- integer()
    ends <- integer()
    j <- 1L
    nb <- nrow(bd)

    for (i in seq_len(nrow(ad))) {
      s <- ad$start[[i]]
      e <- ad$end[[i]]

      while (j <= nb && bd$end[[j]] <= s) {
        j <- j + 1L
      }

      cur <- s
      k <- j

      while (k <= nb && bd$start[[k]] < e) {
        if (bd$start[[k]] > cur) {
          starts <- c(starts, cur)
          ends <- c(ends, min(bd$start[[k]], e))
        }
        cur <- max(cur, bd$end[[k]])
        if (cur >= e) {
          break
        }
        k <- k + 1L
      }

      if (cur < e) {
        starts <- c(starts, cur)
        ends <- c(ends, e)
      }
    }

    out[[idx]] <- if (length(starts) == 0) {
      NULL
    } else {
      data.table(chrom = chrom_id, start = starts, end = ends)
    }
  }

  merged <- rbindlist(out, use.names = TRUE, fill = TRUE)
  if (nrow(merged) == 0) {
    return(empty_bed())
  }
  merge_intervals(merged)
}

read_exclusion_beds <- function(benchmark_id, resources_dir) {
  exclusion_dir <- file.path(resources_dir, "exclusions", benchmark_id)
  if (!dir.exists(exclusion_dir)) {
    return(list())
  }

  files <- list.files(
    exclusion_dir,
    pattern = "\\.bed$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(list())
  }

  exclusion_names <- sub("_[0-9]+\\.bed$", "", basename(files))
  grouped <- split(files, exclusion_names)

  lapply(grouped, function(paths) {
    merged <- rbindlist(
      lapply(paths, function(p) read_bed(p, merge = FALSE)),
      use.names = TRUE
    )
    merge_intervals(merged)
  })
}

read_context_beds <- function(ref, resources_dir) {
  context_names <- c("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
  beds <- purrr::map(
    context_names,
    function(ctx) {
      read_bed(
        file.path(resources_dir, "stratifications", paste0(ref, "_", ctx, ".bed.gz"))
      )
    }
  )
  names(beds) <- context_names
  beds
}

subset_variants_by_regions <- function(variants_dt, regions_dt) {
  if (nrow(variants_dt) == 0 || nrow(regions_dt) == 0) {
    return(variants_dt[0])
  }

  reg <- copy(regions_dt)[, .(
    chrom,
    start,
    end = end - 1L
  )][end >= start]

  var_iv <- copy(variants_dt)[, .(
    idx = .I,
    chrom,
    start = pmax(as.integer(pos) - 1L, 0L),
    end = pmax(as.integer(pos), as.integer(end)) - 1L
  )][end >= start]

  setkey(reg, chrom, start, end)
  setkey(var_iv, chrom, start, end)
  hits <- foverlaps(var_iv, reg, nomatch = 0L)
  variants_dt[sort(unique(hits$idx))]
}

read_variant_table_parquet <- function(path, bench_type, benchmark_only = TRUE) {
  if (!file.exists(path)) {
    stop(glue::glue("Variant table not found: {path}"), call. = FALSE)
  }

  df <- arrow::read_parquet(path)
  dt <- as.data.table(df)

  dt[, chrom := normalize_chrom(chrom)]
  dt[, pos := as.integer(pos)]
  dt[, end := as.integer(end)]

  if (!"context_ids" %in% names(dt)) {
    dt[, context_ids := "."]
  }
  if (!"region_ids" %in% names(dt)) {
    dt[, region_ids := "."]
  }
  dt[, context_ids := as.character(context_ids)]
  dt[, region_ids := as.character(region_ids)]

  if (bench_type == "smvar") {
    dt[,
      variant_type := dplyr::case_when(
        var_type %in% c("SNP", "SNV") ~ "SNV",
        var_type == "INS" ~ "INS",
        var_type == "DEL" ~ "DEL",
        var_type == "INDEL" & alt_len > ref_len ~ "INS",
        var_type == "INDEL" & alt_len < ref_len ~ "DEL",
        var_type == "INDEL" ~ "INDEL",
        TRUE ~ "OTHER"
      )
    ]
    dt <- dt[abs(var_size) < 50]
  } else {
    dt[,
      variant_type := dplyr::case_when(
        var_type == "INS" ~ "SV_INS",
        var_type == "DEL" ~ "SV_DEL",
        alt_len > ref_len ~ "SV_INS",
        alt_len < ref_len ~ "SV_DEL",
        TRUE ~ "SV_OTHER"
      )
    ]
    dt <- dt[abs(var_size) >= 50]
  }

  dt[, in_benchmark := grepl("(^|,)BMKREGIONS(,|$)", region_ids)]
  if (benchmark_only) {
    dt <- dt[in_benchmark == TRUE]
  }

  dt[, .(
    chrom,
    pos,
    end,
    variant_type,
    var_size = abs(var_size),
    context_ids,
    region_ids
  )]
}

normalize_old_only_variants <- function(path, bench_type) {
  if (!file.exists(path)) {
    return(data.table())
  }
  dt <- fread(path, sep = "\t", showProgress = FALSE)
  if (nrow(dt) == 0) {
    return(dt)
  }

  dt[, chrom := normalize_chrom(chrom)]
  dt[, pos := as.integer(pos)]
  dt[, end := as.integer(end)]
  dt[, ref := as.character(ref)]
  dt[, alt := as.character(alt)]
  dt[, status := as.character(status)]
  dt[, exclusion_ids := as.character(exclusion_ids)]
  dt[is.na(exclusion_ids), exclusion_ids := ""]

  dt[, ref_len := nchar(ref)]
  dt[, alt_len := nchar(alt)]
  dt[, size_len := abs(alt_len - ref_len)]
  dt[, size_coord := pmax(end - pos + 1L, 0L)]
  dt[, var_size := pmax(size_len, size_coord, na.rm = TRUE)]

  dt[,
    variant_type := dplyr::case_when(
      var_type %in% c("SNP", "SNV") ~ "SNV",
      var_type == "INDEL" & alt_len > ref_len ~ "INS",
      var_type == "INDEL" & alt_len < ref_len ~ "DEL",
      var_type == "INDEL" ~ "INDEL",
      TRUE ~ "OTHER"
    )
  ]
  dt[variant_type == "SNV", var_size := 1L]

  if (bench_type == "smvar") {
    dt <- dt[var_size < 50]
  } else {
    dt <- dt[var_size >= 50]
  }

  dt[, .(
    chrom,
    pos,
    end,
    variant_type,
    var_size,
    status,
    exclusion_ids
  )]
}

size_bin_for <- function(var_size, variant_type, bench_type) {
  if (bench_type == "smvar") {
    dplyr::case_when(
      variant_type == "SNV" ~ "SNV (1 bp)",
      var_size <= 5 ~ "2-5 bp",
      var_size <= 14 ~ "6-14 bp",
      var_size <= 49 ~ "15-49 bp",
      TRUE ~ ">=50 bp"
    )
  } else {
    dplyr::case_when(
      var_size <= 99 ~ "50-99 bp",
      var_size <= 299 ~ "100-299 bp",
      var_size <= 999 ~ "300-999 bp",
      var_size <= 4999 ~ "1-4.9 kb",
      TRUE ~ ">=5 kb"
    )
  }
}

summarize_variant_type_size <- function(
  dt,
  bench_type,
  comparison_id,
  comparison_label,
  region_group
) {
  if (nrow(dt) == 0) {
    return(tibble::tibble(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      region_group = region_group,
      variant_type = character(),
      size_bin = character(),
      variant_count = integer()
    ))
  }

  tibble::as_tibble(dt) %>%
    dplyr::mutate(size_bin = size_bin_for(var_size, variant_type, bench_type)) %>%
    dplyr::count(variant_type, size_bin, name = "variant_count") %>%
    dplyr::mutate(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      region_group = region_group,
      .before = 1
    ) %>%
    dplyr::arrange(comparison_label, region_group, dplyr::desc(variant_count))
}

summarize_variant_context <- function(dt, comparison_id, comparison_label, region_group) {
  if (nrow(dt) == 0) {
    return(tibble::tibble(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      region_group = region_group,
      context_name = character(),
      variant_count = integer()
    ))
  }

  tibble::as_tibble(dt) %>%
    dplyr::transmute(context_ids) %>%
    dplyr::filter(!is.na(context_ids), context_ids != ".", context_ids != "") %>%
    tidyr::separate_rows(context_ids, sep = ",") %>%
    dplyr::filter(context_ids != "", !is.na(context_ids)) %>%
    dplyr::count(context_ids, name = "variant_count") %>%
    dplyr::rename(context_name = context_ids) %>%
    dplyr::mutate(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      region_group = region_group,
      .before = 1
    ) %>%
    dplyr::arrange(comparison_label, dplyr::desc(variant_count))
}

summarize_region_context <- function(
  regions_dt,
  context_beds,
  comparison_id,
  comparison_label,
  region_group
) {
  total_bp <- interval_bases(regions_dt)
  tibble::tibble(
    context_name = names(context_beds),
    overlap_bp = purrr::map_dbl(
      context_beds,
      ~ interval_bases(intersect_intervals(regions_dt, .x))
    )
  ) %>%
    dplyr::mutate(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      region_group = region_group,
      total_region_bp = total_bp,
      pct_of_region = if (total_bp > 0) {
        100 * overlap_bp / total_bp
      } else {
        0
      },
      .before = 1
    ) %>%
    dplyr::arrange(comparison_label, dplyr::desc(overlap_bp))
}

get_exclusion_descriptions <- function(config_path, benchmark_id) {
  cfg <- yaml::read_yaml(config_path)
  exclusions <- cfg$benchmarksets[[benchmark_id]]$exclusions
  if (is.null(exclusions) || length(exclusions) == 0) {
    return(tibble::tibble(exclusion = character(), description = character()))
  }
  tibble::tibble(
    exclusion = purrr::map_chr(exclusions, "name"),
    description = purrr::map_chr(exclusions, ~ .x$description %||% "")
  )
}

analyze_comparison <- function(comp, results_dir, resources_dir, config_path, context_cache) {
  old_bed <- read_bed(file.path(
    resources_dir,
    "benchmarksets",
    paste0(comp$old_benchmark, "_benchmark.bed")
  ))
  new_bed <- read_bed(file.path(
    resources_dir,
    "benchmarksets",
    paste0(comp$new_benchmark, "_benchmark.bed")
  ))
  new_dip_bed <- read_bed(file.path(
    resources_dir,
    "benchmarksets",
    paste0(comp$new_benchmark, "_dip.bed")
  ))

  old_only <- subtract_intervals(old_bed, new_bed)
  new_only <- subtract_intervals(new_bed, old_bed)
  shared <- intersect_intervals(old_bed, new_bed)

  exclusion_beds <- read_exclusion_beds(comp$new_benchmark, resources_dir)
  all_exclusions <- if (length(exclusion_beds) == 0) {
    empty_bed()
  } else {
    merge_intervals(rbindlist(exclusion_beds, use.names = TRUE))
  }

  old_only_in_dip <- intersect_intervals(old_only, new_dip_bed)
  old_only_excluded <- if (nrow(all_exclusions) > 0) {
    intersect_intervals(old_only_in_dip, all_exclusions)
  } else {
    empty_bed()
  }
  old_only_not_in_dip <- subtract_intervals(old_only, new_dip_bed)
  old_only_in_dip_not_excluded <- if (nrow(all_exclusions) > 0) {
    subtract_intervals(old_only_in_dip, all_exclusions)
  } else {
    old_only_in_dip
  }

  exclusion_desc <- get_exclusion_descriptions(config_path, comp$new_benchmark)

  old_exclusion_base <- purrr::imap_dfr(
    exclusion_beds,
    function(excl_bed, excl_name) {
      tibble::tibble(
        exclusion = excl_name,
        bases_bp = interval_bases(intersect_intervals(old_only_excluded, excl_bed))
      )
    }
  )

  old_excluded_bp <- interval_bases(old_only_excluded)

  old_exclusion_base <- old_exclusion_base %>%
    dplyr::left_join(exclusion_desc, by = "exclusion") %>%
    dplyr::mutate(
      comparison_id = comp$comp_id,
      comparison_label = comp$comparison_label,
      pct_of_excluded_bases = if (old_excluded_bp > 0) {
        100 * bases_bp / old_excluded_bp
      } else {
        0
      },
      .before = 1
    ) %>%
    dplyr::arrange(comparison_label, dplyr::desc(bases_bp))

  old_variant_dt <- normalize_old_only_variants(
    file.path(results_dir, "exclusions", comp$comp_id, "old_only_variants.tsv"),
    comp$bench_type
  )

  old_variant_status <- tibble::as_tibble(old_variant_dt) %>%
    dplyr::count(status, name = "variant_count") %>%
    dplyr::mutate(
      pct_of_old_only_variants = if (sum(variant_count) > 0) {
        100 * variant_count / sum(variant_count)
      } else {
        0
      }
    ) %>%
    dplyr::mutate(
      comparison_id = comp$comp_id,
      comparison_label = comp$comparison_label,
      .before = 1
    )

  old_variant_exclusion <- tibble::as_tibble(old_variant_dt) %>%
    dplyr::filter(status == "excluded", !is.na(exclusion_ids), exclusion_ids != "") %>%
    tidyr::separate_rows(exclusion_ids, sep = ",") %>%
    dplyr::count(exclusion_ids, name = "variant_count") %>%
    dplyr::rename(exclusion = exclusion_ids) %>%
    dplyr::left_join(exclusion_desc, by = "exclusion") %>%
    dplyr::mutate(
      comparison_id = comp$comp_id,
      comparison_label = comp$comparison_label,
      .before = 1
    ) %>%
    dplyr::arrange(comparison_label, dplyr::desc(variant_count))

  new_variant_dt <- read_variant_table_parquet(
    file.path(results_dir, "variant_tables", comp$new_benchmark, "variants.parquet"),
    bench_type = comp$bench_type,
    benchmark_only = TRUE
  )
  new_only_variant_dt <- subset_variants_by_regions(new_variant_dt, new_only)

  context_beds <- context_cache[[comp$ref]]

  summary_tbl <- tibble::tibble(
    comparison_id = comp$comp_id,
    comparison_label = comp$comparison_label,
    old_benchmark = comp$old_benchmark,
    new_benchmark = comp$new_benchmark,
    old_total_bp = interval_bases(old_bed),
    new_total_bp = interval_bases(new_bed),
    shared_bp = interval_bases(shared),
    old_only_bp = interval_bases(old_only),
    new_only_bp = interval_bases(new_only),
    old_only_variants = nrow(old_variant_dt),
    new_only_variants = nrow(new_only_variant_dt)
  ) %>%
    dplyr::mutate(
      pct_old_only_of_old = dplyr::if_else(old_total_bp > 0, 100 * old_only_bp / old_total_bp, 0),
      pct_new_only_of_v5 = dplyr::if_else(new_total_bp > 0, 100 * new_only_bp / new_total_bp, 0)
    )

  old_region_reason <- tibble::tibble(
    comparison_id = comp$comp_id,
    comparison_label = comp$comparison_label,
    category = c("not_in_dipbed", "excluded", "in_v5_dipbed_not_excluded"),
    bases_bp = c(
      interval_bases(old_only_not_in_dip),
      interval_bases(old_only_excluded),
      interval_bases(old_only_in_dip_not_excluded)
    )
  ) %>%
    dplyr::mutate(
      pct_of_old_only = if (sum(bases_bp) > 0) {
        100 * bases_bp / sum(bases_bp)
      } else {
        0
      }
    )

  list(
    summary = summary_tbl,
    old_region_reason = old_region_reason,
    old_exclusion_base = old_exclusion_base,
    old_variant_status = old_variant_status,
    old_variant_exclusion = old_variant_exclusion,
    old_variant_type_size = summarize_variant_type_size(
      old_variant_dt,
      comp$bench_type,
      comp$comp_id,
      comp$comparison_label,
      "old_only_regions"
    ),
    old_excluded_region_context = summarize_region_context(
      old_only_excluded,
      context_beds,
      comp$comp_id,
      comp$comparison_label,
      "old_only_excluded_regions"
    ),
    new_variant_type_size = summarize_variant_type_size(
      new_only_variant_dt,
      comp$bench_type,
      comp$comp_id,
      comp$comparison_label,
      "v5_only_regions"
    ),
    new_variant_context = summarize_variant_context(
      new_only_variant_dt,
      comp$comp_id,
      comp$comparison_label,
      "v5_only_regions"
    ),
    new_region_context = summarize_region_context(
      new_only,
      context_beds,
      comp$comp_id,
      comp$comparison_label,
      "v5_only_regions"
    )
  )
}
