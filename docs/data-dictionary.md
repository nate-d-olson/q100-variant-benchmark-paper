# Data Dictionary

This document is the single reference for all data structures in the Q100 variant benchmark project. It covers:

1. **R Loading Function Schemas** — column types and descriptions for data frames returned by `R/data_loading.R`
2. **Pipeline Metric Definitions** — detailed interpretations of metrics, genomic contexts, exclusions, and variant classifications

> **Maintenance Note:** Update this document whenever changes are made to `R/data_loading.R`, `R/schemas.R`, or pipeline output formats.

---

## R Loading Function Schemas

The functions in `R/data_loading.R` load and standardize pipeline outputs. This section documents the data frames they return. The pipeline output files (CSVs, TSVs, BEDs) are documented in [Pipeline Outputs Reference](pipeline-outputs.md).

### Genomic Context Metrics

**Function:** `load_genomic_context_metrics()`

| Column                   | Type    | Description                                                              |
| :----------------------- | :------ | :----------------------------------------------------------------------- |
| `bench_version`          | Factor  | Benchmark version (`v0.6`, `v4.2.1`, `v5.0q`)                            |
| `ref`                    | Factor  | Reference genome (`GRCh37`, `GRCh38`, `CHM13v2.0`)                       |
| `bench_type`             | Factor  | Benchmark set type (`smvar`, `stvar`)                                    |
| `context_name`           | Factor  | Genomic context identifier (e.g., `HP`, `TR`, `SD`)                      |
| `context_bp`             | Numeric | Total size of the genomic context in bases                               |
| `intersect_bp`           | Numeric | Number of bases where the genomic context overlaps the benchmark regions |
| `pct_of_context`         | Numeric | Percentage of the genomic context covered by the benchmark               |
| `pct_of_bench`           | Numeric | Percentage of the benchmark covered by this genomic context              |
| `total_variants`         | Integer | Total number of variants in this context                                 |
| `snv_count`              | Integer | Count of Single Nucleotide Variants (SNV)                                |
| `indel_count`            | Integer | Count of Insertions/Deletions (INDEL, smvar only)                        |
| `del_count`              | Integer | Count of Deletions (DEL)                                                 |
| `ins_count`              | Integer | Count of Insertions (INS)                                                |
| `variant_density_per_mb` | Numeric | Number of variants per Megabase                                          |

### Exclusion Metrics

**Function:** `load_exclusion_metrics()`

Available for v5.0q benchmarks only.

| Column             | Type      | Description                                          |
| :----------------- | :-------- | :--------------------------------------------------- |
| `bench_version`    | Factor    | Benchmark version (e.g., `v5.0q`)                    |
| `ref`              | Factor    | Reference genome                                     |
| `bench_type`       | Factor    | Benchmark set type                                   |
| `exclusions`       | Character | Name of the exclusion region                         |
| `exclusion_bp`     | Numeric   | Total size of the exclusion region                   |
| `intersect_bp`     | Numeric   | Overlap with benchmark                               |
| `pct_of_exclusion` | Numeric   | Percentage of exclusion region overlapping benchmark |
| `pct_of_dip`       | Numeric   | Percentage of diploid genome overlapping benchmark   |
| `total_variants`   | Integer   | Total variants in this exclusion                     |

### Reference Genome Sizes

**Function:** `load_reference_sizes()`

| Column   | Type    | Description                                      |
| :------- | :------ | :----------------------------------------------- |
| `ref`    | Factor  | Reference genome name                            |
| `chrom`  | Factor  | Chromosome name (standardized with "chr" prefix) |
| `length` | Integer | Total length of the chromosome                   |
| `ns`     | Integer | Number of 'N' bases (gaps)                       |
| `asm_bp` | Integer | Assembled bases (`length` - `ns`)                |

### Variant Table

**Function:** `load_variant_table()` or `load_all_variant_tables()`

Large datasets. Columns may vary slightly between benchmark versions.

| Column          | Type      | Description                                                     |
| :-------------- | :-------- | :-------------------------------------------------------------- |
| `bench_version` | Factor    | Benchmark version                                               |
| `ref`           | Factor    | Reference genome                                                |
| `bench_type`    | Factor    | Benchmark set type                                              |
| `chrom`         | Factor    | Chromosome                                                      |
| `pos`           | Integer   | 1-based start position                                          |
| `end`           | Integer   | End position                                                    |
| `gt`            | Character | Genotype                                                        |
| `vkx`           | Character | Variant class                                                   |
| `var_type`      | Character | Variant type (SNV, INDEL, DEL, INS, COMPLEX, OTHER)             |
| `len_ref`       | Integer   | Length of reference allele                                      |
| `len_alt`       | Integer   | Length of alternate allele                                      |
| `var_size`      | Integer   | Size of the variant                                             |
| `region_ids`    | Character | Comma-separated IDs of regions overlapping the variant          |
| `context_ids`   | Character | Comma-separated IDs of genomic contexts overlapping the variant |

### Genomic Context Coverage

**Function:** `load_genomic_context_coverage()`

Base-level coverage metrics for difficult genomic contexts.

| Column         | Type    | Description                                            |
| :------------- | :------ | :----------------------------------------------------- |
| `context_name` | Factor  | Genomic context identifier                             |
| `chrom`        | Factor  | Chromosome                                             |
| `start`        | Integer | Start position (0-based)                               |
| `end`          | Integer | End position                                           |
| `n_overlap`    | Integer | Number of overlapping intervals                        |
| `bases_cov`    | Integer | Number of bases covered                                |
| `ivl_len`      | Integer | Length of the interval                                 |
| `frac_cov`     | Numeric | Fraction of interval covered (`bases_cov` / `ivl_len`) |

### Benchmark Regions

**Function:** `load_benchmark_regions()`

| Column          | Type    | Description                            |
| :-------------- | :------ | :------------------------------------- |
| `bench_version` | Factor  | Benchmark version                      |
| `ref`           | Factor  | Reference genome                       |
| `bench_type`    | Factor  | Benchmark set type                     |
| `chrom`         | Factor  | Chromosome                             |
| `start`         | Integer | Start position (0-based)               |
| `end`           | Integer | End position                           |
| `interval_size` | Integer | Size of the interval (`end` - `start`) |

### HG002 Q100 Assembly Size

**Function:** `load_hg002q100_size()`

Returns a single numeric value representing the total base pairs of the HG002 Q100 maternal assembly (v1.1).

---

## Pipeline Metric Definitions

### Core Overlap Metrics

#### context_bp

**Definition:** Total number of bases in a genomic context region

**Calculation:** Sum of all bases in the BED file defining that genomic context

**Units:** Bases (bp)

**Range:** Millions to billions depending on genomic context

**Interpretation:**

- Larger values = genomic context covers more of the genome
- Used as denominator for calculating `pct_of_context`
- Example: HP genomic context has ~84 million bases in GRCh38

**Examples by Stratification:**

| Stratification | GRCh37 | GRCh38 | CHM13v2.0 |
| -------------- | ------ | ------ | --------- |
| HP             | 81.2M  | 83.9M  | 85.4M     |
| MAP            | 245.3M | 248.9M | 251.2M    |
| SD             | 159.4M | 166.8M | 170.1M    |
| SD10kb         | 141.2M | 151.8M | 155.3M    |
| TR             | 234.5M | 241.2M | 244.8M    |
| TR10kb         | 198.3M | 205.4M | 208.6M    |

---

#### intersect_bp

**Definition:** Number of bases in a genomic context that overlap with benchmark regions

**Calculation:**

```
bedtools intersect -a genomic context.bed -b benchmark.bed | awk '{sum += $3-$2} END {print sum}'
```

**Units:** Bases (bp)

**Range:** 0 to context_bp (cannot exceed total genomic context size)

**Interpretation:**

- Shows how much of the difficult region is actually included in the benchmark
- Larger values = benchmark better covers this difficult region
- Zero means genomic context not covered by benchmark

---

#### pct_of_context

**Definition:** Percentage of genomic context region included in benchmark

**Formula:**

```
pct_of_context = (intersect_bp / context_bp) × 100
```

**Units:** Percentage (0-100%)

**Range:** 0-100%

**Interpretation:**

- Shows benchmark's completeness for this difficult region type
- **High values (>80%):** Benchmark includes most of this difficult region
- **Medium values (40-80%):** Benchmark includes part of this difficult region
- **Low values (<40%):** Benchmark avoids this difficult region
- Useful for deciding which benchmarks to use for specific analyses

**Example Interpretation:**

```
If pct_of_context = 45% for MAP region:
→ Only 45% of all low-complexity regions are in the benchmark
→ Benchmark may not be suitable for analyzing complexity-related variants
```

---

#### pct_of_bench

**Definition:** Percentage of benchmark regions overlapping with a specific genomic context

**Formula:**

```
pct_of_bench = (intersect_bp / total_benchmark_bp) × 100
```

**Units:** Percentage (0-100%)

**Range:** 0-100%, but typically sums to ~80-120% (regions overlap)

**Interpretation:**

- Shows how much of the benchmark consists of a particular difficult region
- Larger values = this difficult region is a major component of the benchmark
- Values add to ~100% when summing across all genomic contexts (but overlap)
- Useful for understanding benchmark composition

**Example Interpretation:**

```
For HP genomic context: pct_of_bench = 2.92%
→ About 2.92% of the entire benchmark overlaps with homopolymer regions
→ This is the fraction of benchmark bases users will be working with
```

---

#### variant_density_per_mb

**Definition:** Number of variants per megabase in the benchmark-genomic context overlap region

**Formula:**

```
variant_density_per_mb = (total_variants / (intersect_bp / 1,000,000))
```

**Units:** Variants per megabase (v/Mb)

**Range:** Typically 0 - 50,000 (highest in most difficult regions)

**Interpretation:**

- Shows variant frequency in difficult regions
- **Higher values:** Indicates more variant-dense regions
- **Compare across genomic contexts:** Identifies which regions have highest variant density
- **Use for:** Understanding which difficult regions contribute most variants to benchmark

**Example Interpretation:**

```
If variant_density_per_mb = 9310 for HP in small variants:
→ For every million bases in homopolymer regions, there are ~9,310 variants
→ This is high compared to genome average (~10-50 v/Mb)
→ Suggests homopolymers are challenging for variant calling
```

**Density Comparisons:**

- Genome average small variant density: ~10-50 variants/Mb
- Low-complexity regions (MAP): 5,000-8,000 v/Mb
- Homopolymers (HP): 8,000-10,000 v/Mb
- Tandem repeats (TR): 3,000-6,000 v/Mb
- Structural variants: 10-100 v/Mb

---

### Variant Counts

#### total_variants

**Definition:** Total count of all variants in the benchmark-genomic context overlap region

**Units:** Integer count

**Calculation:** Sum of all `*_count` columns (snv_count + indel_count + del_count + ins_count + any additional type counts)

**Variant Classification Logic:**

- **Small variants (≤50bp or missing SVTYPE)**: Classified by REF/ALT length comparison
- **Structural variants (>50bp with SVTYPE)**: Uses INFO/SVTYPE field from VCF

**Interpretation:**

- Represents the total number of variants within this genomic context-benchmark overlap
- Useful for understanding how many variants tools must process in difficult regions
- Larger counts = more challenging region for variant calling

**Context:**

- For small variants: typically thousands to hundreds of thousands per genomic context
- For structural variants: typically tens to thousands per genomic context
- Varies significantly by genomic context type and reference genome

---

#### snv_count

**Definition:** Number of single nucleotide variants (SNVs)

**Note:** Pipeline outputs use `SNP` internally; the R loading function renames to `snv_count`.

**Definition of SNV:**

- Single base substitution
- Reference length = 1 bp, Alternate length = 1 bp
- Example: A → G

**Classification:**

- Classified by Truvari's `VariantRecord.var_type()` as `SV.SNP`
- Applied to all benchmarks (small variants and structural variants)

**Interpretation:**

- SNVs are most common variant type
- Generally easiest to genotype accurately
- Usually comprise 30-50% of total variants in benchmarks

**Context:**

- Small variants: 180,000-700,000 SNVs per benchmark-genomic context
- Structural variants: fewer SNVs (focus is on larger variants ≥50bp)

---

#### indel_count

**Definition:** Number of insertion/deletion variants (INDELs) in small variant benchmarks

**Definition of INDEL:**

- Variant where reference and alternate alleles differ in length
- For small variant benchmarks: variants < 50bp
- Classified by Truvari's `VariantRecord.var_type()` and reclassified as INDEL for smvar benchmarks

**Interpretation:**

- Represents combined insertions and deletions in small variant benchmarks
- Heavily concentrated in homopolymer regions (HP) and tandem repeats (TR)
- Usually 40-65% of total small variants

---

#### del_count

**Definition:** Number of deletion variants (all sizes)

**Definition of Deletion:**

- Removal of bases from reference
- Reference length > Alternate length
- Applies to small deletions (<50bp) and structural deletions (≥50bp)
- Examples: ATG → A (2bp deletion), 1000bp deletion

**Classification:**

- **Small variants (≤50bp)**: Classified when REF_len > ALT_len
- **Structural variants (>50bp with SVTYPE=DEL)**: Uses SVTYPE field from VCF

**Interpretation:**

- Represents loss-of-sequence variants of any size
- Includes both small indels and large structural deletions
- Important for detecting copy number losses
- Distribution varies by genomic context and benchmark type

**Context:**

- Small variant benchmarks: includes deletions 1-49bp
- Structural variant benchmarks: focuses on deletions ≥50bp
- Harder to sequence due to missing coverage for large deletions

---

#### ins_count

**Definition:** Number of insertion variants (all sizes)

**Definition of Insertion:**

- Addition of bases in alternate vs. reference
- Alternate length > Reference length
- Applies to small insertions (<50bp) and structural insertions (≥50bp)
- Examples: A → ATG (2bp insertion), 1000bp insertion

**Classification:**

- **Small variants (≤50bp)**: Classified when ALT_len > REF_len
- **Structural variants (>50bp with SVTYPE=INS)**: Uses SVTYPE field from VCF

**Interpretation:**

- Represents gain-of-sequence variants of any size
- Includes both small indels and large structural insertions
- Important for detecting copy number gains and novel sequences
- Usually 10-50% of structural variants
- May be harder to sequence than deletions

**Context:**

- Structural variants: 0-10,000 insertions per benchmark
- Often more common than deletions in complex regions

---

#### complex_count

**Definition:** Number of complex variants

**Definition of Complex Variant:**

- Both reference and alternate are >1 bp in length
- Cannot be classified as simple SNV or INDEL
- Example: ACGT → TT (deletion + substitution)

**Interpretation:**

- Indicates variants with simultaneous substitution and length change
- Harder to genotype due to multiple differences
- Often results of incorrect alignment or genuine complex events
- Usually <5% of total variants

**Context:**

- Small variants: 0-50,000 complex variants per benchmark
- Structural variants: 0-5,000 complex variants per benchmark
- More common in difficult regions like segmental duplications

---

#### other_count

**Definition:** Number of other variant types, including structural variant types from SVTYPE

**Included Variant Types:**

- **Structural variant types (from SVTYPE field)**: DUP (duplications), INV (inversions), BND (breakends), CNV (copy number variants)
- **Unclassified variants**: Variants that don't fit standard categories

**Classification:**

- For variants >50bp with SVTYPE annotation, includes non-DEL/INS structural variant types
- May also include variants with data quality issues or unusual patterns
- Includes variants of unknown significance

**Interpretation:**

- Usually small proportion of total variants (<5%)
- Investigate if significantly high in a region
- May indicate data processing artifacts

**Context:**

- Small variants: typically 0-50,000 other variants
- Structural variants: typically <1% of calls

---

## Stratification Regions

Stratifications represent different categories of "difficult" genomic regions - areas where variant calling is challenging.

### HP - Homopolymer Regions

**Definition:** Runs of identical nucleotides (≥7 bp)

**Examples:** AAAAAAA, TTTTTTTT, GGGGGGGG, CCCCCCCC

**Biological Significance:**

- Sequencing and alignment errors common in homopolymers
- Indels frequently misaligned due to lack of sequence context
- Particularly problematic for short-read sequencing

**Source & Version:**

- GRCh37/38 definitions from GIAB
- CHM13v2.0 adapted from pangenome references

**Coverage in Benchmarks:**

- ~95-98% of homopolymer regions covered
- Usually represents 2-3% of benchmark bases
- Relatively well-covered across all benchmarks

**Variant Characteristics:**

- Dominated by INDELs (70-80% of variants)
- High density: ~8,000-10,000 variants/Mb
- Most challenging region type for variant calling

**Typical Values (v5.0q_GRCh38_smvar):**

```
context_bp = 83,977,437
intersect_bp = 80,057,238
pct_of_context = 95.3% (well covered)
pct_of_bench = 2.92% (3% of benchmark)
variant_density = 9,310 v/Mb (very high)
```

---

### MAP - Low-Complexity/Mapping Regions

**Definition:** Low-complexity regions with mapping quality issues

**Characteristics:**

- Tandem repeats with short periodicity
- Low sequence complexity (high AT or GC content)
- Regions with multiple optimal alignments
- Often associated with mapping artifacts

**Biological Significance:**

- Reads align poorly due to multiple valid positions
- Consensus mappers produce ambiguous alignments
- Variants in these regions have high false-positive rate

**Coverage in Benchmarks:**

- ~45-50% of mapping regions covered
- Usually represents 4-5% of benchmark bases
- Partial coverage indicates intentional focus on solvable cases

**Variant Characteristics:**

- High SNV proportion (85% SNVs, 15% INDELs)
- Medium density: ~7,000-8,000 variants/Mb
- Many variants result from alignment errors

**Typical Values (v5.0q_GRCh38_smvar):**

```
context_bp = 248,876,839
intersect_bp = 113,370,415
pct_of_context = 45.6% (partially covered)
pct_of_bench = 4.14% (4% of benchmark)
total_variants = 851,471
```

---

### SD - Segmental Duplications

**Definition:** Segmental duplications (>1kb, ≥90% sequence identity)

**Biological Significance:**

>

- > 90% similar sequences elsewhere in genome
- Impossible to uniquely map short reads
- Classic "copy-number variation" regions
- High false-positive and false-negative rate

**Characteristics:**

- Very long duplications (often >100kb)
- Found throughout genome, often pericentromeric
- Rapidly evolving (sequence diverges over evolutionary time)

**Source & Version:**

- UCSC genome browser segmental duplication track
- Consistent across GRCh37/38
- CHM13v2.0 updated with new duplications

**Coverage in Benchmarks:**

- ~48-50% of segmental duplications covered
- Usually represents 3% of benchmark bases
- Most conservative genomic context (heavily avoiding)

**Variant Characteristics:**

- High proportion SNVs (80-85%)
- Medium density: ~6,000-6,500 variants/Mb
- Difficult to validate due to ambiguous mapping

**Typical Values (v5.0q_GRCh38_smvar):**

```
context_bp = 166,860,344
intersect_bp = 81,066,975
pct_of_context = 48.6% (half covered)
pct_of_bench = 2.96% (3% of benchmark)
```

---

### SD10kb - Segmental Duplication Flanks ±10kb

**Definition:** Regions ±10kb flanking segmental duplications

**Purpose:** Captures edge effects near segmental duplications

**Biological Significance:**

- Spillover mapping errors from nearby duplications
- Still challenging but more tractable than exact duplications
- Represents transition zone between difficult and easy regions

**Coverage in Benchmarks:**

- ~44-45% coverage
- Usually represents 2.5% of benchmark bases
- Often included when SD regions are avoided

**Variant Characteristics:**

- Similar profile to SD regions but somewhat less dense
- Density: ~5,500-7,000 variants/Mb
- More amenable to variant calling than exact SD

---

### TR - Tandem Repeats

**Definition:** Tandem repeats of any length and periodicity

**Examples:**

- Satellite DNA (large repeats, >100bp period)
- Microsatellites (short repeats, 1-6bp period)
- Minisatellites (medium repeats, 10-100bp period)

**Biological Significance:**

- Highly polymorphic in human populations
- Difficult to align and phase
- Common cause of failed variant calling
- Important for population genetics and identification

**Source & Version:**

- RepeatMasker annotations
- Consistent across GRCh37/38
- CHM13v2.0 includes additional repeats

**Coverage in Benchmarks:**

- ~40-50% of tandem repeats covered
- Usually represents 4-7% of benchmark bases
- Moderate coverage reflecting difficulty

**Variant Characteristics:**

- Mixed SNVs and INDELs
- Medium density: ~5,000-7,000 variants/Mb
- Many variants are repeat-size variations

**Typical Values (v5.0q_GRCh38_smvar):**

```
context_bp = 241,276,584
intersect_bp = varies by region
pct_of_bench = 4-7% (large contributor to benchmark)
```

---

### TR10kb - Tandem Repeat Flanks ±10kb

**Definition:** Regions ±10kb flanking tandem repeats

**Purpose:** Captures sequence context around tandem repeats

**Biological Significance:**

- Edge effects from repeat polymorphism
- Can affect calling in flanking sequences
- More tractable than exact repeat regions

**Coverage in Benchmarks:**

- Similar to TR coverage (40-50%)
- Usually represents 5-6% of benchmark bases
- Often included alongside TR regions

---

## Exclusion Categories (v5.0q only)

Exclusions represent regions explicitly removed from v5.0q benchmarks.

### consecutive-svs

**Definition:** Regions with multiple consecutive structural variants in source data

**Rationale:**

- Complex regions with uncertain true variant structure
- High risk of incorrect SV breakpoint determination
- Better to exclude than include unreliable calls

**Typical Impact:**

- Usually 10-20% of SV benchmark bases excluded
- Creates gaps in benchmark coverage
- Concentrated in pericentromeric and telomeric regions

---

### flanks

**Definition:** Flanking regions around large features (excluded)

**Rationale:**

- Capture-edge effects near major variants
- Sequencing depth variations at boundaries
- Alignment complexity near breakpoints

**Typical Impact:**

- Smaller exclusion, ~5-10% of benchmark
- Affects edges of large variants
- Reduces false positives at variant boundaries

---

### satellites

**Definition:** Satellite DNA and alpha-satellite regions

**Examples:**

- Alpha satellite DNA (pericentromeric)
- Beta satellite DNA
- Satellite repeats (classical satellites)

**Rationale:**

- Highly repetitive, unmappable with short reads
- No unique mapping possible
- Variants fundamentally uncallable

**Typical Impact:**

- 5-15% of benchmark bases
- Concentrated at centromeres and heterochromatin
- Most conservative exclusion (rightfully excluded)

---

### segdups

**Definition:** Segmental duplications (>90% identity, >1kb)

**Rationale:**

- Same as SD genomic context
- Duplicate regions unmappable with short reads
- Excluded from benchmark as uncallable

**Typical Impact:**

- 20-30% of benchmark bases
- Largest single exclusion category
- But still 50% of segdups are included in benchmark

---

### LCR-unique

**Definition:** Low-complexity repetitive regions (unique mapping)

**Rationale:**

- Some LCR regions have unique mapping sites
- May be callable despite low complexity
- Excluded only if truly problematic

**Typical Impact:**

- 10-20% of benchmark bases
- Overlaps with MAP genomic context
- Partially included to retain challenging-but-solvable regions

---

### tandem-repeats

**Definition:** All tandem repeats (including satellites)

**Rationale:**

- Tandem repeats difficult to align and phase
- Variants ambiguous due to repeat structure
- Excluded from most benchmarks historically

**Typical Impact:**

- 15-25% of benchmark bases
- But 40-50% still included (TR genomic context)
- Reflects modern ability to call some repeat variants

---

## Variant Classifications

### SNV (Single Nucleotide Variant)

**Definition:** Single base substitution

**Characteristics:**

- Reference allele: 1 base (A, C, G, or T)
- Alternate allele: 1 different base
- No length change

**Examples:**

- A → G (transition)
- C → T (transition)
- A → C (transversion)
- G → T (transversion)

**Prevalence:** 30-50% of variants in small variant benchmarks

**Callability:** Easiest variant type to call accurately

**Stratification Pattern:**

- Higher density in difficult regions
- Particularly common in low-complexity regions (MAP)
- Genome average: 3-5 SNVs per 10,000 bases

---

### INDEL (Insertion or Deletion)

**Definition:** Insertion or deletion of bases (size 1-49 bp for small variants)

**Characteristics:**

- Net change in sequence length
- Reference length ≠ Alternate length
- Typically 1-49 bp for small variant benchmarks
- 50+ bp classified as structural variants

**Examples:**

- AT → A (deletion of T)
- ACGT → ACGATGT (insertion of AT)
- TTT → T (deletion of TT)

**Size Distribution:**

- Median size: 2-5 bp
- Range: 1-49 bp for small variants
- Longer indels ~50+ bp classified as structural variants

**Prevalence:** 40-65% of variants in small variant benchmarks

**Callability:** Harder than SNVs due to alignment complexity, easier than complex variants

**Stratification Pattern:**

- Heavily concentrated in homopolymer regions (HP)
- Also common in tandem repeats (TR)
- Rarer in low-complexity regions (MAP)

**Special Cases:**

- **Homopolymer indels:** Most common cause of alignment errors
- **Repeat indels:** Often represent repeat-size variations
- **Complex indels:** May be misaligned single variants

---

### DEL (Deletion, Structural Variants)

**Definition:** Deletion of ≥50 bp (structural variant classification)

**Characteristics:**

- Reference longer than alternate
- Size: ≥50 bp
- Often megabase-scale

**Examples:**

- Deletion of 1 kb: 1000 bp sequence → missing in alternate
- Deletion of 1 Mb: chromosome segment missing

**Prevalence:** 10-50% of structural variants

**Callability:** Depends on size and flanking sequence; generally harder than small variants

**Detection Methods:**

- Read depth analysis (coverage drop)
- Split reads (paired reads spanning breakpoints)
- Assembly-based methods

**Biological Significance:**

- Copy number losses
- Gene disruptions
- Regulatory element deletions

---

### INS (Insertion, Structural Variants)

**Definition:** Insertion of ≥50 bp (structural variant classification)

**Characteristics:**

- Alternate longer than reference
- Size: ≥50 bp
- Often megabase-scale

**Examples:**

- Insertion of 1 kb sequence
- Insertion of transposable element
- Duplication of chromosome segment

**Prevalence:** 10-50% of structural variants

**Callability:** Often harder than deletions because sequence must be present in reads

**Detection Methods:**

- Split reads (reads spanning insertion)
- Pair-end mapping (insert size anomalies)
- Assembly-based methods

**Biological Significance:**

- Copy number gains
- Transposable element insertions
- Regulatory element duplications

---

### COMPLEX (Complex Variants)

**Definition:** Variants with simultaneous substitution and length change

**Characteristics:**

- Both reference and alternate >1 bp
- Cannot be decomposed to simple SNV + INDEL

**Examples:**

- ACGT → TT (deletion + substitution)
- A → CGT (insertion + substitution)
- GATC → TC (deletion + substitution)

**Prevalence:** <5% of variants typically

**Callability:** Most difficult variant type due to ambiguity

**Causes:**

- True complex events (rare)
- Alignment artifacts (common)
- Overlapping variants incorrectly merged
- VCF representation choices

**Handling Strategies:**

- Some tools decompose to SNV + INDEL
- Others treat as single variant
- Benchmark representation varies

---

### OTHER (Unclassified)

**Definition:** Variants not fitting standard categories

**When Occurs:**

- Variants with no clear classification
- Data processing artifacts
- Potentially poor-quality calls

**Prevalence:** <5% of variants

**Interpretation:**

- Should be minimal
- Investigate if percentage is high
- May indicate upstream quality issues

---

## Chromosome Naming Conventions

### Standard Naming

- **Autosomes:** chr1 - chr22 (GRCh convention with "chr" prefix)
- **Sex chromosomes:** chrX, chrY (not chr23/24)
- **Mitochondrial:** chrM (not chrMT or MT)

### Sources of Variation

- Some resources use "1" instead of "chr1"
- Older data may use haplotype-specific names (e.g., "haplotype1.1")
- Mitochondrial sometimes written as "MT" or "chrMT"

### Our Convention

- All files standardized to use "chr" prefix
- All chromosomes in reference size files have "chr" prefix
- Variant files inherit from source but usually include "chr"

---

## Summary Statistics by Benchmark

### v5.0q_GRCh38 Example

**Small Variants (smvar):**

- Total variants: ~4.9M
- Average density: 6,500 variants/Mb
- Benchmark coverage: 3.0 Gb (94% of assembled genome)
- Composition: ~35% SNVs, ~60% INDELs, ~5% complex

**Structural Variants (stvar):**

- Total variants: ~14,000
- Density: ~4.7 variants/Mb
- Composition: ~50% deletions, ~45% insertions, ~5% complex

### v5.0q_GRCh37 Example

**Small Variants:**

- Total variants: ~4.6M
- Similar distribution to GRCh38
- Benchmark coverage: 2.9 Gb

**Structural Variants:**

- Total variants: ~13,000
- Similar density and composition to GRCh38

---

## Related Documentation

- **[Pipeline Outputs Reference](pipeline-outputs.md)** - File formats and loading instructions
- **[Output Relationships Diagram](diagrams/output-relationships.mmd)** - Visual guide to data relationships
- **[Architecture Overview](architecture.md)** - How metrics are calculated in the pipeline
