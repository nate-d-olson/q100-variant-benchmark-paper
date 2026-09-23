# Figure 1 legend — concept-d rev 3.1

Drafted 2026-08-18 against `HG002v5-variant-benchmark-manuscript_20260814_WERB.docx`.
Replaces the previous caption ("DeFrABB pipeline workflow overview with specific
configuration used to create the v5.0q benchmark set…"). All terms verified against
the manuscript: component names (Results/Methods §DeFrABB), Miqa curation details
(Methods §External Evaluation Manual Curation), platforms (Table 2), exclusions
(Table 3), curation classes (Figure 6), comparison-callset error framing (Figure 7).

---

**Figure 1: Development of the GIAB HG002 v5.0q whole-genome variant benchmarks.**
Top: the GIAB benchmark development cycle. The HG002 Q100 diploid assembly and
three reference genomes (GRCh37, GRCh38, CHM13) are inputs to the automated
DeFrABB pipeline (panel a), which produces draft small-variant and structural
variant (SV) benchmark sets. Drafts are publicly released for community use,
formally evaluated by external groups (panel b), and assessed against the RIDE
principle (panel c). Red dashed paths mark checkpoints where community feedback,
benchmark errors, or a failed assessment lead to revised parameters and
exclusions and a DeFrABB re-run; benchmarks that pass become official GIAB
releases. **(a)** The three DeFrABB components. (1) Assembly-based variant
calling: each haplotype is aligned to the reference; ribbons connect matched
regions and narrow where large insertions (INS) or deletions (DEL) change region
size, yielding variant calls (VCF) and diploid regions with 1:1 alignment of
both haplotypes (BED). (2) Draft benchmark generation: benchmark regions are
diploid regions minus exclusions (assembly gaps, alignment breaks in large
repeats, complex SVs, and known assembly and variant-calling errors); phased
variants of all sizes and benchmark regions are produced for small variants and
SVs on all three references. (3) Evaluation: automated comparisons to existing
callsets generate an analysis report used as internal QC for parameter
optimization and rapid iteration. **(b)** External evaluation. Draft benchmarks
(GRCh38) are compared with callsets contributed by external groups from multiple
sequencing technologies using hap.py/vcfeval (small variants) and truvari (SVs).
A stratified random subset of discrepancies is manually curated by each
submitting group in the Miqa platform (Magna Labs) using pre-configured IGV
sessions showing multi-technology read evidence on GRCh38 and both Q100
haplotypes; each discrepancy is classified as correct, incorrect, or unsure in
the benchmark. Incorrect and unsure calls are re-curated at NIST, and confirmed
benchmark errors feed back into the development cycle. **(c)** The RIDE
("reliable identification of errors") principle (Olson, Nat. Rev. Genet. 24,
2023). A benchmark is fit for purpose when discrepancies with high-quality
comparison callsets are reliably errors in the callsets rather than in the
benchmark. Within benchmark regions (BED), calls unique to the benchmark are
false negatives (FNs) and calls unique to a comparison callset are false
positives (FPs); curated discrepancy outcomes determine whether the benchmark
passes to official release. Schematic; genomic features are not to scale.
