"""
Chr8 synteny figure rules.

ref/mat/pat must be keys from config.references; the existing prepare_reference
rule handles downloading, bgzip conversion, and indexing automatically.

Pipeline steps:
  1. (prepare_reference)  — download + index ref/mat/pat FASTAs [downloads.smk]
  2. chr8_extract_contig  — extract the chr8 contig from each FASTA
  3. chr8_index_fasta     — samtools index + write chromosome-length (.cl) file
  4. chr8_align           — align mat/pat chr8 to reference chr8 (minimap2)
  5. chr8_syri            — structural rearrangement identification (SyRI)
  6. chr8_find_inversion  — detect inversion coordinates for the zoom panel
  7. chr8_make_figure     — assemble multi-panel synteny figure (plotsr)
"""

_CHR8_CFG = config.get("chr8_synteny", {})
_CHR8_CHROM = _CHR8_CFG.get("chrom", "chr8")
_CHR8_CONTIG = {
    "ref": _CHR8_CFG.get("ref_contig", _CHR8_CHROM),
    "mat": _CHR8_CFG.get("mat_contig", f"{_CHR8_CHROM}_MATERNAL"),
    "pat": _CHR8_CFG.get("pat_contig", f"{_CHR8_CHROM}_PATERNAL"),
}
_CHR8_THREADS = int(_CHR8_CFG.get("threads", 8))
_CHR8_ZOOM_PADDING = int(_CHR8_CFG.get("zoom_padding", 2_000_000))
_CHR8_MIN_INV_SIZE = int(_CHR8_CFG.get("min_inversion_size", 0))


rule chr8_extract_contig:
    """Extract chr8 contig from each configured full-assembly FASTA."""
    input:
        fasta=lambda w: f"resources/references/{_CHR8_CFG[w.sample]}.fa.gz",
        fai=lambda w: f"resources/references/{_CHR8_CFG[w.sample]}.fa.gz.fai",
    output:
        fa="results/chr8_synteny/fasta/{sample}_chr8.fa",
    params:
        contig=lambda wildcards: _CHR8_CONTIG[wildcards.sample],
        target_chrom=_CHR8_CHROM,
    log:
        "logs/chr8_synteny/extract_{sample}.log",
    resources:
        mem_mb=4096,
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        sample="ref|mat|pat",
    shell:
        """
        set -euo pipefail
        echo "Extracting {params.contig} → {params.target_chrom} ({wildcards.sample})" > {log}
        echo "Input: {input.fasta}" >> {log}
        echo "Started at $(date)" >> {log}

        samtools faidx {input.fasta} {params.contig} 2>> {log} \
            | sed "1s/^>.*/>{params.target_chrom}/" > {output.fa}

        echo "Completed at $(date)" >> {log}
        """


rule chr8_index_fasta:
    """Index extracted chr8 FASTA and write chromosome-length (.cl) file for plotsr."""
    input:
        fa="results/chr8_synteny/fasta/{sample}_chr8.fa",
    output:
        fai="results/chr8_synteny/fasta/{sample}_chr8.fa.fai",
        cl="results/chr8_synteny/fasta/{sample}_chr8.cl",
    log:
        "logs/chr8_synteny/index_{sample}.log",
    resources:
        mem_mb=2048,
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        sample="ref|mat|pat",
    shell:
        """
        echo "Indexing {input.fa}" > {log}
        samtools faidx {input.fa} 2>> {log}
        cut -f1,2 {output.fai} > {output.cl} 2>> {log}
        echo "Completed at $(date)" >> {log}
        """


rule chr8_align:
    """Align query chr8 to reference chr8 with minimap2 (asm5 preset).

    Supports any ref/query pair: ref_mat, ref_pat, mat_pat.
    mat_pat is needed so plotsr receives consecutive-genome SyRI files
    (REF↔MAT, MAT↔PAT) rather than two REF-anchored files.
    """
    input:
        ref="results/chr8_synteny/fasta/{ref_samp}_chr8.fa",
        qry="results/chr8_synteny/fasta/{qry_samp}_chr8.fa",
    output:
        bam="results/chr8_synteny/alignments/{ref_samp}_{qry_samp}.bam",
        bai="results/chr8_synteny/alignments/{ref_samp}_{qry_samp}.bam.bai",
    log:
        "logs/chr8_synteny/align_{ref_samp}_{qry_samp}.log",
    threads: _CHR8_THREADS
    resources:
        mem_mb=16384,
    conda:
        "../envs/plotsr.yaml"
    wildcard_constraints:
        ref_samp="ref|mat",
        qry_samp="mat|pat",
    shell:
        """
        set -euo pipefail
        echo "Aligning {wildcards.qry_samp} chr8 to {wildcards.ref_samp} chr8" > {log}
        echo "Started at $(date)" >> {log}

        minimap2 -ax asm5 --eqx -t {threads} {input.ref} {input.qry} 2>> {log} \
            | samtools sort -@ {threads} -O BAM -o {output.bam} - 2>> {log}
        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}

        echo "Completed at $(date)" >> {log}
        """


rule chr8_syri:
    """Run SyRI structural rearrangement identification.

    Supports any ref/query pair matching chr8_align: ref_mat, ref_pat, mat_pat.
    """
    input:
        bam="results/chr8_synteny/alignments/{ref_samp}_{qry_samp}.bam",
        ref="results/chr8_synteny/fasta/{ref_samp}_chr8.fa",
        qry="results/chr8_synteny/fasta/{qry_samp}_chr8.fa",
    output:
        syri="results/chr8_synteny/syri/{ref_samp}_{qry_samp}syri.out",
        summary="results/chr8_synteny/syri/{ref_samp}_{qry_samp}syri.summary",
    params:
        outdir="results/chr8_synteny/syri",
        prefix="{ref_samp}_{qry_samp}",
    log:
        "logs/chr8_synteny/syri_{ref_samp}_{qry_samp}.log",
    threads: _CHR8_THREADS
    resources:
        mem_mb=8192,
    conda:
        "../envs/plotsr.yaml"
    wildcard_constraints:
        ref_samp="ref|mat",
        qry_samp="mat|pat",
    shell:
        """
        set -euo pipefail
        echo "Running SyRI: {wildcards.ref_samp} vs {wildcards.qry_samp}" > {log}
        echo "Started at $(date)" >> {log}

        syri -c {input.bam} -r {input.ref} -q {input.qry} -F B \
            --dir {params.outdir} --prefix {params.prefix} --nc {threads} 2>> {log} || {{
            echo "ERROR: syri failed with exit code $?" >> {log}
            exit 1
        }}

        echo "Completed at $(date)" >> {log}
        """


rule chr8_find_inversion:
    """Detect largest inversion coordinates from SyRI output for the zoom panel.

    Parses the ref-vs-pat SyRI output to find the largest inversion on the
    target chromosome and writes reference/query coordinates plus zoom region
    boundaries to a JSON file.  This file drives the zoom region in the figure.
    """
    input:
        syri_rp="results/chr8_synteny/syri/ref_patsyri.out",
    output:
        coords="results/chr8_synteny/inversion_coords.json",
    params:
        chrom=_CHR8_CHROM,
        min_inv_size=_CHR8_MIN_INV_SIZE,
        zoom_padding=_CHR8_ZOOM_PADDING,
    log:
        "logs/chr8_synteny/find_inversion.log",
    resources:
        mem_mb=2048,
    conda:
        "../envs/plotsr.yaml"
    shell:
        """
        echo "Detecting largest inversion on {params.chrom}" > {log}
        echo "Min inversion size: {params.min_inv_size} bp" >> {log}
        echo "Zoom padding: {params.zoom_padding} bp" >> {log}
        echo "Started at $(date)" >> {log}

        python workflow/scripts/find_chr8_inversion.py \
            --syri-rp {input.syri_rp} \
            --chrom {params.chrom} \
            --min-inv-size {params.min_inv_size} \
            --zoom-padding {params.zoom_padding} \
            --output {output.coords} >> {log} 2>&1

        echo "Inversion coords: $(cat {output.coords})" >> {log}
        echo "Completed at $(date)" >> {log}
        """


rule chr8_make_figure:
    """Generate the chr8 multi-panel synteny figure (full view + zoomed inversion)."""
    input:
        ref_cl="results/chr8_synteny/fasta/ref_chr8.cl",
        mat_cl="results/chr8_synteny/fasta/mat_chr8.cl",
        pat_cl="results/chr8_synteny/fasta/pat_chr8.cl",
        syri_rm="results/chr8_synteny/syri/ref_matsyri.out",
        syri_mp="results/chr8_synteny/syri/mat_patsyri.out",
        coords="results/chr8_synteny/inversion_coords.json",
    output:
        pdf="results/chr8_synteny/chr8_figure.pdf",
        png="results/chr8_synteny/chr8_figure.png",
    params:
        chrom=_CHR8_CHROM,
        out_base=lambda w, output: os.path.splitext(output.pdf)[0],
    log:
        "logs/chr8_synteny/figure.log",
    resources:
        mem_mb=8192,
    conda:
        "../envs/plotsr.yaml"
    shell:
        """
        echo "Generating chr8 synteny figure" > {log}
        echo "Started at $(date)" >> {log}

        python workflow/scripts/make_chr8_figure.py \
            --ref {input.ref_cl} \
            --mat {input.mat_cl} \
            --pat {input.pat_cl} \
            --rm {input.syri_rm} \
            --mp {input.syri_mp} \
            --coords {input.coords} \
            --chrom {params.chrom} \
            --out {params.out_base} >> {log} 2>&1

        echo "Completed at $(date)" >> {log}
        """


rule chr8_synteny:
    """Convenience target: generate all chr8 synteny outputs."""
    input:
        "results/chr8_synteny/chr8_figure.pdf",
        "results/chr8_synteny/chr8_figure.png",
