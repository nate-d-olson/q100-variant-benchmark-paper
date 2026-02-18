"""
chr8_synteny.smk — Rules for the Chr8 multi-panel synteny figure

Produces:
  results/chr8_synteny/chr8_figure.pdf   <- publication-quality 3-panel figure
  results/chr8_synteny/chr8_figure.png

Pipeline
--------
1. Index FASTAs (samtools faidx) and extract .cl chromosome-length files
2. Pairwise alignment: minimap2 (REF vs MAT, MAT vs PAT)
3. Structural annotation: syri
4. make_chr8_figure.py: run plotsr and assemble 3-panel figure

Configuration
-------------
Set file paths and inversion coordinates in config/config.yaml under
the `chr8_synteny` key. Example:

  chr8_synteny:
    ref: "assemblies/hg002_ref_chr8.fa"
    mat: "assemblies/hg002_mat_chr8.fa"
    pat: "assemblies/hg002_pat_chr8.fa"
    chrom: "chr8"
    inv_start: 50000000
    inv_end: 58000000
    zoom_padding: 2000000
    threads: 16
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

_CHR8_CFG = config.get("chr8_synteny", {})

_CHR8_REF = _CHR8_CFG.get("ref", "")
_CHR8_MAT = _CHR8_CFG.get("mat", "")
_CHR8_PAT = _CHR8_CFG.get("pat", "")
_CHR8_CHROM = _CHR8_CFG.get("chrom", "chr8")
_CHR8_INV_START = _CHR8_CFG.get("inv_start", 50_000_000)
_CHR8_INV_END = _CHR8_CFG.get("inv_end", 58_000_000)
_CHR8_ZOOM_PAD = _CHR8_CFG.get("zoom_padding", 2_000_000)
_CHR8_THREADS = _CHR8_CFG.get("threads", 8)


# ---------------------------------------------------------------------------
# Chromosome-length files (.cl) — faster than re-parsing FASTA in plotsr
# ---------------------------------------------------------------------------


rule chr8_faidx:
    """Index a FASTA assembly."""
    input:
        "{sample}.fa",
    output:
        "{sample}.fa.fai",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"


rule chr8_make_cl:
    """Convert .fai to 2-col chromosome-length file for plotsr (ft:cl)."""
    input:
        "{sample}.fa.fai",
    output:
        "{sample}.cl",
    shell:
        "cut -f1,2 {input} > {output}"


# ---------------------------------------------------------------------------
# Minimap2 pairwise alignments
# ---------------------------------------------------------------------------


rule chr8_align_ref_mat:
    """Align REF vs MAT chromosome assemblies."""
    input:
        ref=_CHR8_REF,
        qry=_CHR8_MAT,
    output:
        bam="results/chr8_synteny/alignments/ref_mat.bam",
        bai="results/chr8_synteny/alignments/ref_mat.bam.bai",
    threads: _CHR8_THREADS
    conda:
        "../envs/plotsr.yaml"
    shell:
        """
        minimap2 -ax asm5 -t {threads} --eqx {input.ref} {input.qry} \
          | samtools sort -O BAM - > {output.bam}
        samtools index {output.bam}
        """


rule chr8_align_mat_pat:
    """Align MAT vs PAT chromosome assemblies."""
    input:
        ref=_CHR8_MAT,
        qry=_CHR8_PAT,
    output:
        bam="results/chr8_synteny/alignments/mat_pat.bam",
        bai="results/chr8_synteny/alignments/mat_pat.bam.bai",
    threads: _CHR8_THREADS
    conda:
        "../envs/plotsr.yaml"
    shell:
        """
        minimap2 -ax asm5 -t {threads} --eqx {input.ref} {input.qry} \
          | samtools sort -O BAM - > {output.bam}
        samtools index {output.bam}
        """


# ---------------------------------------------------------------------------
# SyRI structural variant annotation
# ---------------------------------------------------------------------------


rule chr8_syri_ref_mat:
    """Run SyRI on REF vs MAT alignment."""
    input:
        bam=rules.chr8_align_ref_mat.output.bam,
        ref=_CHR8_REF,
        qry=_CHR8_MAT,
    output:
        syri="results/chr8_synteny/syri/ref_matsyri.out",
        summary="results/chr8_synteny/syri/ref_matsyri.summary",
    params:
        prefix="results/chr8_synteny/syri/ref_mat",
    conda:
        "../envs/plotsr.yaml"
    shell:
        """
        syri -c {input.bam} -r {input.ref} -q {input.qry} \
             -F B --prefix {params.prefix} --nc 8
        """


rule chr8_syri_mat_pat:
    """Run SyRI on MAT vs PAT alignment."""
    input:
        bam=rules.chr8_align_mat_pat.output.bam,
        ref=_CHR8_MAT,
        qry=_CHR8_PAT,
    output:
        syri="results/chr8_synteny/syri/mat_patsyri.out",
        summary="results/chr8_synteny/syri/mat_patsyri.summary",
    params:
        prefix="results/chr8_synteny/syri/mat_pat",
    conda:
        "../envs/plotsr.yaml"
    shell:
        """
        syri -c {input.bam} -r {input.ref} -q {input.qry} \
             -F B --prefix {params.prefix} --nc 8
        """


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------


rule chr8_make_figure:
    """Generate 3-panel Chr8 synteny figure (Panel A: full, B: zoom, C: legend)."""
    input:
        ref_cl=_CHR8_REF.replace(".fa", ".cl") if _CHR8_REF.endswith(".fa") else _CHR8_REF,
        mat_cl=_CHR8_MAT.replace(".fa", ".cl") if _CHR8_MAT.endswith(".fa") else _CHR8_MAT,
        pat_cl=_CHR8_PAT.replace(".fa", ".cl") if _CHR8_PAT.endswith(".fa") else _CHR8_PAT,
        syri_rm=rules.chr8_syri_ref_mat.output.syri,
        syri_mp=rules.chr8_syri_mat_pat.output.syri,
        script=workflow.source_path("../scripts/make_chr8_figure.py"),
    output:
        pdf="results/chr8_synteny/chr8_figure.pdf",
        png="results/chr8_synteny/chr8_figure.png",
    params:
        inv_start=_CHR8_INV_START,
        inv_end=_CHR8_INV_END,
        zoom_padding=_CHR8_ZOOM_PAD,
        chrom=_CHR8_CHROM,
        out_base="results/chr8_synteny/chr8_figure",
    conda:
        "../envs/plotsr.yaml"
    shell:
        """
        python {input.script} \
            --ref   {input.ref_cl} \
            --mat   {input.mat_cl} \
            --pat   {input.pat_cl} \
            --rm    {input.syri_rm} \
            --mp    {input.syri_mp} \
            --chrom {params.chrom} \
            --inv-start    {params.inv_start} \
            --inv-end      {params.inv_end} \
            --zoom-padding {params.zoom_padding} \
            --out   {params.out_base}
        """
