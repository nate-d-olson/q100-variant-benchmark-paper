"""
VCF processing rules for variant characterization
"""


rule generate_sv_len:
    input:
        vcf=config["vcf_files"]["chm13_stvar"],
        bed=config["bed_files"]["chm13_stvar"],
    output:
        tsv=config["outputs"]["sv_len"],
    log:
        "logs/vcf_processing/sv_len.log",
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools query -HH \
            -f '%CHROM\\t%POS0\\t%END\\t%SVTYPE\\t%VKX\\t%SVLEN\\t%TRF\\t%TRFstart\\t%TRFend\\t%RM_clsfam\\n' \
            -R {input.bed} {input.vcf} \
            > {output.tsv} 2> {log}
        
        if [ ! -s {output.tsv} ]; then
            echo "ERROR: Output file is empty" >> {log}
            exit 1
        fi
        """
