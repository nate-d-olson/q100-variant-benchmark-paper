#!/usr/bin/env python3
"""
Merge exclusion configuration into config.yaml.
"""

import json
import sys


def merge_configs(output_path):
    """Read JSON exclusion config and write merged YAML to file."""
    # Read the JSON exclusion config
    with open("/tmp/exclusion_config.json") as f:
        exclusion_config = json.load(f)
    
    lines = []
    lines.append("# Q100 Variant Benchmark Analysis Pipeline Configuration")
    lines.append("")
    lines.append("# Benchmark sets for RTG vcfstats analysis")
    lines.append("benchmarksets:")
    
    # Add exclusion entries for each benchmark set that has them
    benchmark_order = [
        "v5q_chm13_smvar",
        "v5q_chm13_stvar", 
        "v5q_grch37_smvar",
        "v5q_grch37_stvar",
        "v5q_grch38_smvar",
        "v5q_grch38_stvar",
    ]
    
    # VCF/BED URLs for each benchmark
    benchmark_urls = {
        "v5q_chm13_smvar": {
            "vcf_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz",
            "vcf_checksum": "401375f823c4add943ebac3602c8c753951d98363b8624b831ee2682cc791fad",
            "bed_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            "bed_checksum": "224a6b65769765700d6333a432ba21ff303e54cd9475b79fbe1efea7e1786417",
        },
        "v5q_chm13_stvar": {
            "vcf_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz",
            "vcf_checksum": "011a604daa45c2a5c13b3403436ab16185be64405870d9636a84d55e19736db7",
            "bed_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed",
            "bed_checksum": "6b9cc176bbfda0375a3979d3f00524514f7dfbff96b7948425634ae9a0eb81ff",
        },
        "v5q_grch37_smvar": {
            "vcf_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh37_HG2-T2TQ100-V1.1_smvar.vcf.gz",
            "vcf_checksum": "972bc22ebf7ad23c645ac1c7535f0cfb5c1895f7c5378a12c2da3fddf5e6d7e6",
            "bed_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh37_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            "bed_checksum": "91d0d56fb9a0fb954c6b97a92f914fe6eb25d5718c8de9d4907e0f20e371722a",
        },
        "v5q_grch37_stvar": {
            "vcf_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh37_HG2-T2TQ100-V1.1_stvar.vcf.gz",
            "vcf_checksum": "edb582ceec508f6d745acd0a8c522ee4f235ff5bf4754ba6bdc4c3abee122b5c",
            "bed_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh37_HG2-T2TQ100-V1.1_stvar.benchmark.bed",
            "bed_checksum": "fbdec183e5b0ba83efdb5880ae038db98cd12cd9715e28debb9ed96541c10f40",
        },
        "v5q_grch38_smvar": {
            "vcf_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz",
            "vcf_checksum": "c7f9d7a4781f1f0a8a42ec54b051d6b95db1935c5e12318bae219740f9c50dc8",
            "bed_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            "bed_checksum": "ef7d898dd421ea357bb9db915894ab3e5c6c1689e36ccf619befda4bbf1cef5b",
        },
        "v5q_grch38_stvar": {
            "vcf_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz",
            "vcf_checksum": "d66b2d2496ff5418763813d3195599007dbff950514d10e5bc27b6c8b76e34b8",
            "bed_url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed",
            "bed_checksum": "2f75ce942e1dd9e1a443e4e04fac640104aae5e23fba24d7d3f2c1f1f9985b00",
        },
    }
    
    for benchmark_id in benchmark_order:
        if benchmark_id not in exclusion_config:
            continue
        
        data = exclusion_config[benchmark_id]
        urls = benchmark_urls[benchmark_id]
        
        lines.append(f"  {benchmark_id}:")
        lines.append(f"    vcf:")
        lines.append(f"      url: {urls['vcf_url']}")
        lines.append(f"      checksum: {urls['vcf_checksum']}")
        lines.append(f"    bed:")
        lines.append(f"      url: {urls['bed_url']}")
        lines.append(f"      checksum: {urls['bed_checksum']}")
        lines.append(f"    dip_bed:")
        lines.append(f"      path: \"{data['dip_bed']['path']}\"")
        lines.append(f"      sha256: \"{data['dip_bed']['sha256']}\"")
        lines.append(f"    exclusions:")
        
        for excl in data['exclusions']:
            lines.append(f"      - name: \"{excl['name']}\"")
            lines.append(f"        type: \"{excl['type']}\"")
            lines.append(f"        files:")
            for file_info in excl['files']:
                lines.append(f"          - path: \"{file_info['path']}\"")
                lines.append(f"            sha256: \"{file_info['sha256']}\"")
    
    # Add historical benchmarks (without exclusions)
    lines.append("  v421_grch38_smvar:")
    lines.append("    vcf:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
    lines.append("      checksum: adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c")
    lines.append("    bed:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
    lines.append("      checksum: fba9a57c36ec88e5d14ea3e259c8866c7935f4998d3ec0fa2d6c3962da5b5575")
    
    lines.append("  cmrg_grch38_smvar:")
    lines.append("    vcf:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz")
    lines.append("      checksum: 80f18efa283d953cfc615d30dd01b1ff80debcc46f52042b942da783acd48da2")
    lines.append("    bed:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed")
    lines.append("      checksum: dc77803bb474802fd3213d2f3649076dc14675d3b3c0ae3284036be09d2813f7")
    
    lines.append("  cmrg_grch38_stvar:")
    lines.append("    vcf:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz")
    lines.append("      checksum: 06a8740ad93dd0ede325fcf04c85f687a835699da715c72bdcfade5d588fd19e")
    lines.append("    bed:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed")
    lines.append("      checksum: 3306319595ffa472d8ba9a8ca8932304b01ec663e102900a61604410737ddf81")
    
    lines.append("  v06_grch37_stvar:")
    lines.append("    vcf:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz")
    lines.append("      checksum: 1bc8e515189e32eb02a8a2050a4f70e88a80d9ad1d9fafe65ca1931d6f7b7f13")
    lines.append("    bed:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed")
    lines.append("      checksum: 962da1439091836346e3c76b2cb4c1191586ed2f8ce49bc3ea49eeec75b4ab2b")
    
    lines.append("  tr_grch38_trvar:")
    lines.append("    vcf:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/TandemRepeats_v1.0/GRCh38/HG002_GRCh38_TandemRepeats_v1.0.1.vcf.gz")
    lines.append("      checksum: ea8e1805c1b615ce7f121da5ea1695a3566f26658763ab17a40338d93c8cbf86")
    lines.append("    bed:")
    lines.append("      url: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/TandemRepeats_v1.0/GRCh38/HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz")
    lines.append("      checksum: bd4be0fca501b023ae87fefc1708303a9dd899ce9c0bcb972395b92987a5c0fc")
    
    # References section
    lines.append("")
    lines.append("references:")
    lines.append("  hg002v1.1_fai:")
    lines.append("    path: \"data/hg002v1.1.mat_Y_EBV_MT.fasta.gz.fai\"")
    lines.append("    url: \"https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.mat_Y_EBV_MT.fasta.gz.fai\"")
    lines.append("    checksum: \"c84f852b1cd4a00b12e6439ae7a2dd87\"")
    lines.append("  grch37:")
    lines.append("    path: \"data/hs37d5.fa.gz\"")
    lines.append("    url: \"https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/hs37d5.fa.gz\"")
    lines.append("    checksum: \"d10eebe06c0dbbcb04253e3294d63efc\"")
    lines.append("  grch38:")
    lines.append("    path: \"data/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz\"")
    lines.append("    url: \"https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz\"")
    lines.append("    checksum: \"939ce19062d1462c09b88c55faca4d76\"")
    lines.append("  chm13:")
    lines.append("    path: \"data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz\"")
    lines.append("    url: \"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz\"")
    lines.append("    checksum: \"9e6bf6b586bc8954208d1cc1d5f2fc99\"")
    
    # Outputs section
    lines.append("")
    lines.append("outputs:")
    lines.append("  sv_len: \"sv_len.tsv\"")
    lines.append("  rtg_stats_dir: \"results/vcfstats\"")
    lines.append("  subset_vcf_dir: \"results/subset_vcfs\"")
    
    # Write to file
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))
        f.write('\n')
    
    print(f"Config written to {output_path}")


if __name__ == "__main__":
    output = sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"
    merge_configs(output)
