#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from modules.scripts.config_setup import updateConfig
from modules.scripts.metasheet_setup import updateMeta

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
config = updateConfig(config)
config = updateMeta(config)
#-----------------------------------------

def getTargetInfo(config):
    ls = []
    for sample in config['ordered_sample_list']:
        ls.append("analysis/gatk4/%s/%s.recaltable" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.bam" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.haplotypecall.g.vcf" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.haplotypecall.vcf" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.haplotypecall.final.vcf" %(sample, sample))
    return ls

rule target:
    input: getTargetInfo(config)
    message: "Compiling all output"



rule gatk_bqsr:
    input:
#         bam="analysis/picard/{sample}/{sample}.dedup.bam",
        bam="analysis/BWA/{sample}/{sample}.dedup.bam",
    output:
        recal_table="analysis/gatk4/{sample}/{sample}.recaltable"
#     wildcard_constraints:
#           sample="[^/]+"
#     log:
#         "logs/gatk/bqsr/{sample}.log"
    params:
        extra="",  # optional
        java_opts="", # optional
        tmp=temp("analysis/gatk4/{sample}/temp"),
    shell:
        "mkdir -p {params.tmp} && "
        "gatk  BaseRecalibrator"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_fasta]}"
        " -I {input.bam}"
        " -O {output.recal_table}"
        " --known-sites {config[dbsnp_gatk]}"


rule gatk_bqsr_apply:
    input:
        bam="analysis/BWA/{sample}/{sample}.dedup.bam",
        recal_table="analysis/gatk4/{sample}/{sample}.recaltable"
    output:
        bam="analysis/gatk4/{sample}/{sample}.recal.bam",
#     wildcard_constraints:
#           sample="[^/]+"
#     log:
#         "logs/gatk/bqsr/{sample}.log"
    params:
        extra="",  # optional
        java_opts="", # optional
        tmp=temp("analysis/gatk4/{sample}/temp"),
    shell:
        "mkdir -p {params.tmp} && "
        "gatk ApplyBQSR"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_fasta]}"
        " -I {input.bam} "
        " --bqsr-recal-file {input.recal_table}"
        " -O {output.bam}"


rule haplotype_caller_GVCF:
    input:
        bam="analysis/gatk4/{sample}/{sample}.recal.bam",
    output:
        gvcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.g.vcf",
#     wildcard_constraints:
#           sample="[^/]+"
    params:
        extra="",  # optional
        java_opts="", # optional
        tmp=temp("analysis/gatk4/{sample}/temp"),
    shell:
        "mkdir -p {params.tmp} && "
        "gatk HaplotypeCaller"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_fasta]}"
        " -I {input.bam} "
        " -O {output.gvcf}"
        " --intervals {config[variant_target_bed]}"
        " -ERC GVCF"
        " --dbsnp {config[dbsnp_gatk]}"


rule haplotype_caller_VCF:
    input:
        gvcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.g.vcf",
    output:
        vcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.vcf",
#     wildcard_constraints:
#           sample="[^/]+"
    params:
        extra="",  # optional
        java_opts="", # optional
        tmp=temp("analysis/gatk4/{sample}/temp"),
    shell:
        "mkdir -p {params.tmp} && "
        "gatk GenotypeGVCFs"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_fasta]}"
        " --intervals {config[variant_target_bed]}"
        " -V {input.gvcf} "
        " -O {output.vcf}"

rule gatk_filter:
    input:
        vcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.vcf",
    output:
        finalvcf="analysis/gatk4/{sample}/{sample}.recal.haplotypecall.final.vcf",
    params:
        tmp=temp("analysis/gatk4/{sample}/temp"),
#     wildcard_constraints:
#           sample="[^/]+"
#     log:
#         "logs/gatk/filter/snvs.log"
#     params:
#         filters={"snp_hardfilters": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", 
#                 "indel_hardfilters":"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"}
    shell:
        "mkdir -p {params.tmp} && "
        "gatk VariantFiltration"
        " --tmp-dir {params.tmp}"
        " -V {input.vcf} "
        " -O {output.finalvcf}"

