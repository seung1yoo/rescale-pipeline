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
        ls.append("analysis/BWA/%s/%s.bam" %(sample, sample))
        ls.append("analysis/BWA/%s/%s.sort.bam" %(sample, sample))
        ls.append("analysis/BWA/%s/%s.dedup.bam" %(sample, sample))
    return ls

rule target:
    input: getTargetInfo(config)
    message: "Compiling all output"

rule run_BWA:
    input:
        fastq_1="analysis/cutadapt/{sample}/{sample}_R1.fastq.gz",
        fastq_2="analysis/cutadapt/{sample}/{sample}_R2.fastq.gz"
    output:
        bam="analysis/BWA/{sample}/{sample}.bam",
    wildcard_constraints:
          sample="[^/]+"
    params:
        seed=config["bwa_seed"],
        tmp="analysis/BWA/{sample}/temp",
        readgroup=lambda wildcards: "\'@RG\\tID:{sample}\\tPL:illumina\\tSM:{sample}\'".format(sample=wildcards.sample)
    threads: 8
    message: "Running BAM Alignment on {wildcards.sample}"
    priority: 10
    benchmark:
        "benchmarks/{sample}/{sample}.run_BWA.txt"
#     log:
#         err="logs/BWA/{sample}/{sample}.err.log",
#         out="logs/BWA/{sample}/{sample}.out.log"
    shell:
        "mkdir -p {params.tmp} && "
        "bwa mem -M"
        " -R {params.readgroup}"
        " -t {threads}"
        " -k {params.seed}"
        " {config[ref_fasta]}"
        " {input.fastq_1} {input.fastq_2} |" 
        " samtools view -Sb -F 0x100 > {output.bam} "
#         " 2> {log.err} 1> {log.out}"

rule sortbam:
    input:
        bam="analysis/BWA/{sample}/{sample}.bam",
    output:
        bam="analysis/BWA/{sample}/{sample}.sort.bam",
    wildcard_constraints:
          sample="[^/]+"
    params:
        tmp="analysis/BWA/{sample}/temp",
    threads: 8
    message: "Running BAM Alignment on {wildcards.sample}"
    priority: 10
    benchmark:
        "benchmarks/{sample}/{sample}.run_BWA.txt"
#     log:
#         err="logs/BWA/{sample}/{sample}.err.log",
#         out="logs/BWA/{sample}/{sample}.out.log"
    shell:
        "mkdir -p {params.tmp} && "
        "picard SortSam"
        " TMP_DIR={params.tmp}"
        " I={input.bam}"
        " O={output.bam}"
        " SO=coordinate"


rule run_rmdup:
    input:
        bam="analysis/BWA/{sample}/{sample}.sort.bam",
    output:
        bam="analysis/BWA/{sample}/{sample}.dedup.bam",
        metrics="analysis/BWA/{sample}/{sample}.dedup.txt",
    wildcard_constraints:
          sample="[^/]+"
    params:
        tmp="analysis/BWA/{sample}/temp",
    threads: 8
    message: "Running dedup on {wildcards.sample}"
    priority: 10
    benchmark:
        "benchmarks/{sample}/{sample}.run_dedup.txt"
    shell:
        "picard MarkDuplicates"
        " -Xmx50g"
        " TMP_DIR={params.tmp}"
        " I={input.bam}"
        " O={output.bam}"
        " M={output.metrics}"
        " CREATE_INDEX=true"
        " CREATE_MD5_FILE=true"
        " REMOVE_SEQUENCING_DUPLICATES=true"

