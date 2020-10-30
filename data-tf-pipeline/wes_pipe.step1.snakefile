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
        ls.append("analysis/cutadapt/%s/%s_R1.fastq.gz" %(sample, sample))
        ls.append("analysis/cutadapt/%s/%s_R2.fastq.gz" %(sample, sample))
    return ls

rule target:
    input: getTargetInfo(config)
    message: "Compiling all output"

rule run_cutadapt:
    input:
#         fastq_1="analysis/seqtk/{sample}/{sample}_R1.fastq.gz",
#         fastq_2="analysis/seqtk/{sample}/{sample}_R2.fastq.gz"
        fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
        fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
    output:
        fastq_1="analysis/cutadapt/{sample}/{sample}_R1.fastq.gz",
        fastq_2="analysis/cutadapt/{sample}/{sample}_R2.fastq.gz",
        log_file="analysis/cutadapt/{sample}/{sample}_cutadapt.log"
#     wildcard_constraints:
#         sample = "{sample}".split("/")[0]
    params:
        basequality=20,
        minimumLen=50,
        left_adapter ="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        right_adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
#     log:
#         err="logs/cutadpat/{sample}/{sample}.err.log",
#         out="logs/cutadpat/{sample}/{sample}.out.log"
    threads: 4
    message: "Running cutadapt on {wildcards.sample}"
    priority: 10
    benchmark:
        "benchmarks/{sample}/{sample}.run_cutadapt.txt"
    shell:
        "cutadapt -q {params.basequality}"
        " -a {params.left_adapter}"
        " -A {params.right_adapter}"
        " --minimum-length {params.minimumLen}"
        " -o {output.fastq_1}"
        " -p {output.fastq_2}"
        " {input} > {output.log_file}"
#         " 2> {log.err} 1> {log.out}"
