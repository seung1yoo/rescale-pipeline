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
        ls.append("analysis/picard/%s/%s.dedup.HSmetrics.txt" %(sample, sample))
        ls.append("analysis/picard/%s/%s.dedup.HSmetrics_pertarget.txt" %(sample, sample))
        #ls.append("analysis/picard/%s/%s.insertsize.txt" %(sample, sample))
        #ls.append("analysis/picard/%s/%s.insertsize.pdf" %(sample, sample))
        ls.append("analysis/statistics/%s/%s.coverage.cov" %(sample, sample)),
    return ls

rule target:
    input: getTargetInfo(config)
    message: "Compiling all output"

rule CollectHsMetrics:
    input:
        bam="analysis/BWA/{sample}/{sample}.dedup.bam",
    output:
        metrics="analysis/picard/{sample}/{sample}.dedup.HSmetrics.txt",
        per_t_metrics="analysis/picard/{sample}/{sample}.dedup.HSmetrics_pertarget.txt",
        cov="analysis/statistics/{sample}/{sample}.coverage.cov",
#     wildcard_constraints:
#         sample = "{sample}".split("/")[0]
#     log:
#         err="logs/picard/annotate_metrics/{sample}.err.log",
#         out="logs/picard/annotate_metrics/{sample}.out.log"
    shell:
        "picard CollectHsMetrics "
        " -Xmx50g"
        " I={input.bam} O={output.metrics} "
        " R={config[ref_fasta]} "
        " TARGET_INTERVALS={config[interval_list]} "
        " BAIT_INTERVALS={config[interval_list]} "
        " PER_TARGET_COVERAGE={output.per_t_metrics} "
        " COVERAGE_CAP=100000 && "
        "cp {output.per_t_metrics} {output.cov} "
#         " 2> {log.err} 1> {log.out}"


rule CollecInsertsize:
    input:
        bam="analysis/BWA/{sample}/{sample}.dedup.bam",
    output:
        metrics="analysis/picard/{sample}/{sample}.insertsize.txt",
        pdf="analysis/picard/{sample}/{sample}.insertsize.pdf",
#     wildcard_constraints:
#         sample = "{sample}".split("/")[0]
#     log:
#         err="logs/picard/filtercons_hs_metrics/{sample}.err.log",
#         out="logs/picard/filtercons_hs_metrics/{sample}.out.log"
    shell:
        "picard CollectInsertSizeMetrics "
        " -Xmx50g"
        " I={input.bam}"
        " O={output.metrics}"
        " H={output.pdf}"
        " M=0.5"
