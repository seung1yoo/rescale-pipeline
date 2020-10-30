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
        ls.append("analysis/mutect2/%s/%s.gatk4.metrics.pre_adapter_detail_metrics" %(sample, sample))
        if config['assembly'] == "hg19" or config['assembly'] == "hg38":
            ls.append("analysis/mutect2/%s/%s.gatk4.table" %(sample, sample))
            ls.append("analysis/mutect2/%s/%s.gatk4.contamination.table" %(sample, sample))
        ls.append("analysis/mutect2/%s/%s.mutect2.vcf" %(sample, sample))
        ls.append("analysis/mutect2/%s/%s.mutect2.filt.vcf" %(sample, sample))
        ls.append("analysis/mutect2/%s/%s.mutect2.filt.OxoG.vcf" %(sample, sample))
        ls.append("analysis/mutect2/%s/%s.mutect2.final.vcf" %(sample, sample))
    return ls

rule target:
    input: getTargetInfo(config)
    message: "Compiling all output"


rule gatk4_metrics:
    input:
        bam="analysis/gatk4/{sample}/{sample}.recal.bam",
    output:
        metrics="analysis/mutect2/{sample}/{sample}.gatk4.metrics",
        adapter_metrics="analysis/mutect2/{sample}/{sample}.gatk4.metrics.pre_adapter_detail_metrics",
    message: "Running collectSequencing Artifact Metrics on {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.gatk4_metrics.txt"
    shell:
        "picard CollectSequencingArtifactMetrics"
        " I={input.bam}"
        " O={output.metrics}"
        " R={config[ref_fasta]} && "
        "touch {output.metrics}"


rule getpileupsummary:
    input:
        bam="analysis/gatk4/{sample}/{sample}.recal.bam",
    output:
        table="analysis/mutect2/{sample}/{sample}.gatk4.table",
    benchmark:
        "benchmarks/{sample}/{sample}.getpileupsummary.txt"
    shell:
        "gatk GetPileupSummaries"
        " -I {input.bam}"
        " -O {output.table}"
        " -V {config[gnomad_mutect2]}"
        " -R {config[ref_fasta]}"
        " -L {config[variant_target_bed]}"

rule CalcContamination:
    input:
        table="analysis/mutect2/{sample}/{sample}.gatk4.table",
    output:
        table="analysis/mutect2/{sample}/{sample}.gatk4.contamination.table",
    benchmark:
        "benchmarks/{sample}/{sample}.CalcContamination.txt"
    params:
        tmp="analysis/mutect2/{sample}/temp",
    shell:
        "mkdir -p {params.tmp} && "
        "gatk CalculateContamination"
        " -I {input.table}"
        " -O {output.table}"

rule Mutect2:
    input:
        bam="analysis/gatk4/{sample}/{sample}.recal.bam",
    output:
        vcf="analysis/mutect2/{sample}/{sample}.mutect2.vcf",
        bam="analysis/mutect2/{sample}/{sample}.mutect2.bam"
    params:
        tmp="analysis/mutect2/{sample}/temp",
        filt_normal_af="0.000025",
        tumorlod="2",
        sampleid="{sample}"
    message: "Running Mutect2 on {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.Mutect2.txt"
    shell:
        "mkdir -p {params.tmp} && "
        "gatk Mutect2"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_fasta]}"
        " --germline-resource {config[gnomad_mutect2]}"
        " -af-of-alleles-not-in-resource {params.filt_normal_af}"
        " -I {input.bam}"
        " -tumor {params.sampleid}"
        " --intervals {config[variant_target_bed]}"
        " --output {output.vcf}"
        " -bamout {output.bam}"
        " --tumor-lod-to-emit {params.tumorlod}"
        " --active-probability-threshold 0.001"
        " --max-reads-per-alignment-start 0"
        " --annotation AlleleFraction"
        " --annotation AS_FisherStrand"
        " --annotation TandemRepeat"

rule FilterMutect2Var:
    input:
        vcf="analysis/mutect2/{sample}/{sample}.mutect2.vcf",
#         table="analysis/mutect2/{sample}/{sample}.gatk4.contamination.table",
    output:
        filtvcf="analysis/mutect2/{sample}/{sample}.mutect2.filt.vcf",
    params:
        tmp="analysis/mutect2/{sample}/temp",
    benchmark:
        "benchmarks/{sample}/{sample}.FilterMutect2Var.txt"
    shell:
        "mkdir -p {params.tmp} && "
        "gatk FilterMutectCalls"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_fasta]}"
#         " --contamination-table {input.table}"
        " -V {input.vcf}"
        " -O {output.filtvcf}"


rule FilterMutect2Bias:
    input:
        vcf="analysis/mutect2/{sample}/{sample}.mutect2.filt.vcf",
        metrics="analysis/mutect2/{sample}/{sample}.gatk4.metrics.pre_adapter_detail_metrics",
    output:
        filtvcf="analysis/mutect2/{sample}/{sample}.mutect2.filt.OxoG.vcf",
    params:
        tmp=temp("analysis/mutect2/{sample}/temp"),
    benchmark:
        "benchmarks/{sample}/{sample}.FilterMutect2Bias.txt"
    shell:
        "mkdir -p {params.tmp} && "
        "gatk FilterByOrientationBias"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_fasta]}"
        " -P {input.metrics}"
        " -L {config[variant_target_bed]}"
        " -V {input.vcf}"
        " -O {output.filtvcf}"
        " --artifact-modes G/T"


rule FilterMutect2Final:
    input:
        vcf="analysis/mutect2/{sample}/{sample}.mutect2.filt.OxoG.vcf",
    output:
        vcf="analysis/mutect2/{sample}/{sample}.mutect2.final.vcf",
    benchmark:
        "benchmarks/{sample}/{sample}.FilterMutect2Final.txt"
    shell:
        "python {config[pipedir]}/modules/scripts/normal_variants_filter.py -i {input.vcf} -o {output.vcf}"


