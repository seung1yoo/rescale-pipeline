#!/usr/bin/env python


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
