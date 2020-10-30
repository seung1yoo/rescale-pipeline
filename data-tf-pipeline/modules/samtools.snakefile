#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab


rule mpileup:
    input:
        bam="analysis/gatk4/{sample}/{sample}.recal.bam",
    output:
        mpileup="analysis/samtools/{sample}/{sample}.mpileup.txt",
    wildcard_constraints:
          sample="[^/]+"
#     params:
#         frag_len=config["frag_len"]
    threads: 4
    message: "Running samtools mpileup on {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.run_mpileup.txt"
    shell:
        "samtools mpileup -d100000 -C50 -l {config[variant_target_bed]} -f {config[ref_fasta]} {input.bam} >  {output.mpileup}"
