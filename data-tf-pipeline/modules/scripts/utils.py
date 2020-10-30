#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import sys

def getTargetInfo(config):
    targetFiles = []
    targetFiles.extend([
                        _getcutadapt(config),
                        _getbwa(config),
                        _getmutect2(config),
                        _getgatk4(config),
                        _getpicard(config)
                        ])
    return targetFiles


def _getcutadapt(config):
    ls = []
    for sample in config['ordered_sample_list']:
        ls.append("analysis/cutadapt/%s/%s_R1.fastq.gz" %(sample, sample))
        ls.append("analysis/cutadapt/%s/%s_R2.fastq.gz" %(sample, sample))
    return ls

def _getbwa(config):
    ls = []
    for sample in config['ordered_sample_list']:
        ls.append("analysis/BWA/%s/%s.bam" %(sample, sample))
        ls.append("analysis/BWA/%s/%s.sort.bam" %(sample, sample))
        ls.append("analysis/BWA/%s/%s.dedup.bam" %(sample, sample))
    return ls


def _getmutect2(config):
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

def _getpicard(config):
    ls = []
    for sample in config['ordered_sample_list']:
        ls.append("analysis/picard/%s/%s.dedup.HSmetrics.txt" %(sample, sample))
        ls.append("analysis/picard/%s/%s.dedup.HSmetrics_pertarget.txt" %(sample, sample))
        ls.append("analysis/statistics/%s/%s.coverage.cov" %(sample, sample)),
    return ls

def _getgatk4(config):
    ls = []
    for sample in config['ordered_sample_list']:
        ls.append("analysis/gatk4/%s/%s.recaltable" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.bam" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.haplotypecall.g.vcf" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.haplotypecall.vcf" %(sample, sample))
        ls.append("analysis/gatk4/%s/%s.recal.haplotypecall.final.vcf" %(sample, sample))
    return ls
