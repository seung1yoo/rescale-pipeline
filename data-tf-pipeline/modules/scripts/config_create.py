#!/usr/bin/env python

import sys
import os
import glob
import argparse
from collections import Counter


# target_bed_name="YONSEI"
#target_bed_name="TMB356"


def write_target_bed_file(output_fp, target_bed_name, project_path):
    if target_bed_name=="V5":
        print("target_bed:         %s/References/homo_sapiens/hg19/targetkit/SureSelect_Human_All_Exon_V5.bed" %project_path, file=output_fp)
        print("interval_list:      %s/References/homo_sapiens/hg19/targetkit/SureSelect_Human_All_Exon_V5.interval_list" %project_path, file=output_fp)
        print("sv_target_bed:      %s/References/homo_sapiens/hg19/targetkit/SureSelect_Human_All_Exon_V5.bed" %project_path, file=output_fp)
        print("cnv_target_bed:     %s/References/homo_sapiens/hg19/targetkit/SureSelect_Human_All_Exon_V5.bed" %project_path, file=output_fp)
        print("variant_target_bed: %s/References/homo_sapiens/hg19/targetkit/SureSelect_Human_All_Exon_V5.bed" %project_path, file=output_fp)
    else:
        print("target_bed:     ", file=output_fp)
        print("interval_list:  ", file=output_fp)
        print("sv_target_bed:  ", file=output_fp)
        print("target_tran:    ", file=output_fp)
    return

def create_configure(project_path, sample_file, target_bed_name):
    input_rawdata_path = "%s/rawdata" %project_path
    output_rawdata_path = "%s/rawdata_merge" %project_path
    os.system("mkdir -p %s" %output_rawdata_path)

    dic_sample_id = {}
    with open(sample_file) as input_fp:
        for line in input_fp:
            units = line.strip().split()
            sample_id = units[0]
            delivery_id = units[1]
            dic_sample_id[sample_id] = delivery_id

    list_fastq_file = glob.glob("%s/*/*.gz" %(input_rawdata_path))
    list_fastq_file.sort()
    if len(list_fastq_file) == 0:
        list_fastq_file = glob.glob("%s/*.gz" %(input_rawdata_path))

    set_fastq_sample =set()
    total_fastq_file = {}
    for fastq_file in list_fastq_file:
        fastq_file_name = os.path.basename(fastq_file)
        sample_id = fastq_file_name.split("--")[0]
        if  sample_id in dic_sample_id:
            print("KNOWN : %s" %fastq_file)
        else:
            sample_id = fastq_file_name.split("_")[0]
            if sample_id in dic_sample_id:
                print("KNOWN : %s" %fastq_file)
            else:
                sample_id = fastq_file_name.split("_R")[0]
                if sample_id in dic_sample_id:
                    print("KNOWN : %s" %fastq_file)
                else:
                    print("UNKNOWN : %s" %fastq_file)
                    continue

        set_fastq_sample.add(sample_id)

        if sample_id in total_fastq_file:
            if "_R1" in fastq_file_name:
                total_fastq_file[sample_id]["R1"].append(fastq_file)
            elif "_R2" in fastq_file_name:
                total_fastq_file[sample_id]["R2"].append(fastq_file)
            else:
                if "_1" in fastq_file_name:
                    total_fastq_file[sample_id]["R1"].append(fastq_file)
                if "_2" in fastq_file_name:
                    total_fastq_file[sample_id]["R2"].append(fastq_file)
        else:
            total_fastq_file[sample_id] = {}
            total_fastq_file[sample_id]["R1"] = []
            total_fastq_file[sample_id]["R2"] = []
            if "_R1" in fastq_file_name:
                total_fastq_file[sample_id]["R1"].append(fastq_file)
            elif "_R2" in fastq_file_name:
                total_fastq_file[sample_id]["R2"].append(fastq_file)
            else:
                if "_1" in fastq_file_name:
                    total_fastq_file[sample_id]["R1"].append(fastq_file)
                if "_2" in fastq_file_name:
                    total_fastq_file[sample_id]["R2"].append(fastq_file)


    # sample ID check
    print("========================================================================================================================================================")
    print("========================================================================================================================================================")
    print("# Sample ID check")
    for sample_id in dic_sample_id:
        if sample_id in set_fastq_sample:
            print("KNOWN : %s" %sample_id)
        else:
            print("UNKNOWN : %s" %sample_id)
    print("========================================================================================================================================================")
    print("========================================================================================================================================================")


    # Merge fastq file
    os.system("mkdir -p %s" %output_rawdata_path)

    output_fp = open("%s/merge.sh" %project_path, 'w')
    for sample_id in total_fastq_file.keys():
        output_R1_file = '%s/%s_R1.fastq.gz' %(output_rawdata_path, sample_id)
        output_R2_file = '%s/%s_R2.fastq.gz' %(output_rawdata_path, sample_id)
        if len(total_fastq_file[sample_id]["R1"]) == 1:
            print('ln -s %s %s'  %(' '.join(total_fastq_file[sample_id]["R1"]), output_R1_file), file=output_fp)
            print('ln -s %s %s'  %(' '.join(total_fastq_file[sample_id]["R2"]), output_R2_file), file=output_fp)
        else:
            print('cat %s > %s'  %(' '.join(total_fastq_file[sample_id]["R1"]), output_R1_file), file=output_fp)
            print('cat %s > %s'  %(' '.join(total_fastq_file[sample_id]["R2"]), output_R2_file), file=output_fp)
    output_fp.close()

    output_fp = open("%s/config.yaml" %project_path, "w")
    print("---", file=output_fp)
    print("metasheet: metasheet.csv", file=output_fp)
    print("ref:     %s/data-tf-pipeline/ref.yaml" %project_path, file=output_fp)
    print("pipedir: %s/data-tf-pipeline" %project_path, file=output_fp)
    if target_bed_name == "mm10":
        print("assembly: 'mm10'", file=output_fp)
    else:
        print("assembly: 'hg19'", file=output_fp)
    print("cluster: qsub", file=output_fp)
    print("workdir: %s" %project_path, file=output_fp)
    print("", file=output_fp)
    print("#Target Region BED file", file=output_fp)

    ## target bed file 
    write_target_bed_file(output_fp, target_bed_name, project_path)

    print("", file=output_fp)
    print("pvacseq: False", file=output_fp)
    print("", file=output_fp)
    print("#PARAMS", file=output_fp)
    print("threads: 8", file=output_fp)
    print("memory: 50", file=output_fp)
    print("bwa_seed: 18", file=output_fp)
    print("", file=output_fp)
    print("samples:", file=output_fp)
    for sample_id in dic_sample_id:
        delivery_id = dic_sample_id[sample_id]
        output_R1_file = '%s/%s_R1.fastq.gz' %(output_rawdata_path, sample_id)
        output_R2_file = '%s/%s_R2.fastq.gz' %(output_rawdata_path, sample_id)
        print("  %s:" %delivery_id, file=output_fp)
        print("    fastq_1: %s" %output_R1_file, file=output_fp)
        print("    fastq_2: %s" %output_R2_file, file=output_fp)
    output_fp.close()

    output_fp = open("%s/metasheet.csv" %project_path, "w")
    print("comp_id,normal,tumor", file=output_fp)
#     for sample_id in dic_sample_id:
#         delivery_id = dic_sample_id[sample_id]
#         print("%s" %delivery_id, file=output_fp)
    output_fp.close()

    output_fp = open("%s/run.sh" %project_path, "w")
    print("#!/bin/bash", file=output_fp)
    print("", file=output_fp)
    print("mkdir %s/analysis" %project_path, file=output_fp)
    print("mkdir %s/logs" %project_path, file=output_fp)
    print("", file=output_fp)
    print("snakemake \\", file=output_fp)
    print("\t-j 150 \\", file=output_fp)
    print("\t-s %s/data-tf-pipeline/wes_pipe.snakemake \\" %project_path, file=output_fp)
    print("\t-p \\", file=output_fp)
    print("\t--local-cores 8 \\", file=output_fp)
    print("\t--latency-wait=90 \\", file=output_fp)
    print("\t--max-jobs-per-second 5 \\", file=output_fp)
    print("\t--notemp \\", file=output_fp)
    print("\t--keep-target-files \\", file=output_fp)
    print("\t--keep-going \\", file=output_fp)
    print("\t--cluster-config %s/data-tf-pipeline/cluster.json \\" %project_path, file=output_fp)
    print("\t--cluster \"qsub -V -o {config[workdir]}/{cluster.output} -e {config[workdir]}/{cluster.error} -pe smp {cluster.threads} -N {cluster.jobName} -S /bin/bash\"", file=output_fp)
    output_fp.close()
    return

def usage():
    message='''
python %s  --project_path  /TBI/Share/BioPeople/shsong/Projects/WES/BioProjects/TBD190181-MedPacto-20190910-UMI -s ./sample.xls  -b V5

-i, --project_path : input data path
-s, --sample_file : sample id
-b, --bed_file_name : [V5, V6, mm10]
    ''' %sys.argv[0]
    print(message)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--project_path')
    parser.add_argument('-s', '--sample_file')
    parser.add_argument('-b', '--bedfile_name', default="V5")
    args = parser.parse_args()
    try:
        len(args.project_path) > 0
        len(args.sample_file) > 0

    except:
        usage()
        sys.exit(2)

    create_configure(args.project_path, args.sample_file, args.bedfile_name)

if __name__ == '__main__':
    main()

