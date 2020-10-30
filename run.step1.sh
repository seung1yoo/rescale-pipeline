#!/bin/bash

mkdir analysis
mkdir logs

snakemake \
	-j 150 \
    -s data-tf-pipeline/wes_pipe.step1.snakefile \
	-p \
	--local-cores 8 \
	--latency-wait=90 \
	--max-jobs-per-second 5 \
	--notemp \
	--keep-target-files \
	--keep-going \
	--cluster-config data-tf-pipeline/cluster.json \
	--cluster "qsub -V -o {config[workdir]}/{cluster.output} -e {config[workdir]}/{cluster.error} -pe smp {cluster.threads} -N {cluster.jobName} -S /bin/bash"

