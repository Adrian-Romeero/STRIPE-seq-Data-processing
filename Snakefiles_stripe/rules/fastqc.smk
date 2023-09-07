rule fastqc_raw:
    input:
        "resources/fastq/{sample}.fastq.gz"
    output:
        "results/qc/fastqc_raw/{sample}_fastqc.zip",
        "results/qc/fastqc_raw/{sample}_fastqc.html"
    params:
        outdir=lambda wildcards, output: path.dirname(output[0])
    log:
        "logs/fastqc/{sample}.log"
    conda:
        "../envs/fastqc.yml"
    threads: 8
    shell:
        "fastqc -o {params.outdir} {input} &> {log}"

use rule fastqc_raw as fastqc_umi with:
    input:
        "results/umi_tools/{sample}.umi_tools.fastq.gz"
    output:
        "results/qc/qc_umi/{sample}.umi_tools_fastqc.zip",
        "results/qc/qc_umi/{sample}.umi_tools_fastqc.html"
    log:
        "logs/fastqc/{sample}.umi.log"

use rule fastqc_raw as qc_first_cutadapt with:
    input:
        "results/cutadapt/{sample}.cutadapt.fastq.gz"
    output:
        "results/qc/qc_cutadapt/{sample}.cutadapt_fastqc.zip",
        "results/qc/qc_cutadapt/{sample}.cutadapt_fastqc.html"
    log:
        "logs/fastqc/{sample}.first_cutadapt.log"

use rule fastqc_raw as fastq_second_pass with:
    input:
        "results/cutadapt_second_pass/{sample}.template_switch.fastq.gz"
    output:
        "results/qc/fastqc_curated/{sample}.template_switch_fastqc.html",
        "results/qc/fastqc_curated/{sample}.template_switch_fastqc.zip"
    log:
        "logs/fastqc/{sample}.curated.log"
