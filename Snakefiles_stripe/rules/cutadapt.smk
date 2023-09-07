def get_cutadapt_options(wildcards):
    if not PAIRED_END:
        return "-g ^TATAG{3} -O 6"
    else:
        return "-g ^N{8}TATAG{3} -O 14"

def get_cutadapt_input(wildcards):
    if PAIRED_END:
        pass
    else:
        return "results/umi_tools/{sample}.umi_tools.fastq.gz"

rule cutadapt:
    input:
        "results/umi_tools/{sample}.umi_tools.fastq.gz"
    output:
        "results/cutadapt/{sample}.cutadapt.fastq.gz"
    params:
        options=get_cutadapt_options
    conda:
        "../envs/cutadapt.yml"
    threads: 8
    log:
        "logs/cutadapt/{sample}.log"
    shell:
        (
            "cutadapt -j {threads} "
            "{params.options} -e 1 --discard-untrimmed "
            "-o {output} {input}  &> {log}"
        )

rule cutadapt_second_pass:
    input:
        rules.cutadapt.output
    output:
        "results/cutadapt_second_pass/{sample}.template_switch.fastq.gz"
    params:
        sequence=lambda wildcards: '^CCTACACGACGCTCTTCCGATCTN{8}TATAG{3}',
    conda:
        "../envs/cutadapt.yml"
    threads: 8
    log:
        "logs/cutadapt/{sample}.template_switch.log"
    shell:
        (
            "cutadapt -g {params.sequence} "
            "-j {threads} "
            "-O 21 -e 2 -m 15 --discard-untrimmed "
            "-o {output} {input} &> {log}"
        )
