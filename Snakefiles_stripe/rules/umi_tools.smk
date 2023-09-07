def get_fastq(wildcards):
    read1 = SAMPLES[SAMPLES.Name == wildcards.sample].Read1.values[0]

    return "resources/fastq/{sample}.fastq.gz".format(
        sample=read1
    )

rule umi_tools:
    input:
        get_fastq
    output:
        "results/umi_tools/{sample}.umi_tools.fastq.gz"
    params:
        umi = "NNNNNNNN"
    log:
        "logs/umi_tools/{sample}.log"
    conda:
        "../envs/umi_tools.yml"
    threads: 8
    shell:
        "umi_tools extract -I {input} -S {output} -p {params.umi} -L {log} -E {log}"

rule umi_tools_second_pass:
    input:
        bam="results/samtools/{sample}.out.bam",
        bai="results/samtools/{sample}.out.bam.bai"
    output:
        "results/umi_tools/{sample}.deduplicated.bam"
    log:
        "logs/umi_tools/{sample}.dedup.log"
    conda:
        "../envs/umi_tools.yml"
    threads: 8
    shell:
        "umi_tools dedup -I {input.bam} -S {output} &> {log}"
        