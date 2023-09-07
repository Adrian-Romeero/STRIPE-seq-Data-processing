rule samtools_sort_to_bam:
    input:
        "results/star/{sample}/Aligned.out.sam"
    output:
        "results/samtools/{sample}.out.bam"
    log:
        "logs/samtools/{sample}.log"
    conda:
        "../envs/samtools.yml"
    threads: 8
    shell:
        "samtools sort {input} | samtools view -F 2820 -O BAM -o {output} &> {log}"

rule samtools_index:
    input:
        "results/samtools/{sample}.out.bam"
    output:
        "results/samtools/{sample}.out.bam.bai"
    log:
        "logs/samtools/{sample}.index.log"
    conda:
        "../envs/samtools.yml"
    threads: 8
    shell:
        "samtools index {input} &> {log}"

