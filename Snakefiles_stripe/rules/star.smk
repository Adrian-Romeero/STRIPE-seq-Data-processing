rule index:
    input:
        fas ="resources/reference/TAIR10_chr_all.fas",
        gtf="resources/reference/Araport11_GTF_genes.gffread.gtf"
    output:
        directory("resources/reference/TAIR10_STAR_index/")

    conda:
        "../envs/star.yml"
    threads: 4
    shell:
        (
        "STAR --runThreadN {threads}"
        "--runMode genomeGenerate"
        "--genomeSAindexNbases 12"
        "--genomeDir {output}"
        "--genomeFastaFiles {input.fas}"
        "--sjdbGTFfile {input.gtf}"
        " --sjdbOverhang 60"
        )

rule align:
    input:
        rules.cutadapt.output,
        refdir = "resources/reference/TAIR10_STAR_index/"
    output:
        sam="results/star/{sample}/Aligned.out.sam"
    log:
        "logs/star/{sample}.log"
    params:
        outdir=lambda wildcards, output: path.dirname(output.sam)
    conda:
        "../envs/star.yml"
    threads: 4 
    shell:
        (
            "STAR --runThreadN {threads}"
            "--genomeDir {input.refdir}"
            "--readFilesIn {input}"
            "--readFilesCommand gunzip -c"
            " --outFileNamePrefix {params.outdir} "
            "--outSAMtype SAM &> {log}"
        )