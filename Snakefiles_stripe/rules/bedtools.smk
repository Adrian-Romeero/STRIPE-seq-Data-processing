def get_strand(wildcards):
    if wildcards.strand == "forward":
        return "+"
    else:
        return "-"

ruleorder: get_bigwigs > expand_bigwigs

rule get_bigwigs:
    input:
        "results/umi_tools/{sample}.deduplicated.bam"
    output:
        "results/bedtools/single/{sample}.{strand}.bg"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/bedtools/{sample}.{strand}.log"
    params:
        strand=get_strand

    wildcard_constraints:
        sample="[a-z0-9]+_rep[123]"
    threads: 4
    shell:
        (
            "bedtools genomecov -ibam {input} -bga -5 -strand {params.strand} "
            " |sort -k 1,1 -k 2,2n > {output} 2> {log}"
        )

ruleorder: convert_to_bigwigs > convert_to_bigwigs_expanded
rule convert_to_bigwigs:
    input:
        bg="results/bedtools/single/{sample}.{strand}.bg",
        sizes="resources/reference/TAIR10.chr_sizes"
    wildcard_constraints: 
        strand=r"(forward|reverse)"
    output:
        "results/bedtools/single/{sample}.{strand}.bw"
    log:
        "logs/bedGraphToBigWig/{sample}.{strand}.log"
    conda:
        "../envs/bedgraphtobigwig.yml"
    threads: 4
    shell:
        "bedGraphToBigWig {input.bg} {input.sizes} {output} 2> {log}"

rule expand_bigwigs:
    input:
        "results/bedtools/single/{sample}.{strand}.bg"
    wildcard_constraints:
        strand=r"(forward|reverse)"
    output:
        "results/bedtools/single/{sample}.{strand}.expanded.bg"
    log:
        "logs/bedtools/{sample}.{strand}.expanded.log"
    shell:
        " python3 workflow/scripts/expand_bedgraph.py -i {input} -o {output} 2> {log}"

def get_condition(wildcards):
    samples = SAMPLES.loc[SAMPLES.Condition == wildcards.condition, "Name"]

    return expand(
        "results/bedtools/single/{sample}.{strand}.expanded.bg",
        sample=samples,
        strand=wildcards.strand
    )

use rule convert_to_bigwigs as convert_to_bigwigs_expanded with:
    input:
        bg="results/bedtools/single/{sample}.{strand}.expanded.bg",
        sizes="resources/reference/TAIR10.chr_sizes"
    output:
        "results/bedtools/single/{sample}.{strand}.expanded.bw"
    log:
        "logs/bedGraphToBigWig/{sample}.{strand}.expanded.log"

rule join_bigwigs:
    input:
        get_condition
    wildcard_constraints:
        strand =r"(forward|reverse)"
    output:
        "results/bedtools/join/{condition}.{strand}.bg"
    log:
        "logs/wiggle_tools/{condition}.{strand}.log"
    conda:
        "../envs/wiggle_tools.yml"
    threads: 8
    shell:
        "(wiggletools mean {input} | wiggletools write_bg {output} -) &> {log}"


use rule convert_to_bigwigs as convert_to_bigwigs_join with:
    input:
        bg="results/bedtools/join/{condition}.{strand}.bg",
        sizes="resources/reference/TAIR10.chr_sizes"
    output:
        "results/bedtools/join/{condition}.{strand}.bw"
    log:
        "logs/bedGraphToBigWig/{condition}.{strand}.log"
    wildcard_constraints:
        strand=["forward","reverse"]

rule combine_tracks:
    input:
        "results/bedtools/single/{sample}.forward.expanded.bw",
        "results/bedtools/single/{sample}.reverse.expanded.bw"
    output:
        "results/bedtools/join/{sample}.bg"
    log:
        "logs/wiggle_tools/{sample}.log"
    conda:
        "../envs/wiggle_tools.yml"
    threads: 8
    shell:
        "(wiggletools sum {input} | wiggletools write_bg {output} -) &> {log}"

use rule convert_to_bigwigs as convert_to_bw_joined_tracks with:
    input:
        bg="results/bedtools/join/{sample}.bg",
        sizes="resources/reference/TAIR10.chr_sizes"
    output:
        "results/bedtools/join/{sample}.bw"
    log:
        "logs/bedGraphToBigWig/{sample}.log"

rule log_convert:
    input:
        "results/bedtools/join/means/{sample}_mean.bw"
    output:
        "results/bedtools/join/log2_means/{sample}_mean.log2.bg"
    conda:
        "../envs/wiggle_tools.yml"
    log:
        "logs/wiggle_tools/{sample}_mean.log"
    threads: 8
    shell:
        "(wiggletools offset 1 {input} | wiggletools log 2 - | wiggletools write_bg - - | sort -k 1,1 -k 2,2n > {output}) &> {log}"

use rule convert_to_bigwigs as convert_to_bw_log_bgs with:
    input:
        bg ="results/bedtools/join/log2_means/{sample}_mean.log2.bg",
        sizes ="resources/reference/TAIR10.chr_sizes"
    output:
        "results/bedtools/means/log2_means/{sample}_mean.log2.bw"
    log:
        "logs/bedGraphToBigwig/{sample}_mean.log"

