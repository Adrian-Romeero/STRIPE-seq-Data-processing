SAMPLES = pd.read_csv("workflow/schemas/metadata.csv", sep="\t")

PAIRED_END = not SAMPLES.Read2.isna().all()

def qc_raw_reads():
    return expand(
        "results/qc/fastqc_raw/{sample}_fastqc.{suffix}",
        sample=SAMPLES.Read1,
        suffix=["html", "zip"]
    )
def trim_reads():
    return expand(
        "results/cutadapt/{sample}.cutadapt.fastq.gz",
        sample=SAMPLES.Name
    )
def qc_umi_extracted_reads():
    return expand(
        "results/qc/qc_umi/{sample}.umi_tools_fastqc.{suffix}",
            sample=SAMPLES.Name,
            suffix=["html", "zip"]
    )
def qc_trimmed_reads():
    return expand(
        "results/qc/qc_cutadapt/{sample}.cutadapt_fastqc.{suffix}",
            sample=SAMPLES.Name,
            suffix=["html", "zip"]
    )

def qc_template_switch_trimmed_reads():
    if config["cutadapt"]["template_switch"]:
        return expand(
            "results/qc/fastqc_curated/{sample}.template_switch_fastqc.{suffix}",
                sample=SAMPLES.Name,
                suffix=["html", "zip"]
        )
    else:
        return []

def index_genome():
    return "resources/reference/TAIR10_STAR_index"

def align():
    return expand(
        "results/star/{sample}",
        sample=SAMPLES.Name
        )
def sort_to_bam():
    return expand(
        "results/samtools/{sample}.out.bam",
        sample=SAMPLES.Name
    )
def index_bams():
    return expand(
        "results/samtools/{sample}.out.bam.bai",
        sample=SAMPLES.Name
    )
def remove_duplicates():
    return expand(
        "results/umi_tools/{sample}.deduplicated.bam",
        sample=SAMPLES.Name
    )            
def cutadapt_second_pass():
    if config["cutadapt"]["template_switch"]:
        return expand(
            "results/cutadapt_second_pass/{sample}.template_switch.fastq.gz",
            sample=SAMPLES.Name
        )
    else:
        return []
def umi_tools():
    if PAIRED_END:
        return []
    else:
        return expand(
            "results/umi_tools/{sample}.umi_tools.fastq.gz",
            sample=SAMPLES.Name
        )
def get_bigwigs():
    return expand(
        "results/bedtools/single/{sample}.{strand}.{suffix}",
        sample=SAMPLES.Name,
        strand=["forward", "reverse"],
        suffix=["bg", "bw"]
    )
def expand_bigwigs():
    return expand(
        "results/bedtools/single/{sample}.{strand}.expanded.{suffix}",
        sample=SAMPLES.Name,
        strand=["forward","reverse"],
        suffix=["bg", "bw"]
    )

def join_bigwigs():
    return expand(
        "results/bedtools/join/{condition}.{strand}.{suffix}",
        condition=SAMPLES.Condition.unique(),
        suffix=["bg", "bw"],
        strand=["forward","reverse"]
    )
def join_tracks():
    return expand(
        "results/bedtools/join/{sample}.bw",
        sample=SAMPLES.Name
    )
def log_convert():
    return expand(
        "results/bedtools/means/log2_means/{sample}_mean.log2.bw",
        sample= SAMPLES.Condition.unique()
        )