
chr_sizes = {}

with open("resources/reference/TAIR10_chr_all.fas", "r") as fasta:
    for line in fasta:
        if line.startswith(">"):
            chr_name = line.split(" ")[0].strip(">")
            chr_sizes[chr_name] = 0
        else:
            bases = len(line.strip("\n"))
            chr_sizes[chr_name] += bases

with open("resources/reference/TAIR10.chr_sizes", "w") as out:
    for chrom, size in chr_sizes.items():
        out.write(chrom + "\t" + str(size) + "\n")
