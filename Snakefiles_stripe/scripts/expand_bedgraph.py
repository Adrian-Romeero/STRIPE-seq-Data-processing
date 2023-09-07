import argparse

def get_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    
    out = open(args.output, "w")

    with open(args.input) as f:
        for line in f:
            chrom, start, end, coverage = line.strip().split("\t")
            coverage = int(coverage)
            if coverage <= 1:
                continue

            width = int(end) - int(start)

            if width > 1:
                for i in range(width):
                    new_start = int(start) + i
                    new_end = new_start + 1

                    new_line = "\t".join([chrom, str(new_start), str(new_end), str(coverage)])
                    out.write(new_line + "\n")
            else:
                out.write(line)
    
    out.close()


if __name__ == "__main__":
    main()
