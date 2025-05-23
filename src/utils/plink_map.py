# This script parses the VCF header to generate a dummy map file.
import sys
import re
import os


def parse_vcf_header_from_stdin():
    contigs = []
    for line in sys.stdin:
        if line.startswith("##contig="):
            match = re.search(r"ID=([^,]+),length=(\d+)", line)
            if match:
                contigs.append((match.group(1), int(match.group(2))))
        elif line.startswith("#CHROM"):
            break
    return contigs


def main(rho):
    contigs = parse_vcf_header_from_stdin()
    for chrom, length in contigs:
        cm_length = length * rho * 100
        print(f"{chrom}\trs1\t0\t0")
        print(f"{chrom}\trs2\t{cm_length:.6f}\t{length}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: plink_map.py <RHO> < <INFILE> > <OUTPUT_FILE>")
        sys.exit(1)

    try:
        rho = float(sys.argv[1])
    except ValueError:
        print("RHO must be a numeric value.")
        sys.exit(1)

    main(rho)
