# A plink-friendly script to write VCF files from tree-sequences.
import tskit
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: tskit_vcf_export.py <INFILE> <contig> > <OUTPUT_FILE>")
        sys.exit(1)
    infile = sys.argv[1]
    contig_id = sys.argv[2]
    ts = tskit.load(infile)
    n_dip_indv = int(ts.num_samples / 2)
    indv_names = [f"tsk_{i}indv" for i in range(n_dip_indv)]
    with open("/dev/stdout", "w") as vcf_file:
        ts.write_vcf(
            vcf_file,
            individual_names=indv_names,
            contig_id=contig_id,
            allow_position_zero=True,
        )
