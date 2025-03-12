#!/usr/bin/env python
import tskit
import tszip
import msprime
import gzip
import sys


def main(infile, seed):
    # Read tszip file
    ts = tszip.decompress(infile)

    # Overlay mutations
    mts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)

    # Convert to VCF
    n_dip_indv = int(mts.num_samples / 2)
    indv_names = [f"tsk_{i}indv" for i in range(n_dip_indv)]
    mts.write_vcf(sys.stdout, individual_names=indv_names)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.ts.zip> <seed>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
