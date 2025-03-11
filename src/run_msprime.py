import sys
import msprime
import tskit
import demesdraw
import demes
import numpy as np
import matplotlib.pyplot as plt
import tszip


def main(infile, seed, outfile):
    # Read demography from demes
    graph = demes.load(infile)
    demography = msprime.Demography.from_demes(graph)
    # Simulation
    sequence_length = 1e8
    recombination_rate = 1e-8
    ts = msprime.sim_ancestry(
        samples=100,
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=seed,
    )
    # Compress output with tszip
    tszip.compress(ts, outfile)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_yaml> <seed> <output>")
        sys.exit(1)

    input_yaml = sys.argv[1]
    seed = int(sys.argv[2])  # Ensure seed is an integer
    output = sys.argv[3]

    main(input_yaml, seed, output)
