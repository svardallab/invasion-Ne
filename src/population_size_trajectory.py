import sys
import msprime
import demes
import numpy as np
import pandas as pd


def main(infile):
    # Read demography from demes
    graph = demes.load(infile)
    demography = msprime.Demography.from_demes(graph)
    debug = demography.debug()
    generations = np.arange(0, 600 + 1)
    ne = debug.population_size_trajectory(generations)[:, 0]
    pd.DataFrame({"generations": generations, "Ne": ne}).to_csv(
        "/dev/stdout", index=False
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python population_size_trajectory.py <input_yaml> > <output>")
        sys.exit(1)
    input_yaml = sys.argv[1]
    main(input_yaml)
