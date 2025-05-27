import tskit
import numpy as np
import pandas as pd
import sys
def analyze(infile: str, recombination_rate: float):
    ts = tskit.load(infile)
    return ts.diversity()

def main(infiles, recombination_rate, mutation_rate):
    # Get a ballpark estimate of Ne in windows and concatenate
    diversities = np.array([
        analyze(infile, recombination_rate=recombination_rate) for infile in infiles
    ])
    ne = diversities / 4 / mutation_rate
    df = pd.DataFrame({'Ne': ne, 'File': infiles})
    df.to_csv('/dev/stdout', index=False)

if __name__ == "__main__":
    mutation_rate = float(sys.argv[1])
    recombination_rate = float(sys.argv[2])
    infiles = sys.argv[3:]
    main(infiles, recombination_rate=recombination_rate, mutation_rate=mutation_rate)
