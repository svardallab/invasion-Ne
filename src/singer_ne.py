#!/usr/bin/env python3

import numpy as np
import pandas as pd
import tskit
import sys


def main():
    tree_files = sys.argv[1:]
    if not tree_files:
        print("Error: No input files provided", file=sys.stderr)
        print("Usage: python singer_ne.py file1.trees file2.trees ...", file=sys.stderr)
        sys.exit(1)
    # Define time windows
    time = np.append(np.linspace(0, 600, 20), np.inf)
    # Initialize results list
    results = []
    # Process each tree file
    for tree_file in tree_files:
        try:
            # Load the tree sequence
            ts = tskit.load(tree_file)
            result = 1 / ts.pair_coalescence_rates(time_windows=time) / 2
            results.append(result)
            print(f"Processed: {tree_file}", file=sys.stderr)
        except Exception as e:
            print(f"Error processing {tree_file}: {e}", file=sys.stderr)
    if not results:
        print("No valid tree files processed", file=sys.stderr)
        sys.exit(1)
    # Convert list of results to a numpy array
    results_array = np.array(results)
    # Create dataframe with results
    df = pd.DataFrame(
        {
            "generations": time[:-1],
            "Ne": results_array.mean(axis=0),
            "Ne_lower90": np.quantile(results_array, 0.05, axis=0),
            "Ne_upper90": np.quantile(results_array, 0.95, axis=0),
        }
    )
    df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
