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

    times = []
    # Initialize results list
    results = []
    # Process each tree file
    quantiles = np.linspace(0, 1, 25)
    expected_length = len(quantiles)

    for tree_file in tree_files:
        ts = tskit.load(tree_file)
        t = ts.pair_coalescence_quantiles(quantiles)

        # Check if t has the expected length
        if len(t) != expected_length:
            print(
                f"Warning: {tree_file} has {len(t)} quantiles instead of {expected_length}. Padding with np.inf.",
                file=sys.stderr,
            )
            # Pad with np.inf to reach expected length
            padded_t = np.full(expected_length, np.inf)
            padded_t[: len(t)] = t
            t = padded_t
        t[0] = 0.0
        result = 1 / ts.pair_coalescence_rates(time_windows=np.append(t, np.inf)) / 2
        times.append(t)
        results.append(result)
        print(f"Processed: {tree_file}", file=sys.stderr)

    # Convert list of results to a numpy array
    times_array = np.array(times)
    results_array = np.array(results)

    # Create dataframe with results
    df = pd.DataFrame(
        {
            "generations": times_array.mean(axis=0),
            "Ne": results_array.mean(axis=0),
            "Ne_lower90": np.quantile(results_array, 0.05, axis=0),
            "Ne_upper90": np.quantile(results_array, 0.95, axis=0),
        }
    )
    df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
