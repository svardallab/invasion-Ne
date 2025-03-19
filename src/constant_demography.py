import sys
import msprime
import tskit
import demesdraw
import demes
import numpy as np
import matplotlib.pyplot as plt


def main(base_name):
    # Parameters
    Ne_present = 1_000
    # Define demography
    demography = msprime.Demography()
    demography.add_population(initial_size=Ne_present)
    # Validate demographic trajectory
    debug = demography.debug()
    # Convert to Demes graph
    graph = demography.to_demes()
    demes.dump(graph, f"{base_name}.yaml")

    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # First plot: Demes tubes
    demesdraw.tubes(graph, ax=axes[0])
    axes[0].set_title("Demography Tubes")

    # Second plot: Size history
    demesdraw.size_history(graph, invert_x=True, log_time=True, ax=axes[1])
    axes[1].set_title("Size History")

    # Save the combined figure
    plt.tight_layout()
    fig.savefig(f"{base_name}.pdf")
    print(f"Demography saved as {base_name}.yaml and {base_name}.pdf")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <base_name>")
        sys.exit(1)
    main(sys.argv[1])
