import sys
import msprime
import tskit
import demesdraw
import demes
import numpy as np
import matplotlib.pyplot as plt


def main(base_name):
    # Parameters
    Ne_present = 10_000
    Ne_anc = 15_000
    time_inv = 75
    n_founders = 10
    time_carrying = 25

    # Compute growth rate
    alpha = -np.log(n_founders / Ne_present) / (time_inv - time_carrying)

    # Define demography
    demography = msprime.Demography()
    demography.add_population(initial_size=Ne_present)
    demography.add_population_parameters_change(
        time=time_carrying, initial_size=Ne_present, growth_rate=alpha
    )
    demography.add_population_parameters_change(
        time=time_inv + 1, initial_size=Ne_anc, growth_rate=0.0
    )

    # Validate demographic trajectory
    debug = demography.debug()
    checkpoints = debug.population_size_trajectory(
        [0, time_carrying, time_inv, time_inv + 1]
    ).reshape(-1)
    assert np.allclose(
        checkpoints, [Ne_present, Ne_present, n_founders, Ne_anc], atol=1e-5
    ), "Demography validation failed"

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
