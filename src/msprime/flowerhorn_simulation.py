import numpy as np
import tskit, msprime
import sys
from concurrent.futures import ProcessPoolExecutor
# Define constants of the model
ANCIENT_NE1 = 15_000
ANCIENT_NE2 = 25_000
ANCIENT_ADMIX = 100
CONTIG_LENGTH = 1e8
SPLIT_TIME = 30_000
RECOMBINATION_RATE = 1e-8
MUTATION_RATE = 1e-8

def simulation(seed: int, Ne_modern: int, Ne_founder: int, t_inv: int, n_samples: int, outfile: str)-> None:
    # Define demography
    alpha = (np.log(Ne_modern) - np.log(Ne_founder)) / t_inv
    demography = msprime.Demography()
    demography.add_population(name="p0", initial_size=Ne_modern, growth_rate=alpha)
    demography.add_population(name="p1", initial_size = ANCIENT_NE1)
    demography.add_population(name="p2", initial_size = ANCIENT_NE2)
    demography.add_population(name="p3", initial_size=ANCIENT_NE1+ANCIENT_NE2)
    demography.add_admixture(time = t_inv, derived = "p0", ancestral=["p1", "p2"], proportions=[0.4, 0.6])
    demography.add_population_split(ancestral="p3", derived=["p1", "p2"], time=SPLIT_TIME)
    # Coalescent simulation
    ts = msprime.sim_ancestry(
        samples={"p0" : n_samples},
        recombination_rate=RECOMBINATION_RATE,
        sequence_length=CONTIG_LENGTH,
        random_seed=seed,
        demography=demography
    )
    mts = msprime.sim_mutations(ts, rate = MUTATION_RATE, random_seed=seed)
    mts.dump(outfile)

if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("Usage: python script.py <seed> <recent_ne> <founders_ne> <t_inv> <n_samples> <outfiles>")
        sys.exit(1)
    seed = int(sys.argv[1])
    recent_ne = int(sys.argv[2])
    founders_ne = int(sys.argv[3])
    t_inv = int(sys.argv[4])
    n_samples = int(sys.argv[5])
    rng = np.random.default_rng(seed)
    outfiles = sys.argv[6:]
    seeds = rng.integers(1, 2**32, len(outfiles))
    def worker(seed, outfile, recent_ne, founders_ne, t_inv, n_samples):
        print(f"Simulating chromosome with seed {seed}")
        simulation(seed, recent_ne, founders_ne, t_inv, n_samples, outfile)
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                worker, seed, outfile, recent_ne, founders_ne, t_inv, n_samples
            )
            for seed, outfile in zip(seeds, outfiles)
        ]
        for future in futures:
            future.result()
