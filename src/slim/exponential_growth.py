import numpy as np
from slimwrap import SLiMModel
import tskit, pyslim, msprime
import tempfile, sys, os

# Define constants of the model
SAMPLE_SIZE = 200
ANCIENT_NE = 15_000
CONTIG_LENGTH = 1e8
RECOMBINATION_RATE = 1e-8
MUTATION_RATE = 1e-8
MODEL_CODE = """
initialize() {
	initializeTreeSeq();
	initializeMutationRate(0);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, L);
	initializeRecombinationRate(RHO);
	// This code assumes a POPSIZES vector exists
	defineConstant("RUNTIME", length(POPSIZES));
}
1 early() {
	sim.addSubpop("p0", POPSIZES[0]);
}
2: early() {
	p0.setSubpopulationSize(POPSIZES[sim.cycle]);
}
(RUNTIME-1) late() {
	sim.treeSeqOutput(OUTFILE);
}
"""

def ne_trajectory(Ne_modern: int, Ne_founder: int, runtime: int)-> np.ndarray:
    # We aim to fit a exponential growth
    # Ne(t) = Ne1 * exp(-\alpha * t)
    # Ne_founder = Ne1 * exp(-\alpha * runtime)
    # alpha = (log(Ne1) - log(Ne_founder)) / runtime
    # alpha <- (log(Ne1) - log(Ne_founder)) / runtime
    alpha = (np.log(Ne_modern) - np.log(Ne_founder)) / runtime
    return np.flip(Ne_modern * np.exp(-alpha*np.arange(runtime))).astype(int)

def simulation(seed: int, Ne_modern: int, Ne_founder: int, t_inv: int, outfile: str)-> None:
    # Run SLiM simulation
    model = SLiMModel(model_code=MODEL_CODE)
    with tempfile.NamedTemporaryFile(delete=False) as temp_outfile:
        params = {
            "L" : int(CONTIG_LENGTH),
            "RHO" : RECOMBINATION_RATE,
            "POPSIZES" :ne_trajectory(Ne_modern, Ne_founder, t_inv),
            "OUTFILE" : temp_outfile.name
        }
        model.run(seed=seed, constants=params)
        ts = tskit.load(temp_outfile.name)
    # Take sample of individuals
    rng = np.random.default_rng(seed=seed)
    alive_inds = pyslim.individuals_alive_at(ts, 0)
    keep_indivs = rng.choice(alive_inds, SAMPLE_SIZE, replace=False)
    keep_nodes = []
    for i in keep_indivs:
        keep_nodes.extend(ts.individual(i).nodes)
    sts = ts.simplify(keep_nodes, keep_input_roots=True)
    # Recapitation
    rts = pyslim.recapitate(
            sts, ancestral_Ne=ANCIENT_NE,
            recombination_rate=RECOMBINATION_RATE,
            random_seed=seed
    )
    # Add neutral mutations
    mts = msprime.sim_mutations(rts, rate = MUTATION_RATE)
    mts.dump(outfile)

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python script.py <seed> <recent_ne> <founders_ne> <t_inv> <outfiles>")
        sys.exit(1)
    seed = int(sys.argv[1])
    recent_ne = int(sys.argv[2])
    founders_ne = int(sys.argv[3])
    t_inv = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    outfiles = sys.argv[5:]
    seeds = rng.integers(1, 2**32, len(outfiles))
    for outfile, seed in zip(outfiles, seeds):
        print(f"Simulating chromosome with seed {seed}")
        simulation(seed, recent_ne, founders_ne, t_inv, outfile)
