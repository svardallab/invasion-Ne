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


def simulation(seed: int, recent_ne: int, outfile: str)-> None:
    # Run SLiM simulation
    model = SLiMModel(model_code=MODEL_CODE)
    with tempfile.NamedTemporaryFile(delete=False) as temp_outfile:
        params = {
            "L" : int(CONTIG_LENGTH),
            "RHO" : RECOMBINATION_RATE,
            "POPSIZES" : np.repeat(recent_ne, 100),
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
    if len(sys.argv) < 5:
        print("Usage: python script.py <seed> <recent_ne> <outfiles>")
        sys.exit(1)
    seed = int(sys.argv[1])
    recent_ne = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    outfiles = sys.argv[3:]
    seeds = rng.integers(1, 2**32, len(outfiles))
    for outfile, seed in zip(outfiles, seeds):
        print(f"Simulating chromosome with seed {seed}")
        simulation(seed, recent_ne, outfile)
