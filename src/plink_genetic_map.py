import tskit, tszip
import json, sys


# This file creates a genetic map according to the PLINK .map specification
def main(infile, chromosome):
    # Read tszip file
    ts = tszip.decompress(infile)
    # Get metadata
    provenance = ts.provenance(0)
    record_json = provenance.asdict()["record"]
    metadata = json.loads(record_json).get("parameters")
    assert metadata["command"] == "sim_ancestry"
    recombination_rate = metadata["recombination_rate"]
    sequence_length = metadata["sequence_length"]
    # Print to stdout the genetic map
    print(f"{chromosome}\trs1\t0\t0")
    cm_length = sequence_length * recombination_rate * 100
    print(f"{chromosome}\trs2\t{cm_length}\t{int(sequence_length)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.ts.zip> <chromosome>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
