import tskit, tszip
import json, sys

def main(infiles, chromosome_list):
    chromosomes = chromosome_list.split(',')
    print("CHR\tFROM_BP\tTO_BP\tNAME\tLENGTH")
    for infile, chrom in zip(infiles, chromosomes):
        ts = tszip.decompress(infile)
        # Get metadata
        provenance = ts.provenance(0)
        record_json = provenance.asdict()["record"]
        metadata = json.loads(record_json).get("parameters")
        assert metadata["command"] == "sim_ancestry"
        recombination_rate = metadata["recombination_rate"]
        sequence_length = metadata["sequence_length"]
        cm_length = sequence_length * recombination_rate * 100
        print(f"{chrom}\t0\t{sequence_length}\t{chrom}\t{cm_length}")


if __name__ == "__main__":
    main(sys.argv[1:-1], sys.argv[-1])
