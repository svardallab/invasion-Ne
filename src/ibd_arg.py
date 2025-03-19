import sys
import tskit, tszip
import pandas as pd


def main(infile):
    # Read demography from demes
    ts = tszip.load(infile)
    min_span = 1e6
    ibds = ts.ibd_segments(min_span=1e6, store_pairs=True, store_segments=True)
    table_rows = []
    for pair in ibds:
        is_roh = (pair[0] // 2) == (pair[1] // 2)
        if is_roh:
            continue
        for ibd in ibds.get(pair):
            table_rows.append(
                {
                    "sample1_id": pair[0] // 2,
                    "hap1_idx": (pair[0] % 2) + 1,
                    "sample2_id": pair[1] // 2,
                    "hap2_idx": (pair[1] % 2) + 1,
                    "chromosome": 1,
                    "start_position": int(ibd.left),
                    "end_position": int(ibd.right),
                    "length_cm": (ibd.right - ibd.left) / 1e6,
                }
            )
    pd.DataFrame(table_rows).to_csv("/dev/stdout", sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python ibd_arg.py <input> > <output>")
        sys.exit(1)
    main(sys.argv[1])
