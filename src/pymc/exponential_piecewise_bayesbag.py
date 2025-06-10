import pandas as pd
import numpy as np
VARS = ["Ne1", "Ne2", "t0", "founders", "alpha"]
def summary(infile):
    data = pd.read_csv(infile)
    data= data.set_index(data.columns[0])
    return data.loc[VARS]["mean"].to_numpy()
def main(infiles):
    summaries = np.array([summary(infile) for infile in infiles])
    df = pd.DataFrame({
        "vars" : VARS,
        "mean" : np.mean(summaries, axis=0),
        "median" : np.median(summaries, axis=0),
        "lower95%" : np.quantile(summaries, q=0.025, axis=0),
        "upper95%" : np.quantile(summaries, q=0.975, axis=0),
    })
    print(df)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
