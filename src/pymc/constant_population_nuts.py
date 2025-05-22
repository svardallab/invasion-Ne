import pandas as pd
import numpy as np
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import sys


# Expected r^2 for constant Ne model
def expected_r2(u_i, u_j, Ne):
    u_i = pt.as_tensor_variable(u_i)
    u_j = pt.as_tensor_variable(u_j)
    return (-pt.log(4 * Ne * u_i + 1) + pt.log(4 * Ne * u_j + 1)) / (
        4 * Ne * (u_j - u_i)
    )


# Correction for finite sample size
def correct_r2(mu, sample_size):
    S = sample_size * 2
    beta = 1 / (S - 1) ** 2
    alpha = (S**2 - S + 2) ** 2 / (S**2 - 3 * S + 2) ** 2
    return (alpha - beta) * mu + 4 * beta


def main(
    infile: str, prior_mean: float, prior_sd: float, sample_size: int, outfile: str
) -> None:
    print(f"Running on PyMC v{pm.__version__}")
    print(f"Processing file: {infile}")

    # Read data
    df = pd.read_csv(
        infile,
        delimiter="\t",
        comment="#",
        names=["bin_index", "left_bin", "right_bin", "N", "mean", "var"],
    )
    # Extract data for model
    Nrows = len(df)
    bins = df["bin_index"].values
    # Get unique bin parameters
    df_bins = df.drop_duplicates("bin_index")[["bin_index", "left_bin", "right_bin"]]
    df_bins = df_bins.sort_values("bin_index")
    u_i = df_bins["left_bin"].values
    u_j = df_bins["right_bin"].values
    Nbins = len(df_bins)
    Nchrom = Nrows // Nbins
    # Calculate expected_sigma2 per bin
    sigma2_per_bin = df.groupby("bin_index")["var"].mean().sort_index().values
    # Get bin indices for each row (make sure they match the ordered bins)
    bin_indices = np.array(
        [np.where(df_bins["bin_index"].values == b)[0][0] for b in df["bin_index"]]
    )

    with pm.Model() as model:
        # Prior for Ne
        Ne_raw = pm.TruncatedNormal(
            "Ne_raw", mu=0, sigma=1, lower=-prior_mean / prior_sd
        )
        Ne = pm.Deterministic("Ne", prior_mean + prior_sd * Ne_raw)

        # Calculate expected values per unique bin:
        r2_corrected = pm.Deterministic(
            "LD", correct_r2(expected_r2(u_i, u_j, Ne), sample_size)
        )

        # Map corrected r² values to each data row based on bin index
        r2_mapped = r2_corrected[bin_indices]
        # Map sigma² values to each data row
        sigma2_mapped = pt.as_tensor_variable(sigma2_per_bin[bin_indices])
        # Compute composite log likelihood
        log_lik = -df["N"].values * (
            0.5 * pt.log(sigma2_mapped)
            + (df["var"].values + (df["mean"].values - r2_mapped) ** 2)
            / (2 * sigma2_mapped)
        )
        # Reshape log_lik to compute the pointwise log likelihood per chromosome
        pointwise_loglik = pm.Deterministic(
            "log_likelihood", pt.sum(log_lik.reshape((-1, Nchrom)), axis=0)
        )
        pm.Potential("likelihood", pt.sum(pointwise_loglik))

        # Sample from the posterior
        idata = pm.sample(
            chains=4,
            tune=10_000,
            draws=2000,
        )
        # Add log_likelihood to its own group
        idata.add_groups(log_likelihood=idata.posterior.log_likelihood)
        # and remove it from the posterior group
        idata.posterior = idata.posterior.drop_vars("log_likelihood")

    # Print summary statistics focusing on Ne
    summary = az.summary(idata)
    print(summary)
    loo = az.loo(idata)
    print(loo)
    print("Saving data to NetCDF file...")
    idata.to_netcdf(outfile)
    return idata


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(
            "Usage: python constant_population_nuts.py <input_file> <prior_mean> <prior_sd> <sample_size> <output_file>"
        )
        sys.exit(1)
    infile = sys.argv[1]
    prior_mean = float(sys.argv[2])
    prior_sd = float(sys.argv[3])
    sample_size = int(sys.argv[4])
    outfile = sys.argv[5]
    main(infile, prior_mean, prior_sd, sample_size, outfile)
