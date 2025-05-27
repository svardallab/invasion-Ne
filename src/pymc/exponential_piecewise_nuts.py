import pandas as pd
import numpy as np
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import sys


# From 0 to t0
def S_ut_piece1(alpha, Ne1, t, u):
    t = pt.as_tensor_variable(t)[:, None]  # Shape (n_quad, 1)
    u = pt.as_tensor_variable(u)[None, :]  # Shape (1, n_points)
    # If alpha is not close to zero
    inner1 = (1 - pt.exp(alpha * t)) / (2 * Ne1 * alpha)
    exponent1 = alpha * t - 2 * t * u + inner1
    res1 = pt.exp(exponent1) / (2 * Ne1)  # Shape (n_quad, n_points)
    # If alpha is close to zero we use Taylor series
    numerator = 4 * Ne1 + alpha * t * (4 * Ne1 - t)
    exponent2 = -t * (4 * Ne1 * u + 1) / (2 * Ne1)
    res2 = numerator * pt.exp(exponent2) / (8 * Ne1**2)
    epsilon = 1e-5
    return pt.switch(pt.abs(alpha) < epsilon, res2, res1)


# From t0 to infinity
def S_ut_piece2(alpha, Ne1, Ne2, t0, t, u):
    t = pt.as_tensor_variable(t)[:, None]  # Shape (n_quad, 1)
    u = pt.as_tensor_variable(u)[None, :]  # Shape (1, n_points)
    # If alpha is not close to zero
    inner1 = (Ne1 * alpha * (t0 - t) + Ne2 * (1 - pt.exp(alpha * t0))) / (
        2 * Ne1 * Ne2 * alpha
    )
    exponent1 = -2 * t * u + inner1
    res1 = pt.exp(exponent1) / (2 * Ne2)  # Shape (n_quad, n_points)
    # If alpha is close to zero we use Taylor series
    inner2 = 4 * Ne1 - alpha * t0**2
    exponent2 = (-4 * Ne1 * Ne2 * t * u + Ne1 * (t0 - t) - Ne2 * t0) / (2 * Ne1 * Ne2)
    res2 = inner2 * pt.exp(exponent2) / (8 * Ne1 * Ne2)
    epsilon = 1e-5
    return pt.switch(pt.abs(alpha) < epsilon, res2, res1)


# Numerical integration using Gaussian quadrature rules
def expected_r2(u_col, Ne1, Ne2, alpha, t0, legendre_x, legendre_w):
    u_col = pt.as_tensor_variable(u_col)
    # First integral: [0, t0]
    times1 = (t0 - 0) / 2 * legendre_x + (t0 + 0) / 2
    f_t_piece1 = S_ut_piece1(alpha, Ne1, times1, u_col)  # (n_quad, n_points)
    integral_piece1 = pt.sum(
        f_t_piece1 * legendre_w[:, None] * (t0 - 0) / 2, axis=0
    )  # (n_points,)

    # Second integral: [t0, ∞)
    trans_legendre_x = 0.5 * legendre_x + 0.5
    trans_legendre_w = 0.5 * legendre_w
    times2 = t0 + trans_legendre_x / (1 - trans_legendre_x)
    f_t_piece2 = S_ut_piece2(alpha, Ne1, Ne2, t0, times2, u_col)
    integral_piece2 = pt.sum(
        f_t_piece2 * (trans_legendre_w[:, None] / (1 - trans_legendre_x)[:, None] ** 2),
        axis=0,
    )  # (n_points,)

    return integral_piece1 + integral_piece2  # shape (n_points,)


def correct_r2(mu, sample_size):
    S = sample_size * 2  # diploid assumption
    beta = 1 / (S - 1) ** 2
    alpha = ((S**2 - S + 2) ** 2) / ((S**2 - 3 * S + 2) ** 2)
    return (alpha - beta) * mu + 4 * beta


# Helper function to re-scale Legendre Gaussian quadrature rules
def gauss(a, b, n=10):
    x, w = np.polynomial.legendre.leggauss(n)
    w = (b - a) / 2 * w
    x = (b - a) / 2 * x + (a + b) / 2
    return x, w


def main(
    ld_file: str,
    ne_anc_file: str,
    ne1_prior_mean: float,
    ne1_prior_sd: float,
    t0_prior_mean: float,
    t0_prior_sd: float,
    alpha_logfold_prior_sd: float,
    sample_size: int,
    seed: int,
    outfile: str,
) -> None:
    print(f"Running on PyMC v{pm.__version__}")
    print(f"Processing file: {ld_file}")

    # Read data
    df = pd.read_csv(
        ld_file,
        delimiter="\t",
        comment="#",
        names=["bin_index", "left_bin", "right_bin", "N", "mean", "var"],
    )
    # Extract data for model
    Nrows = len(df)
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
    print("Processing file:", ne_anc_file)
    ne_df = pd.read_csv(ne_anc_file)
    # Calcula mean and std of the Ne values across all chromosomes
    ne2_prior_mean = ne_df["Ne"].mean()
    ne2_prior_sd = ne_df["Ne"].std()
    print("Ne2 prior mean:", ne2_prior_mean)
    print("Ne2 prior std:", ne2_prior_sd)
    with pm.Model() as model:
        # Prior for Ne1
        # We use a truncated normal with variance to make things easier for NUTS
        # but restrict Ne(t) to be positive
        Ne1_raw = pm.TruncatedNormal(
            "Ne1_raw", mu=0, sigma=1, lower=-ne1_prior_mean / ne1_prior_sd
        )
        Ne1 = pm.Deterministic("Ne1", ne1_prior_mean + ne1_prior_sd * Ne1_raw)
        # Prior for Ne2
        Ne2_raw = pm.TruncatedNormal(
            "Ne2_raw", mu=0, sigma=1, lower=-ne2_prior_mean / ne2_prior_sd
        )
        Ne2 = pm.Deterministic("Ne2", ne2_prior_mean + ne2_prior_sd * Ne2_raw)
        # Prior for t0
        t0_raw = pm.TruncatedNormal(
            "t0_raw", mu=0, sigma=1, lower=-t0_prior_mean / t0_prior_sd
        )
        t0 = pm.Deterministic("t0", t0_prior_mean + t0_prior_sd * t0_raw)
        # Prior for alpha
        # We have to restrict combination that lead to a Ne(t0) < 1
        alpha_raw = pm.TruncatedNormal(
            "alpha_raw", mu=0, sigma=1, upper=pt.log(Ne1) / alpha_logfold_prior_sd
        )
        alpha = pm.Deterministic("alpha", alpha_raw * alpha_logfold_prior_sd / t0)
        pm.Deterministic("founders", Ne1 * pt.exp(-alpha * t0))

        # Numerical integration across both time (0->Inf) and bin
        # Per timepoint points and weights
        legendre_x, legendre_w = np.polynomial.legendre.leggauss(100)
        # Per bin quadrature points and weights
        u_points = np.array([gauss(a, b, 10)[0] for (a, b) in zip(u_i, u_j)])
        u_weights = np.array([gauss(a, b, 10)[1] / (b - a) for (a, b) in zip(u_i, u_j)])

        # Compute R2 for different lengths in one pass
        r2_flat = expected_r2(
            u_points.flatten(), Ne1, Ne2, alpha, t0, legendre_x, legendre_w
        )
        # Reshape it into a matrix and compute integral
        r2_matrix = r2_flat.reshape(u_points.shape)
        r2_per_bin = pt.sum(r2_matrix * u_weights, axis=1)  # shape (n_bins,)
        r2_corrected = pm.Deterministic("r2", correct_r2(r2_per_bin, sample_size))

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
            chains=4, tune=60_000, draws=40_000, target_accept=0.90, random_seed=seed
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
    if len(sys.argv) != 11:
        print(
            "Usage: python exponential_piecewise_nuts.py <ld_file> <ne_anc_file> <ne1_prior_mean> <ne1_prior_sd> <t0_prior_mean> <t0_prior_sd> <alpha_logfold_prior_sd> <sample_size> <seed> <output_file>"
        )
        sys.exit(1)
    ld_file = sys.argv[1]
    ne_anc_file = sys.argv[2]
    ne1_prior_mean = float(sys.argv[3])
    ne1_prior_sd = float(sys.argv[4])
    t0_prior_mean = float(sys.argv[5])
    t0_prior_sd = float(sys.argv[6])
    alpha_logfold_prior_sd = float(sys.argv[7])
    sample_size = int(sys.argv[8])
    seed = int(sys.argv[9])
    outfile = sys.argv[10]
    main(
        ld_file,
        ne_anc_file,
        ne1_prior_mean,
        ne1_prior_sd,
        t0_prior_mean,
        t0_prior_sd,
        alpha_logfold_prior_sd,
        sample_size,
        seed,
        outfile,
    )
