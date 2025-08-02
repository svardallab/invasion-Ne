using CSV, DataFrames, Statistics
using ArviZ, NCDatasets
basedir = "../../../steps/"

function estimate_own(Ne_c, Ne_f, t_inv, seed)
    method = "Own"
    infile = basedir * "inference/exponential_piecewise_model/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed).nc"
    idata = ArviZ.from_netcdf(infile)
    post = idata.posterior
    lower_hdi(x)= hdi(x, prob=0.83).lower
    upper_hdi(x)= hdi(x, prob=0.83).upper
    DataFrame(
        Ground_truth = [Ne_c, Ne_f, t_inv, 15_000],
        Mean = [mean(post.Ne1), mean(post.founders), mean(post.t0), mean(post.Ne2)],
        Lower = [lower_hdi(post.Ne1), lower_hdi(post.founders), lower_hdi(post.t0), lower_hdi(post.Ne2)],
        Upper = [upper_hdi(post.Ne1), upper_hdi(post.founders), upper_hdi(post.t0), upper_hdi(post.Ne2)],
        Method = method,
        Param= ["Ne_c", "Ne_f", "t_inv", "Ne2"]
    )
end

function gather_estimates()
    Ne_c = 10_000
    Ne_f_vals = [10, 100, 1000]
    t_inv_vals = [25, 50, 75]
    seeds = 100:125

    df = DataFrame()

    for Ne_f in Ne_f_vals, t_inv in t_inv_vals, seed in seeds
        append!(df, estimate_own(Ne_c, Ne_f, t_inv, seed))
    end

    return df
end

df = gather_estimates()
CSV.write("../inference_coverage.csv", df)
