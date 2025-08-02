using CSV, DataFrames, Statistics
using ArviZ, NCDatasets
using Distributions, Optim, Random
import NaNMath

times = 0:100
basedir = "../../../steps/"

function trajectory(Ne_c, Ne_f, t_inv)
    Ne_anc = 15_000
    alpha = (log(Ne_c) - log(Ne_f)) / t_inv
    ifelse.(times .<= t_inv, Ne_c * exp.(-alpha*times), Ne_anc)
end

function model(t, p)
    Ne_c, Ne_f, t_inv, Ne_anc = p
    alpha = (NaNMath.log(Ne_c) - NaNMath.log(Ne_f)) / t_inv
    return @. ifelse(t <= t_inv, Ne_c * NaNMath.exp(-alpha * t), Ne_anc)
end

function fit(x::AbstractVector, y::AbstractVector)
    Random.seed!(123) # for reproducibility
    p0 = [
        rand(Truncated(Normal(10_000, 5000), 0.0, Inf)), # Ne_c
        rand(Truncated(Normal(1000, 500), 0.0, Inf)),   # Ne_f
        rand(Truncated(Normal(50, 25), 0.0, Inf)),      # t_inv (time must be positive)
        rand(Truncated(Normal(15_000, 5000), 0.0, Inf)),# Ne_anc
    ]
    objective(p) = sum(abs2.(model(x, p) .- y))
    lower = [0.0, 0.0, 0.0, 0.0]
    upper = [Inf, Inf, Inf, Inf]
    inner_optimizer = GradientDescent()
    results = optimize(objective, lower, upper, p0, Fminbox(inner_optimizer))
    fitted_params = Optim.minimizer(results)
    return fitted_params
end

function estimate_own(Ne_c, Ne_f, t_inv, seed)
    method = "Own"
    infile = basedir * "inference/exponential_piecewise_model/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed).nc"
    idata = ArviZ.from_netcdf(infile)
    post = idata.posterior
    DataFrame(
        Ground_truth = [Ne_c, Ne_f, t_inv, 15_000],
        Estimates = [mean(post.Ne1), mean(post.founders), mean(post.t0), mean(post.Ne2)],
        Method = method,
        Param= ["Ne_c", "Ne_f", "t_inv", "Ne2"]
    )
end

function estimate_gone2(Ne_c, Ne_f, t_inv, seed)
    method = "GONE2"
    infile = basedir * "gone2/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed)_GONE2_Ne"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(Generation = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :Generation)
    params = fit(df_merged.Generation, df_merged.Ne_diploids)
    DataFrame(
        Ground_truth = [Ne_c, Ne_f, t_inv, 15_000],
        Estimates = params,
        Method = method,
        Param= ["Ne_c", "Ne_f", "t_inv", "Ne2"]
    )
end

function estimate_hapne_ibd(Ne_c, Ne_f, t_inv, seed)
    method = "HapNe_IBD"
    infile = basedir * "hapne/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed)_ibd_hapne_estimate.csv"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(TIME = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :TIME)
    params = fit(df_merged.TIME, df_merged[!, Symbol("Q0.5")]./ 2)
    DataFrame(
        Ground_truth = [Ne_c, Ne_f, t_inv, 15_000],
        Estimates = params,
        Method = method,
        Param= ["Ne_c", "Ne_f", "t_inv", "Ne2"]
    )
end

function estimate_hapne_ld(Ne_c, Ne_f, t_inv, seed)
    method = "HapNe_LD"
    infile = basedir * "hapne/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed)_ld_hapne_estimate.csv"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(TIME = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :TIME)
    params = fit(df_merged.TIME, df_merged[!, Symbol("Q0.5")]./ 2)
    DataFrame(
        Ground_truth = [Ne_c, Ne_f, t_inv, 15_000],
        Estimates = params,
        Method = method,
        Param= ["Ne_c", "Ne_f", "t_inv", "Ne2"]
    )
end

function estimate_ibdne(Ne_c, Ne_f, t_inv, seed)
    method = "IBDNe"
    infile = basedir * "ibdne/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed).ne"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(GEN = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :GEN)
    params = fit(df_merged.GEN, df_merged.NE)
    DataFrame(
        Ground_truth = [Ne_c, Ne_f, t_inv, 15_000],
        Estimates = params,
        Method = method,
        Param= ["Ne_c", "Ne_f", "t_inv", "Ne2"]
    )
end

function gather_estimates()
    Ne_c = 10_000
    Ne_f_vals = [10, 100, 1000]
    t_inv_vals = [25, 50, 75]
    seeds = 100:125
    methods = [estimate_gone2, estimate_hapne_ibd, estimate_hapne_ld, estimate_ibdne, estimate_own]

    df = DataFrame()

    for f in methods, Ne_f in Ne_f_vals, t_inv in t_inv_vals, seed in seeds
        try
            append!(df, f(Ne_c, Ne_f, t_inv, seed))
        catch err
            @warn "Error in method=$(nameof(f)) with args: Ne_c=$Ne_c, Ne_f=$Ne_f, t_inv=$t_inv, seed=$seed"
            @warn err
        end
    end

    return df
end

df = gather_estimates()
CSV.write("../methods_parametric_point_estimates.csv", df)
