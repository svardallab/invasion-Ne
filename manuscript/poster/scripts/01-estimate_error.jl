using CSV, DataFrames, Statistics
using ArviZ, NCDatasets

times = 0:100
basedir = "../../../steps/"

function trajectory(Ne_c, Ne_f, t_inv)
    Ne_anc = 15_000
    alpha = (log(Ne_c) - log(Ne_f)) / t_inv
    ifelse.(times .<= t_inv, Ne_c * exp.(-alpha*times), Ne_anc)
end

function error_own(Ne_c, Ne_f, t_inv, seed)
    method = "Own"
    infile = basedir * "inference/exponential_piecewise_model/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed).nc"
    idata = ArviZ.from_netcdf(infile)
    df2 = DataFrame(Generation = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))

    function prediction(Ne_c, Ne_f, t_inv, Ne_anc, gen)
        alpha = (log(Ne_c) - log(Ne_f)) / t_inv
        ifelse.(gen .<= t_inv, Ne_c * exp.(-alpha * gen), Ne_anc)
    end

    post = idata.posterior
    Ne_c_draw = vec(post.Ne1)
    Ne_f_draw = vec(post.founders)
    t_inv_draw = vec(post.t0)
    Ne_anc_draw = vec(post.Ne2)

    preds = reduce(hcat, [
        prediction(Ne_c_draw[i], Ne_f_draw[i], t_inv_draw[i], Ne_anc_draw[i], df2.Generation)
        for i in 1:length(Ne_c_draw)
    ])

    estimate = vec(mean(preds, dims=2))
    errs = abs.(df2.Ground_truth .- estimate) ./ df2.Ground_truth

    DataFrame(
        Generation = df2.Generation,
        Method = method,
        Seed = seed,
        Ne_c = Ne_c, Ne_f = Ne_f, t_inv = t_inv,
        RelativeError = errs
    )
end

function error_own_informed(Ne_c, Ne_f, t_inv, seed)
    method = "Own_informed"
    infile = basedir * "inference/exponential_piecewise_model_informed/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed).nc"
    idata = ArviZ.from_netcdf(infile)
    df2 = DataFrame(Generation = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))

    function prediction(Ne_c, Ne_f, t_inv, Ne_anc, gen)
        alpha = (log(Ne_c) - log(Ne_f)) / t_inv
        ifelse.(gen .<= t_inv, Ne_c * exp.(-alpha * gen), Ne_anc)
    end

    post = idata.posterior
    Ne_c_draw = vec(post.Ne1)
    Ne_f_draw = vec(post.founders)
    t_inv_draw = vec(post.t0)
    Ne_anc_draw = vec(post.Ne2)

    preds = reduce(hcat, [
        prediction(Ne_c_draw[i], Ne_f_draw[i], t_inv_draw[i], Ne_anc_draw[i], df2.Generation)
        for i in 1:length(Ne_c_draw)
    ])

    estimate = vec(mean(preds, dims=2))
    errs = abs.(df2.Ground_truth .- estimate) ./ df2.Ground_truth

    DataFrame(
        Generation = df2.Generation,
        Method = method,
        Seed = seed,
        Ne_c = Ne_c, Ne_f = Ne_f, t_inv = t_inv,
        RelativeError = errs
    )
end

function error_gone2(Ne_c, Ne_f, t_inv, seed)
    method = "GONE2"
    infile = basedir * "gone2/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed)_GONE2_Ne"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(Generation = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :Generation)
    errs = abs.(df_merged.Ground_truth .- df_merged.Ne_diploids) ./ df_merged.Ground_truth
    DataFrame(
        Generation = df_merged.Generation,
        Method = method,
        Seed = seed,
        Ne_c = Ne_c, Ne_f = Ne_f, t_inv = t_inv,
        RelativeError = errs
    )
end

function error_hapne_ibd(Ne_c, Ne_f, t_inv, seed)
    method = "HapNe_IBD"
    infile = basedir * "hapne/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed)_ibd_hapne_estimate.csv"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(TIME = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :TIME)
    errs = abs.(df_merged.Ground_truth .- df_merged[!, Symbol("Q0.5")]./ 2) ./ df_merged.Ground_truth
    DataFrame(
        Generation = df_merged.TIME,
        Method = method,
        Seed = seed,
        Ne_c = Ne_c, Ne_f = Ne_f, t_inv = t_inv,
        RelativeError = errs
    )
end

function error_hapne_ld(Ne_c, Ne_f, t_inv, seed)
    method = "HapNe_LD"
    infile = basedir * "hapne/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed)_ld_hapne_estimate.csv"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(TIME = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :TIME)
    errs = abs.(df_merged.Ground_truth .- df_merged[!, Symbol("Q0.5")] ./ 2) ./ df_merged.Ground_truth
    DataFrame(
        Generation = df_merged.TIME,
        Method = method,
        Seed = seed,
        Ne_c = Ne_c, Ne_f = Ne_f, t_inv = t_inv,
        RelativeError = errs
    )
end

function error_ibdne(Ne_c, Ne_f, t_inv, seed)
    method = "IBDNe"
    infile = basedir * "ibdne/exponential_growth/ne1_$(Ne_c)_ne2_$(Ne_f)_t$(t_inv)/s$(seed).ne"
    df = CSV.read(infile, DataFrame)
    df2 = DataFrame(GEN = times, Ground_truth = trajectory(Ne_c, Ne_f, t_inv))
    df_merged = innerjoin(df, df2, on = :GEN)
    errs = abs.(df_merged.Ground_truth .- df_merged.NE) ./ df_merged.Ground_truth
    DataFrame(
        Generation = df_merged.GEN,
        Method = method,
        Seed = seed,
        Ne_c = Ne_c, Ne_f = Ne_f, t_inv = t_inv,
        RelativeError = errs
    )
end

function gather_errors()
    Ne_c = 10_000
    Ne_f_vals = [10, 100, 1000]
    t_inv_vals = [25, 50, 75]
    seeds = [100]
    methods = [error_gone2, error_hapne_ibd, error_hapne_ld, error_ibdne, error_own, error_own_informed]

    df = DataFrame()

    for f in methods, Ne_f in Ne_f_vals, t_inv in t_inv_vals, seed in seeds
            append!(df, f(Ne_c, Ne_f, t_inv, seed))
    end
    return df
end

df = gather_errors()
CSV.write("../methods_mean_absolute_relative_error.csv", df)
