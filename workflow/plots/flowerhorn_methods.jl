using Pkg
Pkg.activate(".")
using ArviZ, Plots, StatsBase, DataFrames, LaTeXStrings, NCDatasets, CSV

# Posterior prediction
times = 0:125
function ne(Ne1, Ne2, t0, alpha)
	ifelse.(times .<= t0, Ne1 * exp.(-alpha .* times), Ne2)
end
function trajectory(idata)
    Ne1 = vec(idata.posterior[:Ne1])
    Ne2 = vec(idata.posterior[:Ne2])
    t0 = vec(idata.posterior[:t0])
    alpha = vec(idata.posterior[:alpha])
    n_draws = length(Ne1)
    n_times = length(times)
    result = Matrix{Float64}(undef, n_draws, n_times)
    for i in 1:n_draws
        result[i, :] .= ne(Ne1[i], Ne2[i], t0[i], alpha[i])
    end
    return result
end
function summarize_trajectories(idata)
    traj_matrix = trajectory(idata)  # as defined earlier
    # Compute 95% HDIs for each time point (column-wise)
    hdis = hcat([
        collect(hdi(col; prob=0.95))
        for col in eachcol(traj_matrix)
    ]...)'  # Transpose: now each row corresponds to a time point
    # Compute the mean for each column (time point)
    mean_vals = vec(mean(traj_matrix, dims=1))  # Make it a 1D vector
    # Assemble DataFrame
    df = DataFrame(
        time = times,
        lower = hdis[:, 1],
        upper = hdis[:, 2],
        mean = mean_vals
    )
    return df
end

# Input files
bayesbagg = snakemake.input["bayesbagg"]
outfiles = snakemake.output
# Parameters
Ne1 = parse(Float64, snakemake.wildcards["ne1"])
n_samples = parse(Int, snakemake.wildcards["sample_size"])
Ne2_1 = snakemake.params["ne2_1"]
Ne2_2 = snakemake.params["ne2_2"]
founders = parse(Float64, snakemake.wildcards["founders"])
t0 = parse(Float64, snakemake.wildcards["t0"])
alpha = log(Ne1/founders) / t0
# Reading data
# Compute expected trajectories across resampled datasets
trajs_bayesbag = reduce(hcat, [
        summarize_trajectories(from_netcdf(infile)).mean
        for infile in bayesbagg ]
)

# Plot the data
default(
    fontfamily = "Computer Modern",
    size = (900, 600),
    legendfontsize = 10,
    guidefontsize = 14,
    tickfontsize = 12,
    grid = :none
)

lower = [quantile(row, 0.025) for row in eachrow(trajs_bayesbag)]
upper = [quantile(row, 0.975) for row in eachrow(trajs_bayesbag)]
expected = mean(trajs_bayesbag, dims=2)

p = plot(
    times,
    expected,
    ribbon=(expected .- lower, upper .- expected),
    label="95% BayesBagg (own method)",
    xlabel="Time (generations ago)",
    ylabel="Effective population size (Ne)",
    legend=:outerbottom,
    lw=2,
	width=3,
    fillalpha=0.2,
    framestyle=:box,
	dpi=300
)

plot!(
    times,
    ne(Ne1, Ne2_1, t0, alpha),
    label="Ground truth (hybrid ancestry 1)",
    lw=3,
	color=:black,
	width=3,
    linestyle=:dash
)

plot!(
    times,
    ne(Ne1, Ne2_2, t0, alpha),
    label="Ground truth (hybrid ancestry 2)",
    lw=3,
	color=:black,
	width=3,
    linestyle=:dash
)

ylims!(0, 30_000)
xlims!(0, 100)
title!(L"Ne_1=%$Ne1,Ne^1_2=%$Ne2_1,Ne^2_2=%$Ne2_2,Ne_f=%$founders,t_0=%$t0 (n=%$n_samples)")
for outfile in values(outfiles)
	savefig(outfile)
end
