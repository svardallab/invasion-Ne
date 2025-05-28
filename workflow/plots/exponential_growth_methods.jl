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
hapne_ld_file = snakemake.input["hapne_ld_file"]
hapne_ibd_file = snakemake.input["hapne_ibd_file"]
ibdne_file = snakemake.input["ibdne_file"]
gone2_file = snakemake.input["gone2_file"]
bayes_ld_file = snakemake.input["bayes_ld_file"]
outfiles = snakemake.output
# Parameters
Ne1 = parse(Float64, snakemake.wildcards["ne1"])
Ne2 = snakemake.params["ne2"]
founders = parse(Float64, snakemake.wildcards["founders"])
t0 = parse(Float64, snakemake.wildcards["t0"])
alpha = log(Ne1/founders) / t0
# Reading data
df_hapneld = CSV.read(hapne_ld_file, DataFrame)
df_hapneibd = CSV.read(hapne_ibd_file, DataFrame)
df_ibdne = CSV.read(ibdne_file, DataFrame)
df_gone = CSV.read(gone2_file, DataFrame)
# Load from NC dataset
idata = from_netcdf(bayes_ld_file)
df = summarize_trajectories(idata)
# Plot the data
default(
    fontfamily = "Computer Modern",
    size = (900, 600),
    legendfontsize = 10,
    guidefontsize = 14,
    tickfontsize = 12,
    grid = :none
)

p = plot(
    df.time,
    df.mean,
    ribbon=(df.mean .- df.lower, df.upper .- df.mean),
    label="95% HPI (own method)",
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
    df_hapneld.TIME,
    df_hapneld[!, "Q0.5"],
    ribbon=(df_hapneld[!, "Q0.5"] .- df_hapneld[!, "Q0.025"], df_hapneld[!, "Q0.975"] .- df_hapneld[!, "Q0.5"]),
    label="95% bootstrapping CI (HapNe-LD)",
    lw=2,
    fillalpha=0.2
)

plot!(
    df_gone.Generation,
    df_gone.Ne_diploids,
    label="GONE2",
    lw=2,
)

plot!(
    df_hapneibd.TIME,
    df_hapneibd[!, "Q0.5"],
    ribbon=(df_hapneibd[!, "Q0.5"] .- df_hapneibd[!, "Q0.025"], df_hapneibd[!, "Q0.975"] .- df_hapneibd[!, "Q0.5"]),
    label="95% bootstrapping CI (HapNe-IBD)",
    lw=2,
    fillalpha=0.2
)

plot!(
    df_ibdne.GEN,
    df_ibdne.NE,
    ribbon=(df_ibdne.NE .- df_ibdne[!, "LWR-95%CI"], df_ibdne[!, "UPR-95%CI"] .- df_ibdne.NE),
    label="95% bootstrapping CI (IBDNe)",
    lw=2,
    fillalpha=0.2
)

plot!(
    times,
    ne(Ne1, Ne2, t0, alpha),
    label="Ground truth",
    lw=3,
	color=:black,
	width=3,
    linestyle=:dash
)

ylims!(0, 30_000)
xlims!(0, 100)
title!(L"Ne_1=%$Ne1,Ne_2=%$Ne2,Ne_f=%$founders,t_0=%$t0")
for outfile in values(outfiles)
	savefig(outfile)
end
