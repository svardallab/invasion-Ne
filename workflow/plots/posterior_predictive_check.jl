using Pkg
Pkg.activate(".")
using ArviZ, Plots, StatsBase, DataFrames, LaTeXStrings, NCDatasets, CSV, Distributions

# Posterior predictive check
idata = from_netcdf(snakemake.input["inference"])
post_means = mean(idata.posterior)
r2 = idata.posterior[:r2]
r2_mean = let
	# Combine draws and chains into a single dimension
	r2_flat = reshape(r2, size(r2, 1), :)  # shape becomes (19, 8000)
	vec(mean(r2_flat, dims=2))
end
# Observed data
colnames = ["bin_index", "left_bin", "right_bin", "N", "mean", "var" ]
data = CSV.read(snakemake.input["data"], DataFrame; comment="#", header=colnames)
df_bin = combine(groupby(data, :bin_index)) do sdf
	# Assumes left_bin and right_bin are constant within each bin_index group
	left = first(sdf.left_bin)
	right = first(sdf.right_bin)
	midpoint = (left + right) / 2

	(; std_mean = std(sdf.mean), midpoint)
end

# Parameters
Ne1 = parse(Float64, snakemake.wildcards["ne1"])
Ne2 = snakemake.params["ne2"]
founders = parse(Float64, snakemake.wildcards["founders"])
t0 = parse(Float64, snakemake.wildcards["t0"])
alpha = log(Ne1/founders) / t0
# Plot
default(
    fontfamily = "Computer Modern",
    size = (900, 600),
    legendfontsize = 10,
    guidefontsize = 14,
    tickfontsize = 12,
    grid = :none
)
p = plot(
	df_bin.midpoint .* 100,
	r2_mean,
	ribbon = r2_mean .-quantile.(Normal.(r2_mean, df_bin.std_mean), 0.025),
	label="95% HDI posterior predictive distribution",
    xlabel="Genetic distance (cM)",
    ylabel="Mean LD across pairs of loci",
    legend=:outerbottom,
    lw=2,
	width=3,
    fillalpha=0.2,
    framestyle=:box,
	dpi=300
)
scatter!((data.left_bin.+data.right_bin) ./ 2 .* 100, data.mean, label="Observed")
title!(L"Ne_1=%$(round(post_means.Ne1)),Ne_2=%$(round(post_means.Ne2)),Ne_f=%$(round(post_means.founders)),t_0=%$(round(post_means.t0))")
for outfile in values(snakemake.output)
	savefig(outfile)
end
