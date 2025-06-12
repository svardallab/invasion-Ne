rule experiment:
    input:
        expand(
            "manuscript/figures/flowerhorn_methods/ne1_{ne1}_founders_{founders}_t{t_inv}_s{seed}_n{sample_size}.pdf",
            seed=range(100, 105),
            ne1=[10000],
            founders=[10, 100, 1000],
            sample_size = [10],
            t_inv=[5, 10, 15, 20],
        ),

rule plot_method_flowerhorn:
    input:
        bayesbagg=expand(
            "steps/inference/exponential_piecewise_model/flowerhorn/ne1_{{ne1}}_ne2_{{founders}}_t{{t0}}_n{{sample_size}}/s{{seed}}_b{boot}_advi.nc",
            boot=range(50)
        )
    localrule: True
    output:
        multiext(
            "manuscript/figures/flowerhorn_methods/ne1_{ne1}_founders_{founders}_t{t0}_s{seed}_n{sample_size}",
            ".pdf",
            ".png",
            ".svg",
        ),
    envmodules:
        "calcua/2024a",
        "Julia/1.11.1",
    params:
        # It's fixed for now...
        ne2_1=15_000,
        ne2_2=25_000,
    script:
        "plots/flowerhorn_methods.jl"

rule flowerhorn_simulation:
    input:
        "src/msprime/flowerhorn_simulation.py",
    output:
        expand(
            "steps/trees/flowerhorn/ne1_{{ne1}}_ne2_{{ne2}}_t{{t_inv}}_n{{sample_size}}/s{{seed}}_chr{i}.trees",
            i=range(NUM_CHROMOSOMES),
        ),
    resources:
        runtime="12h",
    threads: 4
    conda:
        "../external/conda_env.yaml"
    log:
        "logs/trees/flowerhorn/ne1_{ne1}_ne2_{ne2}_t{t_inv}_n{sample_size}/s{seed}.log",
    shell:
        """
        source {COMMON}
        python {input} {wildcards.seed} {wildcards.ne1} {wildcards.ne2} {wildcards.t_inv} {wildcards.sample_size} {output} 2> {log}
        """

rule fit_exponential_piecewise_model_boot_approx_flowerhorn:
    input:
        "src/pymc/exponential_piecewise_nuts_boot_approx.py",
        "steps/binned_ld/flowerhorn/ne1_{ne1}_ne2_{founders}_t{t_inv}_n{sample_size}/s{seed}.csv",
        "steps/inference/ballpark_ne/flowerhorn/ne1_{ne1}_ne2_{founders}_t{t_inv}_n{sample_size}/s{seed}.csv",
    output:
        "steps/inference/exponential_piecewise_model/flowerhorn/ne1_{ne1}_ne2_{founders}_t{t_inv}_n{sample_size}/s{seed}_b{boot}_advi.nc",
    resources:
        runtime="30min",
    threads: 1
    conda:
        "../external/conda_env.yaml"
    log:
        "logs/inference/exponential_piecewise_model_bagging/flowerhorn/ne1_{ne1}_ne2_{founders}_t{t_inv}_n{sample_size}/s{seed}_b{boot}_advi.log",
    params:
        ne1_prior_sd=10_000,
        t0_prior_mean=50,
        t0_prior_sd=30,
        alpha_logfold_prior_sd=1,
    shell:
        """
        source {COMMON}
        TMP_COMPILEDIR=$(mktemp -d)
        export PYTENSOR_FLAGS="compiledir=${{TMP_COMPILEDIR}}"
        python {input} \
            {params.ne1_prior_sd} {params.t0_prior_mean} \
            {params.t0_prior_sd} {params.alpha_logfold_prior_sd} \
            {wildcards.sample_size} {wildcards.seed} {wildcards.boot} {output} 2>&1 > {log}
        rm -rf "${{TMP_COMPILEDIR}}"
        """
