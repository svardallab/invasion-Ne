rule all:
    input:
        expand(
            "steps/trees/invasion_demography_{seed}.trees.tsz",
            seed=[100 + i for i in range(100)],
        ),


rule model:
    input:
        "src/invasion_demography.py",
    output:
        "simulations/invasion_demography.pdf",
        "simulations/invasion_demography.yaml",
    shell:
        """
        python {input} simulations/invasion_demography
        """


rule coalescence_simulation:
    input:
        "src/run_msprime.py",
        "simulations/invasion_demography.yaml",
    output:
        r"steps/trees/invasion_demography_{seed,\d+}.trees.tsz",
    shell:
        """
        python {input} {wildcards.seed} {output}
        """
