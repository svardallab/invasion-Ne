
rule model:
    input:
        "src/{model}.py",
    output:
        "simulations/{model}.pdf",
        "simulations/{model}.yaml",
    shell:
        """
        source {COMMON}
        python {input} simulations/{wildcards.model}
        """


rule coalescence_simulation:
    input:
        "src/run_msprime.py",
        "simulations/{model}.yaml",
    output:
        expand(
            r"steps/trees/{{model}}/s{{seed,\d+}}_n{{n,\d+}}_{chrom}.trees.tsz",
            chrom=[f"chr{i+1}" for i in range(NUM_CHR)],
        ),
    shell:
        """
        source {COMMON}
        python {input} {wildcards.seed} {wildcards.n} {output}
        """


rule overlay_mutations:
    input:
        "src/overlay_mutations.py",
        "steps/trees/{model}/s{seed}_n{n}_{chrom}.trees.tsz",
    output:
        "steps/vcfs/{model}/s{seed}_n{n}_{chrom}.vcf.gz",
        "steps/vcfs/{model}/s{seed}_n{n}_{chrom}.vcf.gz.tbi",
    shell:
        """
        source {COMMON}
        python {input} {wildcards.seed} | bgzip -c > {output[0]}
        tabix -p vcf {output[0]}
        """

rule plink_genetic_map:
    input:
        "src/plink_genetic_map.py",
        "steps/trees/{model}/s{seed}_n{n}_{chrom}.trees.tsz",
    output:
        "steps/maps/{model}/s{seed}_n{n}_{chrom}.plink.map",
    shell:
        """
        source {COMMON}
        python {input} {wildcards.chrom} > {output}
        """

rule shapeit_genetic_map:
    input:
        "src/shapeit_genetic_map.py",
        "steps/trees/{model}/s{seed}_n{n}_{chrom}.trees.tsz",
    output:
        "steps/maps/{model}/s{seed}_n{n}_{chrom}.shapeit.map",
    shell:
        """
        source {COMMON}
        python {input} {wildcards.chrom} > {output}
        """
