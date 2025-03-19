# Workaround CALCUA VSC requirements about conda environments and containers
COMMON = "calcua.sh"
MODELS = ["constant_demography", "invasion_demography"]
PROGRAMS = ["gone", "ibdne", "singer"]


rule all:
    input:
        expand(
            "steps/{program}/{model}_{seed}_n150.ne.csv",
            seed=[100 + i for i in range(25)],
            model=MODELS,
            program=PROGRAMS,
        ),
        expand(
            "simulations/{model}.csv",
            model=MODELS,
        ),


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
        r"steps/trees/{model}_{seed,\d+}.trees.tsz",
    shell:
        """
    source {COMMON}
        python {input} {wildcards.seed} {output}
        """


rule overlay_mutations:
    input:
        "src/overlay_mutations.py",
        "steps/trees/{model}_{seed}.trees.tsz",
    output:
        r"steps/vcfs/{model}_{seed,\d+}.vcf.gz",
        r"steps/vcfs/{model}_{seed,\d+}.vcf.gz.tbi",
    shell:
        """
    source {COMMON}
        python {input} {wildcards.seed} | bgzip -c > {output[0]}
        tabix -p vcf {output[0]}
    """


rule population_size_trajectory:
    input:
        "src/population_size_trajectory.py",
        "simulations/{file}.yaml",
    output:
        "simulations/{file}.csv",
    shell:
        """
    source {COMMON}
        python {input} > {output}
    """


rule smc_vcf2smc:
    input:
        "steps/vcfs/{model}_{seed}.vcf.gz",
    output:
        r"steps/smcpp/{model}_{seed, \d+}_n150.smc.gz",
    log:
        "steps/smcpp/{model}_{seed}_n150.log",
    container:
        "external/smcpp.sif"
    shell:
        """
        SAMPLES=$(zcat {input} | grep -m1 "^#CHROM" | cut -f10- | tr '\t' ',')
        echo $SAMPLES
        smc++ vcf2smc {input} {output} 1 Pop1:$SAMPLES
        """


rule smcpp_estimate:
    input:
        "steps/smcpp/{file}.smc.gz",
    output:
        "steps/smcpp/{file}.model.final.json",
    log:
        "steps/smcpp/{file}.model.final.log",
    container:
        "external/smcpp.sif"
    threads: 8
    shell:
        """
        smc++ estimate -c {threads} 1e-8 {input} -o steps/smcpp/{wildcards.file}_dir &> {log}
        mv steps/smcpp/{wildcards.file}_dir/model.final.json {output}
        rm -r steps/smcpp/{wildcards.file}_dir
        """


rule smcpp_ne:
    input:
        "steps/smcpp/{file}.model.final.json",
    output:
        "steps/smcpp/{file}.ne.csv",
    container:
        "external/smcpp.sif"
    shell:
        """
        # Generate the full CSV with all columns
        smc++ plot -c temp_{wildcards.file}.pdf -g 1 {input}
        awk -F, 'NR==1 {{print "generations,Ne"}} NR>1 {{print $2","$3}}' temp_{wildcards.file}.csv > {output}
        # Remove the temporary file
        rm temp_{wildcards.file}.*
        """


rule gone_ne:
    input:
        script="external/gone/script_GONE.sh",
        invcf="steps/vcfs/{model}_{seed}.vcf.gz",
    log:
        "steps/gone/{model}_{seed}_n150.log",
    output:
        "steps/gone/{model}_{seed}_n150.ne.csv",
    threads: 4
    shadow:
        "shallow"
    shell:
        """
        source {COMMON}
        export PHASE=1  # Phase = 0 (pseudohaploids), 1 (known phase), 2 (unknown phase)
        export cMMb=1   # CentiMorgans per Megabase (if distance is not available in map file)
        export DIST=1   # none (0), Haldane correction (1) or Kosambi correction (2)
        export NGEN=2000  # Number of generations for which linkage data is obtained in bins
        export NBIN=400   # Number of bins (each bin includes NGEN/NBIN = 5 generations)
        export MAF=0.0    # Minor allele frequency (recommended 0)
        export ZERO=1     # 0: Remove SNPs with zeroes (1: allow them)
        export maxNCHROM=-99  # Max number of chromosomes (-99 = all)
        export maxNSNP=50000  # Max number of SNPs per chromosome (max = 50000)
        export hc=0.05    # Max value of c analyzed (recommended 0.05; max is 0.5)
        export REPS=40    # Number of replicates to run GONE (recommended 40)
        export threads={threads}  # Number of threads (-99 uses all available processors)
        
        plink --vcf {input.invcf} --recode --out myplink &> {log}
        
        bash {input.script} myplink &>> {log}
        
        echo 'generations,Ne' > {output}
        tail -n +3 Output_Ne_myplink | sed 's/\t/,/g' >> {output}
        """


rule ibdne:
    input:
        hapibd="external/hap-ibd.jar",
        ibdne="external/ibdne.jar",
        invcf="steps/vcfs/{model}_{seed}.vcf.gz",
    log:
        "steps/ibdne/{model}_{seed}_n150.log",
    output:
        "steps/ibdne/{model}_{seed}_n150.ne.csv",
    threads: 4
    shadow:
        "shallow"
    envmodules:
        "calcua/2024a",
        "Java/21.0.5",
    shell:
        """
        echo -e "1  rs100  0  1\n1  rs101  1  1000000" > plink.map
        java -Xmx8g -jar {input.hapibd} gt={input.invcf} map=plink.map nthreads={threads} out=out
        cat out.log > {log}
        zcat out.ibd.gz | java -jar {input.ibdne} mincm=1 seed={wildcards.seed} map=plink.map nboots=0 out=ne nthreads={threads}
        ls >> {log}
        cat ne.log >> {log}
        echo 'generations,Ne' > {output}
        tail -n +2 ne.ne | sed 's/\t/,/g' >> {output}
        """


rule singer:
    input:
        singer_smk="external/singer-snakemake",
        singer_config="external/singer-snakemake/config/data_config.yaml",
        invcf="steps/vcfs/{model}_{seed}.vcf.gz",
    log:
        "steps/singer/{model}_{seed}_n150.log",
    output:
        directory("steps/singer/{model}_{seed}_n150"),
    threads: 15
    shadow:
        "shallow"
    shell:
        """
        source {COMMON}
        # Copy the singer snakemake directory
        cp -r {input.singer_smk} .
        # Symlink the VCF file
        mkdir singer-snakemake/data
        ln {input.invcf} singer-snakemake/data/
        # Execute snakemake workflow
        cd singer-snakemake
        snakemake -c {threads} &> ../{log}
        # Copy relevant output files
        cd ..
        mv singer-snakemake/results/{wildcards.model}_{wildcards.seed}/trees {output}
        """


rule singer_ne:
    input:
        script="src/singer_ne.py",
        indir="steps/singer/{model}_{seed}_n150",
    output:
        "steps/singer/{model}_{seed}_n150.ne.csv",
    shell:
        """
        source {COMMON}
        python {input.script} {input.indir}/*.trees > {output}
        """
