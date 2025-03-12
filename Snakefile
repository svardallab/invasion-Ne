# Workaround CALCUA VSC requirements about conda environments and containers
COMMON = "calcua.sh"


rule all:
    input:
        expand(
            "steps/smcpp/invasion_demography_{seed}_n25.ne.csv",
            seed=[100 + i for i in range(15)],
        ),
        expand(
            "steps/gone2/invasion_demography_{seed}_n100.ne.csv",
            seed=[100 + i for i in range(15)],
        ),


rule model:
    input:
        "src/invasion_demography.py",
    output:
        "simulations/invasion_demography.pdf",
        "simulations/invasion_demography.yaml",
    shell:
        """
    source {COMMON}
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
    source {COMMON}
        python {input} {wildcards.seed} {output}
        """


rule overlay_mutations:
    input:
        "src/overlay_mutations.py",
        "steps/trees/invasion_demography_{seed}.trees.tsz",
    output:
        r"steps/vcfs/invasion_demography_{seed,\d+}.vcf.gz",
        r"steps/vcfs/invasion_demography_{seed,\d+}.vcf.gz.tbi",
    shell:
        """
    source {COMMON}
        python {input} {wildcards.seed} | bgzip -c > {output[0]}
        tabix -p vcf {output[0]}
    """


rule smc_vcf2smc:
    input:
        "steps/vcfs/invasion_demography_{seed}.vcf.gz",
    output:
        r"steps/smcpp/invasion_demography_{seed, \d+}_n{n,\d+}.smc.gz",
    log:
        "steps/smcpp/invasion_demography_{seed}_n{n}.log",
    container:
        "external/smcpp.sif"
    shell:
        """
    	SAMPLES=$(printf "tsk_%dindv," $(seq 0 $(({wildcards.n}-1))) | sed 's/,$//')
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


rule gone2_ne:
    input:
        "steps/vcfs/invasion_demography_{seed}.vcf.gz",
    output:
        "steps/gone2/invasion_demography_{seed}_n100.ne.csv",
    threads: 4
    shell:
        """
        # Initialize the output file with headers
        if [ ! -f {output} ]; then
            echo "generations,Ne" > {output}
        fi
    temp="steps/gone2/temp_{wildcards.seed}.vcf"
        gunzip -c {input} > $temp
        for i in {{1..5}}
        do
            seed=$(( {wildcards.seed} + i - 1 ))  # Increment seed value by the loop index
            ./external/gone2 -e -s ${{seed}} -t {threads} -r 1.0 $temp -o steps/gone2/temp_{wildcards.seed}
            awk -F, 'NR > 1 {{print $1","$2}}' steps/gone2/temp_{wildcards.seed}_GONE2_Ne >> {output}
            rm steps/gone2/temp_{wildcards.seed}_GONE2_*
        done
    rm $temp
        """
