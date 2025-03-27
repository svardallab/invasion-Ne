rule run_hapibd:
    input:
        invcf="steps/vcfs/{name}.vcf.gz",
        inmap="steps/maps/{name}.plink.map",
        executable="external/hap-ibd.jar",
    output:
        "steps/hapibd/{name}.ibd.gz",
        "steps/hapibd/{name}.hbd.gz",
    log:
        "logs/hap-ibd/{name}.log",
    shadow:
        "minimal"
    threads: 2
    envmodules:
        "calcua/2024a",
        "Java/21.0.5",
    shell:
        """
        java -jar {input.executable} \
            gt={input.invcf} map={input.inmap} out=steps/hapibd/{wildcards.name} \
            nthreads={threads} &> /dev/null
        mv steps/hapibd/{wildcards.name}.log {log}
        """


rule post_processing_ibd:
    input:
        invcf="steps/vcfs/{name}.vcf.gz",
        inmap="steps/maps/{name}.plink.map",
        ibd="steps/hapibd/{name}.ibd.gz",
        executable="external/merge-ibd-segments.jar",
    output:
        "steps/hapibd/{name}.postprocessed.ibd",
    shadow:
        "minimal"
    params:
        gap=0.6,
        discord=1,
    envmodules:
        "calcua/2024a",
        "Java/21.0.5",
    shell:
        """
        gunzip -c {input.ibd} | \
            java -jar {input.executable} {input.invcf} {input.inmap} \
            {params.gap} {params.discord} > {output}
        """


rule ibdne:
    input:
        ibdne="external/ibdne.jar",
        ibds="steps/hapibd/{model}/s{seed}_n{n}.postprocessed.ibd",
        inmap="steps/maps/{model}/s{seed}_n{n}.plink.map",
    log:
        "logs/ibdne/{model}/s{seed}_n{n}.log",
    output:
        pair_excl="steps/ibdne/{model}/s{seed}_n{n}.pair.excl",
        region_excl="steps/ibdne/{model}/s{seed}_n{n}.region.excl",
        estimate="steps/ibdne/{model}/s{seed}_n{n}.ne",
        boot="steps/ibdne/{model}/s{seed}_n{n}.boot",
    threads: 4
    shadow:
        "minimal"
    envmodules:
        "calcua/2024a",
        "Java/21.0.5",
    shell:
        """
        cat {input.ibds} | java -jar {input.ibdne} out=out \
            seed={wildcards.seed} map={input.inmap} nthreads={threads}
        mv *.log {log}
        mv *.pair.excl {output.pair_excl}
        mv *.region.excl {output.region_excl}
        mv *.ne {output.estimate}
        mv *.boot {output.boot}
        """
