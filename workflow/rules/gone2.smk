rule gone2:
    input:
        invcf="steps/vcfs/{model}/s{seed}_n{n}.vcf.gz",
        inmap="steps/maps/{model}/s{seed}_n{n}.plink.map",
        executable="external/gone2"
    output:
        estimate="steps/gone2/{model}/s{seed}_n{n}_GONE2_Ne",
        d2="steps/gone2/{model}/s{seed}_n{n}_GONE2_d2",
        stats="steps/gone2/{model}/s{seed}_n{n}_GONE2_STATS"
    log:
        "log/gone2/{model}/s{seed}_n{n}.log"
    threads: 4
    shadow: "shallow"
    shell:
        """
        zcat {input.invcf} > out.vcf
        cp {input.inmap} out.map
        {input.executable} -r 1.0 -t {threads} -g 2 out.vcf &> {log}
        mv *_GONE2_Ne {output.estimate}
        mv *_GONE2_d2 {output.d2}
        mv *_GONE2_STATS {output.stats}
        """