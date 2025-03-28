rule run_hapne_ld:
    input:
        invcfs = expand(
            "steps/vcfs/{{name}}_{chrom}.vcf.gz",
            chrom=[f"contig{i+1}" for i in range(NUM_CHR)],
        ),
        inmaps = expand(
            "steps/maps/{{name}}_{chrom}.shapeit.map",
            chrom=[f"contig{i+1}" for i in range(NUM_CHR)],
        ),
        dir="external/hapne-snakemake",
    output:
        table="steps/hapne/{name}_ld_hapne_estimate.csv",
        summary="steps/hapne/{name}_ld_hapne_summary.txt",
        residuals="steps/hapne/{name}_ld_hapne_residuals.png",
        popsize="steps/hapne/{name}_ld_hapne_pop_trajectory.png",
    log:
        "logs/hapne_ld/{name}.log",
    shadow: "shallow"
    threads: 8
    shell:
        """
        source {COMMON}
        cp -r {input.dir} hapne_analysis
        mkdir hapne_analysis/data
        ln {input.invcfs} hapne_analysis/data
        ln {input.inmaps} hapne_analysis/data
        cd hapne_analysis/data 

        toremove=$(basename {wildcards.name})_
        for file in "$toremove"*; do
            new_name="${{file#$toremove}}"
            if [[ "$file" != "$new_name" ]]; then
                mv "$file" "$new_name"
            fi
        done
        
        cd ..
        snakemake -c {threads} --config \
            data_dir='data/' out_dir='results/' \
            map_file_suffix='.shapeit.map' \
            method='ld'  2> ../{log}
        ls .
        cd ..
        ls .
        mv hapne_analysis/results/*_estimate.csv {output.table}
        mv hapne_analysis/results/*_summary.txt {output.summary}
        mv hapne_analysis/results/*_residuals.png {output.residuals}
        mv hapne_analysis/results/*_trajectory.png {output.popsize}
        """

rule run_hapne_ibd:
    input:
        invcfs = expand(
            "steps/vcfs/{{name}}_{chrom}.vcf.gz",
            chrom=[f"contig{i+1}" for i in range(NUM_CHR)],
        ),
        inmaps = expand(
            "steps/maps/{{name}}_{chrom}.plink.map",
            chrom=[f"contig{i+1}" for i in range(NUM_CHR)],
        ),
        build="steps/regions/{name}_region.txt",
        dir="external/hapne-snakemake",
    output:
        table="steps/hapne/{name}_ibd_hapne_estimate.csv",
        summary="steps/hapne/{name}_ibd_hapne_summary.txt",
        residuals="steps/hapne/{name}_ibd_hapne_residuals.png",
        popsize="steps/hapne/{name}_ibd_hapne_pop_trajectory.png",
    log:
        "logs/hapne_ibd/{name}.log",
    shadow: "shallow"
    envmodules:
        "calcua/2024a",
        "Java/21.0.5",
    threads: 8
    shell:
        """
        source {COMMON}
        cp -r {input.dir} hapne_analysis
        mkdir hapne_analysis/data
        ln {input.invcfs} hapne_analysis/data
        ln {input.inmaps} hapne_analysis/data
        ln {input.build} hapne_analysis/data/build.txt
        cd hapne_analysis/data 

        toremove=$(basename {wildcards.name})_
        for file in "$toremove"*; do
            new_name="${{file#$toremove}}"
            if [[ "$file" != "$new_name" ]]; then
                mv "$file" "$new_name"
            fi
        done
        
        cd ..
        snakemake -c {threads} --config \
            data_dir='data/' out_dir='results/' \
            map_file_suffix='.plink.map' \
            method='ibd' genome_build='data/build.txt' 2> ../{log}
        ls .
        ls results
        cd ..
        mv hapne_analysis/results/*_estimate.csv {output.table}
        mv hapne_analysis/results/*_summary.txt {output.summary}
        mv hapne_analysis/results/*_residuals.png {output.residuals}
        mv hapne_analysis/results/*_trajectory.png {output.popsize}
        """
