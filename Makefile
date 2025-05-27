SHELL:=/bin/bash
# Everything that can be installed via conda is containerized here (Calcua VSC limitations)
CONDA_ENV_PREFIX = external/conda_env
CONDA_ENV_YAML = external/conda_env.yaml

# GONE2 compilation settings
GONE2_DIR = external/GONE2/
# See https://github.com/esrud/GONE2/issues/4
GONE2_HASH = 81915701aef063fdaab6e9e7f625098888b7b669
GONE2_BIN = external/gone2

# Compilation flags for increased loci
MAXLOCI = 10000000
MAXIND = 200

# Download java executables from Browning lab
IBDNE_URL = https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
IBDNE_BIN = external/ibdne.jar
HAP_IBD_URL = https://faculty.washington.edu/browning/hap-ibd.jar
HAP_IBD_BIN = external/hap-ibd.jar
MERGE_IBD_URL = https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar
MERGE_IBD_BIN = external/merge-ibd-segments.jar

.PHONY: deps clean gone lint run

deps: $(CONDA_ENV_PREFIX) $(GONE2_BIN) $(IBDNE_BIN) $(HAP_IBD_BIN) $(MERGE_IBD_BIN) $(BINLD_BIN)

run:
	snakemake -c$(cores) --use-envmodules --sdm apptainer

dry-run:
	snakemake -n -c$(cores) --use-envmodules --sdm apptainer

lint:
	./external/conda_env/bin/black src/*.py
	./external/conda_env/bin/snakefmt Snakefile

# Dependencies
$(CONDA_ENV_PREFIX): $(CONDA_ENV_YAML)
	module load hpc-container-wrapper/0.4.0; \
	conda-containerize new --prefix $(CONDA_ENV_PREFIX) $(CONDA_ENV_YAML)

$(GONE2_BIN): $(GONE2_DIR)
	module purge && \
	module load calcua/2024a && \
	module load foss/2024a && \
	cd $(GONE2_DIR) && \
	git checkout $(GONE2_HASH) && \
    make MAXLOCI=$(MAXLOCI) MAXIND=$(MAXIND) gone && \
    mv gone2 ..

$(HAP_IBD_BIN):
	wget $(HAP_IBD_URL) -O $(HAP_IBD_BIN)

$(IBDNE_BIN):
	wget $(IBDNE_URL) -O $(IBDNE_BIN)

$(MERGE_IBD_BIN):
	wget $(MERGE_IBD_URL) -O $(MERGE_IBD_BIN)

# Compile Rust code
BINLD_DIR = external/ld_binning_src/
BINLD_BIN = external/ld_binning

$(BINLD_BIN): $(BINLD_DIR)
	cd $(BINLD_DIR) && \
	module load Rust && \
	module load Clang && \
	module load Perl && \
	export CARGO_HOME=$(mktemp -d /tmp/cargo-home.XXXXXX) && \
	cargo build --release && \
	cp target/release/ld_binning ../../$(BINLD_BIN)

clean:
	rm -rf $(CONDA_ENV_PREFIX) $(GONE2_BIN) $(GONE2_DIR)
