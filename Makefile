# Everything that can be installed via conda is containerized here (Calcua VSC limitations)
CONDA_ENV_PREFIX = external/conda_env
CONDA_ENV_YAML = external/conda_env.yaml

# GONE2 compilation settings
GONE2_REPO = https://github.com/esrud/GONE2
GONE2_DIR = external/GONE2
GONE2_BIN = external/gone2

# Compilation flags for increased loci
MAXLOCI = 4000000
MAXIND = 1000

# Download java executables from Browning lab
IBDNE_URL = https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
IBDNE_BIN = external/ibdne.jar
HAP_IBD_URL = https://faculty.washington.edu/browning/hap-ibd.jar
HAP_IBD_BIN = external/hap-ibd.jar
MERGE_IBD_URL = https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar
MERGE_IBD_BIN = external/merge-ibd-segments.jar

.PHONY: deps clean gone lint run

deps: $(CONDA_ENV_PREFIX) $(GONE2_BIN) $(IBDNE_BIN) $(HAP_IBD_BIN) $(MERGE_IBD_BIN)

run:
	snakemake -c$(cores) --use-envmodules --sdm apptainer

dry-run:
	snakemake -n -c$(cores) --use-envmodules --sdm apptainer

lint:
	./external/conda_env/bin/black src/*.py
	./external/conda_env/bin/snakefmt Snakefile

# Dependencies
$(CONDA_ENV_PREFIX): $(CONDA_ENV_YAML)
	conda-containerize new --prefix $(CONDA_ENV_PREFIX) $(CONDA_ENV_YAML)

$(GONE2_BIN): $(GONE2_DIR)
	cd $(GONE2_DIR) && \
    make MAXLOCI=$(MAXLOCI) MAXIND=$(MAXIND) gone && \
    mv gone2 ..

$(GONE2_DIR):
	mkdir -p external
	git clone $(GONE2_REPO) $(GONE2_DIR)

$(HAP_IBD_BIN):
	wget $(HAP_IBD_URL) -O $(HAP_IBD_BIN)

$(IBDNE_BIN):
	wget $(IBDNE_URL) -O $(IBDNE_BIN)

$(MERGE_IBD_BIN):
	wget $(MERGE_IBD_URL) -O $(MERGE_IBD_BIN)

clean: 
	rm -rf $(CONDA_ENV_PREFIX) $(GONE2_BIN) $(GONE2_DIR)