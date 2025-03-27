# Everything that can be installed via conda is containerized here (Calcua VSC limitations)
CONDA_ENV_PREFIX = external/conda_env
CONDA_ENV_YAML = external/conda_env.yaml
# Maintained smcpp executable is a docker image we have to convert into apptainer
SMCPP_DOCKER = docker://terhorst/smcpp
SMCPP_FILE = external/smcpp.sif
# Download java executables from Browning lab
IBDNE_URL = https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
IBDNE_BIN = external/ibdne.jar
HAP_IBD_URL = https://faculty.washington.edu/browning/hap-ibd.jar
HAP_IBD_BIN = external/hap-ibd.jar
# Custom version of gone 
GONE_DIR = external/gone/PROGRAMMES
GONE_BINS = $(GONE_DIR)/LD_SNP_REAL3 $(GONE_DIR)/SUMM_REP_CHROM3 $(GONE_DIR)/MANAGE_CHROMOSOMES2 $(GONE_DIR)/GONE $(GONE_DIR)/GONEaverage
# Relate
RELATE_STATIC_TAR = external/relate_v1.1.9_x86_64_static.tgz
RELATE_STATIC_BIN = external/relate_v1.1.9_x86_64_static/bin/Relate

.PHONY: deps clean gone lint

deps: $(CONDA_ENV_PREFIX) $(SMCPP_FILE) $(IBDNE_BIN) $(HAP_IBD_BIN) gone $(RELATE_STATIC_DIR)

run:
	snakemake -c$(cores) --use-envmodules --sdm apptainer


dry-run:
	snakemake -n -c$(cores) --use-envmodules --sdm apptainer

lint:
	./external/conda_env/bin/black src/*.py
	./external/conda_env/bin/snakefmt Snakefile

$(CONDA_ENV_PREFIX): $(CONDA_ENV_YAML)
	conda-containerize new --prefix $(CONDA_ENV_PREFIX) $(CONDA_ENV_YAML)

$(SMCPP_FILE):
	apptainer build --fakeroot $(SMCPP_FILE) $(SMCPP_DOCKER)

$(HAP_IBD_BIN):
	wget $(HAP_IBD_URL) -O $(HAP_IBD_BIN)

$(IBDNE_BIN):
	wget $(IBDNE_URL) -O $(IBDNE_BIN)

$(RELATE_STATIC_DIR): $(RELATE_STATIC_TAR)
	tar -xvzf $(RELATE_STATIC_TAR) -C external/

# Compile gone
gone: $(GONE_BINS)

$(GONE_DIR)/LD_SNP_REAL3: $(GONE_DIR)/LD_SNP_REAL3.c $(GONE_DIR)/genlib.c
	cd $(GONE_DIR) && cc -o LD_SNP_REAL3 -w -O -static LD_SNP_REAL3.c genlib.c -lm

$(GONE_DIR)/SUMM_REP_CHROM3: $(GONE_DIR)/SUMM_REP_CHROM3.c $(GONE_DIR)/genlib.c
	cd $(GONE_DIR) && cc -o SUMM_REP_CHROM3 -w -O -static SUMM_REP_CHROM3.c genlib.c -lm

$(GONE_DIR)/MANAGE_CHROMOSOMES2: $(GONE_DIR)/MANAGE_CHROMOSOMES2.c $(GONE_DIR)/genlib.c
	cd $(GONE_DIR) && cc -o MANAGE_CHROMOSOMES2 -w -O -static MANAGE_CHROMOSOMES2.c genlib.c -lm

$(GONE_DIR)/GONE: $(GONE_DIR)/GONE.cpp
	cd $(GONE_DIR) && g++ -std=gnu++0x -static GONE.cpp -o GONE

$(GONE_DIR)/GONEaverage: $(GONE_DIR)/GONEaverage.cpp
	cd $(GONE_DIR) && g++ -std=gnu++0x -static GONEaverage.cpp -o GONEaverage

clean: 
	rm -rf $(CONDA_ENV_PREFIX) $(SMCPP_FILE)
	rm -f $(GONE_BINS)
	rm -rf $(RELATE_STATIC_DIR)
