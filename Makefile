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

.PHONY: deps clean gone

deps: $(CONDA_ENV_PREFIX) $(SMCPP_FILE) $(IBDNE_BIN) $(HAP_IBD_BIN) gone

$(CONDA_ENV_PREFIX): $(CONDA_ENV_YAML)
	conda-containerize new --prefix $(CONDA_ENV_PREFIX) $(CONDA_ENV_YAML)

$(SMCPP_FILE):
	apptainer build --fakeroot $(SMCPP_FILE) $(SMCPP_DOCKER)

$(HAP_IBD_BIN):
	wget $(HAP_IBD_URL) -O $(HAP_IBD_BIN)

$(IBDNE_BIN):
	wget $(IBDNE_URL) -O $(IBDNE_BIN)

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
