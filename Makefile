CONDA_ENV_PREFIX = external/conda_env
CONDA_ENV_YAML = external/conda_env.yaml

SMCPP_DOCKER = docker://terhorst/smcpp
SMCPP_FILE = external/smcpp.sif

IBDNE_URL = https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
IBDNE_BIN = external/ibdne.jar
HAP_IBD_URL = https://faculty.washington.edu/browning/hap-ibd.jar
HAP_IBD_BIN = external/hap-ibd.jar

.PHONY: deps clean

deps: $(CONDA_ENV_PREFIX) $(SMCPP_FILE) $(SINGER_BIN) $(IBDNE_BIN) $(HAP_IBD_BIN)

$(CONDA_ENV_PREFIX):$(CONDA_ENV_YAML)
	conda-containerize new --prefix $(CONDA_ENV_PREFIX) $(CONDA_ENV_YAML)

$(SMCPP_FILE):
	apptainer build --fakeroot $(SMCPP_FILE) $(SMCPP_DOCKER)

$(HAP_IBD_BIN): 
	wget $(HAP_IBD_URL) -O $(HAP_IBD_BIN)
$(IBDNE_BIN): 
	wget $(IBDNE_URL) -O $(IBDNE_BIN)


clean:
	rm -rf $(CONDA_ENV_PREFIX) $(SMCPP_FILE) $(GONE2_BIN)

