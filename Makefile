CONDA_ENV_PREFIX = external/conda_env
CONDA_ENV_YAML = external/conda_env.yaml

SMCPP_DOCKER = docker://terhorst/smcpp
SMCPP_FILE = external/smcpp.sif

GONE2_URL = https://github.com/esrud/GONE2/releases/download/v1.0.1/GONE2_v1.0.1.tar.gz
GONE2_BIN = external/gone2

SINGER_URL = https://github.com/popgenmethods/SINGER/raw/refs/heads/main/releases/singer-0.1.8-beta-linux-x86_64.tar.gz
SINGER_BIN = external/singer-0.1.8-beta-linux-x86_64

.PHONY: deps clean

deps: $(CONDA_ENV_PREFIX) $(SMCPP_FILE) $(GONE2_BIN) $(SINGER_BIN)

$(CONDA_ENV_PREFIX):$(CONDA_ENV_YAML)
	conda-containerize new --prefix $(CONDA_ENV_PREFIX) $(CONDA_ENV_YAML)

$(SMCPP_FILE):
	apptainer build --fakeroot $(SMCPP_FILE) $(SMCPP_DOCKER)

$(GONE2_BIN): 
	wget $(GONE2_URL) -O external/GONE2_v1.0.1.tar.gz
	tar -xzf external/GONE2_v1.0.1.tar.gz -C external
	rm external/GONE2_v1.0.1.tar.gz

$(SINGER_BIN): 
	wget $(SINGER_URL) -O external/singer-0.1.8-beta-linux-x86_64.tar.gz
	tar -xzf external/singer-0.1.8-beta-linux-x86_64.tar.gz -C external
	rm external/singer-0.1.8-beta-linux-x86_64.tar.gz
clean:
	rm -rf $(CONDA_ENV_PREFIX) $(SMCPP_FILE) $(GONE2_BIN)

