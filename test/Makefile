SCRIPT = ../src/InPheRNo_simple.py
RUN_DIR = ./run_dir
RUN_FILES = ../data/run_files
RESULTS_DIR = $(RUN_DIR)/results
INPUT_DATA_DIR = ../data

RMR = rm -rf
MKDIR_P = mkdir -p

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)

# ----------------------------------------------------------------
# Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_InPheRNo_simplified:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file InPheRNo_simplified_TEMPLATE.yml

