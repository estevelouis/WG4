# ------------------------------
# DIVERSITY MEASUREMENT MAKEFILE
# ------------------------------
#
# To compile, run:
# 	make measurement_full
#

# -----------------
# COMPILING OPTIONS
# -----------------

CC = gcc
#CC = clang
CCPP = g++
#CCPP = clang++
C_VERSION = -std=c99
CPP_VERSION = -std=c++11
CFLAGS_PEDANTIC = -Wall -Wextra -Werror -Wformat -Wformat-security -pedantic -march=native
CFLAGS_STANDARD = -Wall -Wextra -Wformat -Wformat-security -march=native
CFLAGS_DEBUG = -Wall -Wextra -Wformat -Wformat-security -g3
CFLAGS_OPTIM = -Wall -Wextra -Wformat -Wformat-security -march=native -ffast-math
CFLAGS_TEST = -Wall -Wextra -Wformat -Wformat-security -march=native
CFLAGS_COVERAGE = -Wall -Wextra -Wformat -Wformat-security -march=native --coverage
CFLAGS_PROFILING = -Wall -Wextra -Wformat -Wformat-security -march=native -pg -g
CPPFLAGS = -Isrc/include -Iudpipe/src_lib_only -I/usr/include -I/usr/include/x86_64-linux-gnu/c++/11 -I/usr/include/c++/11
CPPFLAGS_TEST = $(CPPFLAGS) -Itest/include
PATH_LIB_DIVERSUTILS = $(HOME)/.local/lib/diversutils
LDFLAGS_MINIMUM = -L$(PATH_LIB_DIVERSUTILS) -L/usr/lib/x84_64-linux-gnu
ifeq ($(origin CONDA_PREFIX), undefined)
	LDFLAGS = $(LDFLAGS_MINIMUM)
else
	LDFLAGS = -L$(CONDA_PREFIX)/lib $(LDFLAGS_MINIMUM)
endif
OPT_LEVEL = -O2
OPT_LEVEL_DEBUG = -O0
OPT_LEVEL_OPTIM = -O3
OPT_LEVEL_TEST = -O2
OPT_LEVEL_LINK = -O3
PYTHON_INCLUDE_FLAG = -I$(HOME)/.local/conda/include/python3.11
LINKER_FLAGS = -lc -lm -lrt -pthread
LINKER_FLAGS_EXTRA = -lstdc++ -ludpipe

# ----------------------
# UPDATE LD_LIBRARY_PATH
# ----------------------

export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):$(PATH_LIB_DIVERSUTILS)

# ----------------
# MACRO PARAMETERS
# ----------------

NUM_ROW_THREADS = -1
NUM_MATRIX_THREADS = -1
ROW_GENERATION_BATCH_SIZE = -1

NUM_FILE_READING_THREADS = 4

ENABLE_NON_DISPARITY_MULTITHREADING = 1

ENABLE_AVX256 = 0

ENABLE_DISPARITY_FUNCTIONS = 0

ENABLE_STIRLING = 0
ENABLE_RICOTTA_SZEIDL = 0
ENABLE_PAIRWISE = 0
ENABLE_LEXICOGRAPHIC = 0
ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY = 0
ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY = 0
ENABLE_LEINSTER_COBBOLD_DIVERSITY = 0
ENABLE_FUNCTIONAL_EVENNESS = 0
ENABLE_FUNCTIONAL_DISPERSION = 0
ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED = 0

ENABLE_NON_DISPARITY_FUNCTIONS = 1

ENABLE_SHANNON_WEAVER_ENTROPY = 1
ENABLE_GOOD_ENTROPY = 1
ENABLE_RENYI_ENTROPY = 1
ENABLE_PATIL_TAILLIE_ENTROPY = 0
ENABLE_Q_LOGARITHMIC_ENTROPY = 0
ENABLE_SIMPSON_INDEX = 0
ENABLE_SIMPSON_DOMINANCE_INDEX = 0
ENABLE_HILL_NUMBER_STANDARD = 0
ENABLE_HILL_EVENNESS = 0
ENABLE_BERGER_PARKER_INDEX = 0
ENABLE_JUNGE1994_PAGE20 = 0
ENABLE_JUNGE1994_PAGE22 = 0
ENABLE_BRILLOUIN_DIVERSITY = 0
ENABLE_MCINTOSH_INDEX = 0
ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975 = 0
ENABLE_SW_E_HEIP = 0
ENABLE_SW_E_ONE_MINUS_D = 0
ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964 = 0
ENABLE_SW_E_MINUS_LN_D_PIELOU1977 = 0
ENABLE_SW_F_2_1_ALATALO1981 = 0
ENABLE_SW_G_2_1_MOLINARI1989 = 0
ENABLE_SW_E_BULLA1994 = 0
ENABLE_SW_O_BULLA1994 = 0
ENABLE_SW_E_MCI_PIELOU1969 = 0
ENABLE_SW_E_PRIME_CAMARGO1993 = 0
ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL = 0

ENABLE_TIMINGS = 0
ENABLE_ITERATIVE_DISTANCE_COMPUTATION = 0
ENABLE_MULTITHREADED_ROW_GENERATION = 1
ENABLE_MULTITHREADED_MATRIX_GENERATION = 1

ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP = 1
SENTENCE_COUNT_RECOMPUTE_STEP = 10000
SENTENCE_RECOMPUTE_STEP_USE_LOG10 = 1
SENTENCE_COUNT_RECOMPUTE_STEP_LOG10 = 1.0
ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP = 1
DOCUMENT_COUNT_RECOMPUTE_STEP = 1
DOCUMENT_RECOMPUTE_STEP_USE_LOG10 = 1
DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10 = 0.5

INPUT_PATH = \"measurement_files.txt\"
OUTPUT_PATH = \"measurement_output.tsv\"
OUTPUT_PATH_TIMING = \"measurement_output_timing.tsv\"
OUTPUT_PATH_MEMORY = \"measurement_output_memory.tsv\"

#W2V_PATH = \"path/to/a/word2vec/formatted/file\"
W2V_PATH = \"/home/esteve/Documents/thesis/other_repos/word2vec/bin/MWE_S2S_IT_11GB_100d_skip-gram.bin\"

#TARGET_COLUMN = UD_FORM
TARGET_COLUMN = UD_LEMMA
#TARGET_COLUMN = UD_MWE

# JSONL only
TOKENIZATION_METHOD = 0 # REGEX
#TOKENIZATION_METHOD = 1 # UDPipe
#TOKENIZATION_METHOD = 2 # UDPipe embedded

CPP_MACRO_MULTITHREADING = -DNUM_ROW_THREADS=$(NUM_ROW_THREADS) -DNUM_MATRIX_THREADS=$(NUM_MATRIX_THREADS) -DROW_GENERATION_BATCH_SIZE=$(ROW_GENERATION_BATCH_SIZE) -DNUM_FILE_READING_THREADS=$(NUM_FILE_READING_THREADS) -DENABLE_ITERATIVE_DISTANCE_COMPUTATION=$(ENABLE_ITERATIVE_DISTANCE_COMPUTATION) -DENABLE_MULTITHREADED_ROW_GENERATION=$(ENABLE_MULTITHREADED_ROW_GENERATION) -DENABLE_MULTITHREADED_MATRIX_GENERATION=$(ENABLE_MULTITHREADED_MATRIX_GENERATION) -DENABLE_NON_DISPARITY_MULTITHREADING=$(ENABLE_NON_DISPARITY_MULTITHREADING)

CPP_MACRO_AVX = -DENABLE_AVX256=$(ENABLE_AVX256)

CPP_MACRO_DISPARITY = -DENABLE_DISPARITY_FUNCTIONS=$(ENABLE_DISPARITY_FUNCTIONS) -DENABLE_STIRLING=$(ENABLE_STIRLING) -DENABLE_RICOTTA_SZEIDL=$(ENABLE_RICOTTA_SZEIDL) -DENABLE_PAIRWISE=$(ENABLE_PAIRWISE) -DENABLE_LEXICOGRAPHIC=$(ENABLE_LEXICOGRAPHIC) -DENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY=$(ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY) -DENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY=$(ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY) -DENABLE_LEINSTER_COBBOLD_DIVERSITY=$(ENABLE_LEINSTER_COBBOLD_DIVERSITY) -DENABLE_FUNCTIONAL_EVENNESS=$(ENABLE_FUNCTIONAL_EVENNESS) -DENABLE_FUNCTIONAL_DISPERSION=$(ENABLE_FUNCTIONAL_DISPERSION) -DENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED=$(ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED)

CPP_MACRO_NON_DISPARITY = -DENABLE_NON_DISPARITY_FUNCTIONS=$(ENABLE_NON_DISPARITY_FUNCTIONS) -DENABLE_SHANNON_WEAVER_ENTROPY=$(ENABLE_SHANNON_WEAVER_ENTROPY) -DENABLE_GOOD_ENTROPY=$(ENABLE_GOOD_ENTROPY) -DENABLE_RENYI_ENTROPY=$(ENABLE_RENYI_ENTROPY) -DENABLE_PATIL_TAILLIE_ENTROPY=$(ENABLE_PATIL_TAILLIE_ENTROPY) -DENABLE_Q_LOGARITHMIC_ENTROPY=$(ENABLE_Q_LOGARITHMIC_ENTROPY) -DENABLE_SIMPSON_INDEX=$(ENABLE_SIMPSON_INDEX) -DENABLE_SIMPSON_DOMINANCE_INDEX=$(ENABLE_SIMPSON_DOMINANCE_INDEX) -DENABLE_HILL_NUMBER_STANDARD=$(ENABLE_HILL_NUMBER_STANDARD) -DENABLE_HILL_EVENNESS=$(ENABLE_HILL_EVENNESS) -DENABLE_BERGER_PARKER_INDEX=$(ENABLE_BERGER_PARKER_INDEX) -DENABLE_JUNGE1994_PAGE20=$(ENABLE_JUNGE1994_PAGE20) -DENABLE_JUNGE1994_PAGE22=$(ENABLE_JUNGE1994_PAGE22) -DENABLE_BRILLOUIN_DIVERSITY=$(ENABLE_BRILLOUIN_DIVERSITY) -DENABLE_MCINTOSH_INDEX=$(ENABLE_MCINTOSH_INDEX) -DENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975=$(ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975) -DENABLE_SW_E_HEIP=$(ENABLE_SW_E_HEIP) -DENABLE_SW_E_ONE_MINUS_D=$(ENABLE_SW_E_ONE_MINUS_D) -DENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964=$(ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964) -DENABLE_SW_E_MINUS_LN_D_PIELOU1977=$(ENABLE_SW_E_MINUS_LN_D_PIELOU1977) -DENABLE_SW_F_2_1_ALATALO1981=$(ENABLE_SW_F_2_1_ALATALO1981) -DENABLE_SW_G_2_1_MOLINARI1989=$(ENABLE_SW_G_2_1_MOLINARI1989) -DENABLE_SW_E_BULLA1994=$(ENABLE_SW_E_BULLA1994) -DENABLE_SW_O_BULLA1994=$(ENABLE_SW_O_BULLA1994) -DENABLE_SW_E_MCI_PIELOU1969=$(ENABLE_SW_E_MCI_PIELOU1969) -DENABLE_SW_E_PRIME_CAMARGO1993=$(ENABLE_SW_E_PRIME_CAMARGO1993) -DENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL=$(ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL)

CPP_MACRO_TIMING = -DENABLE_TIMINGS=$(ENABLE_TIMINGS)

CPP_MACRO_RECOMPUTE = -DENABLE_SENTENCE_COUNT_RECOMPUTE_STEP=$(ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP) -DSENTENCE_COUNT_RECOMPUTE_STEP=$(SENTENCE_COUNT_RECOMPUTE_STEP) -DSENTENCE_RECOMPUTE_STEP_USE_LOG10=$(SENTENCE_RECOMPUTE_STEP_USE_LOG10) -DSENTENCE_COUNT_RECOMPUTE_STEP_LOG10=$(SENTENCE_COUNT_RECOMPUTE_STEP_LOG10) -DENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP=$(ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP) -DDOCUMENT_COUNT_RECOMPUTE_STEP=$(DOCUMENT_COUNT_RECOMPUTE_STEP) -DDOCUMENT_RECOMPUTE_STEP_USE_LOG10=$(DOCUMENT_RECOMPUTE_STEP_USE_LOG10) -DDOCUMENT_COUNT_RECOMPUTE_STEP_LOG10=$(DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10)

CPP_MACRO_IO = -DINPUT_PATH=$(INPUT_PATH) -DOUTPUT_PATH=$(OUTPUT_PATH) -DOUTPUT_PATH_TIMING=$(OUTPUT_PATH_TIMING) -DOUTPUT_PATH_MEMORY=$(OUTPUT_PATH_MEMORY) -DW2V_PATH=$(W2V_PATH)

CPP_MACRO_OTHER = -DTARGET_COLUMN=$(TARGET_COLUMN) -DTOKENIZATION_METHOD=$(TOKENIZATION_METHOD)

CPP_MACROS = $(CPP_MACRO_MULTITHREADING) $(CPP_MACRO_AVX) $(CPP_MACRO_DISPARITY) $(CPP_MACRO_NON_DISPARITY) $(CPP_MACRO_TIMING) $(CPP_MACRO_RECOMPUTE) $(CPP_MACRO_IO) $(CPP_MACRO_OTHER)

# JSONL to W2V format variables

TARGET_NUMBER_DOCUMENTS = 1000000

CPP_MACROS_JSONL = -DTARGET_NUMBER_DOCUMENTS=$(TARGET_NUMBER_DOCUMENTS)

################

ensure_bin_dir:
	mkdir -p bin

ensure_build_dir:
	mkdir -p build

ensure_updated_udpipe:
	git submodule update --init --recursive udpipe

ensure_prepared_udpipe: ensure_updated_udpipe
	set -eu
	make --directory=udpipe/src
	make --directory=udpipe/src_lib_only

share_udpipe:
	mkdir -p $(PATH_LIB_DIVERSUTILS)
	$(CCPP) $(CPP_VERSION) $(CFLAGS_OPTIM) $(CPP_FLAGS) $(OPT_LEVEL_LINK) $(LDFLAGS) -shared -fPIC -o $(PATH_LIB_DIVERSUTILS)/libudpipe.so udpipe/src_lib_only/udpipe.cpp $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	export LD_LIBRARY_PATH=$(PATH_LIB_DIVERSUTILS):$$LD_LIBRARY_PATH

obj_udpipe_interface: ensure_bin_dir
	set -eu
	$(CCPP) $(CPP_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_LINK) $(CPP_MACROS) -o build/udpipe_interface.o -c src/target/udpipe_interface.cpp # ok

obj_measurement: ensure_bin_dir
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_LINK) $(CPP_MACROS) -o build/main_measurement.o -c src/target/main_measurement.c # ok

measurement_link: obj_udpipe_interface cpu dfunctions distributions distances stats graph measurement jsonl_parser jsonl_load cupt_parser cupt_load main_measurement sanitize logging sorted_array udpipe_interface_conversion
	set -eu
	$(CCPP) $(CPP_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_LINK) $(CPP_MACROS) $(LDFLAGS) -flto build/udpipe_interface.o build/cpu.o build/dfunctions.o build/distributions.o build/distances.o build/stats.o build/graph.o build/measurement.o build/jsonl_parser.o build/jsonl_load.o build/cupt_parser.o build/cupt_load.o build/main_measurement.o build/sanitize.o build/logging.o build/sorted_array.o build/udpipe_interface_conversion.o -o bin/main_measurement $(LINKER_FLAGS_EXTRA) $(LINKER_FLAGS) # ok?

measurement_full: ensure_prepared_udpipe measurement_link

measurement_standard:
	$(CC) --version
	mkdir -p bin
	$(CC) $(C_VERSION) $(CFLAGS_STANDARD) $(CPPFLAGS) $(OPT_LEVEL) $(CPP_MACROS) -o bin/main_measurement src/main_measurement.c $(LINKER_FLAGS)

measurement_optim:
	$(CC) --version
	mkdir -p bin
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -o bin/main_measurement src/main_measurement.c $(LINKER_FLAGS)

measurement_debug:
	$(CC) --version
	mkdir -p bin
	$(CC) $(C_VERSION) $(CFLAGS_DEBUG) $(CPPFLAGS) $(OPT_LEVEL_DEBUG) $(CPP_MACROS) -o bin/main_measurement src/main_measurement.c $(LINKER_FLAGS)

measurement_pedantic:
	$(CC) --version
	mkdir -p bin
	$(CC) $(C_VERSION) $(CFLAGS_PEDANTIC) $(CPPFLAGS) $(OPT_LEVEL) $(CPP_MACROS) -o bin/main_measurement src/main_measurement.c $(LINKER_FLAGS)

measurement_coverage:
	$(CC) --version
	mkdir -p bin
	$(CC) $(C_VERSION) $(CFLAGS_COVERAGE) $(CPPFLAGS) $(OPT_LEVEL) $(CPP_MACROS)	-o bin/main_measurement src/main_measurement.c $(LINKER_FLAGS)

measurement_profiling:
	$(CC) --version
	mkdir -p bin
	$(CC) $(C_VERSION) $(CFLAGS_PROFILING) $(CPPFLAGS) $(OPT_LEVEL) $(CPP_MACROS) -o bin/main_measurement src/main_measurement.c $(LINKER_FLAGS)

measurement_python_interface:
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_STANDARD) $(CPPFLAGS) $(OPT_LEVEL) $(PYTHON_INCLUDE_FLAG) -shared -fPIC -o diversutils.so src/diversutilsmodule.c $(LINKER_FLAGS)

# ======== multi-step build ========

mkdir_bin:
	set -eu
	mkdir -p bin

mkdir_build:
	set -eu
	mkdir -p build

cpu: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/cpu.o src/target/cpu.c

dfunctions: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/dfunctions.o src/target/dfunctions.c

distances: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/distances.o src/target/distances.c

distributions: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/distributions.o src/target/distributions.c

cupt_parser: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/cupt_parser.o src/target/cupt/parser.c

cupt_load: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/cupt_load.o src/target/cupt/load.c

jsonl_parser: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/jsonl_parser.o src/target/jsonl/parser.c

jsonl_load: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/jsonl_load.o src/target/jsonl/load.c

sanitize: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/sanitize.o src/target/sanitize.c

graph: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/graph.o src/target/graph.c

measurement: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/measurement.o src/target/measurement.c

logging: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/logging.o src/target/logging.c

stats: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/stats.o src/target/stats.c

sorted_array: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/sorted_array.o src/target/sorted_array/array.c

udpipe_interface_conversion: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/udpipe_interface_conversion.o src/target/udpipe_interface/conversion.c

main_measurement: mkdir_build
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -c -o build/main_measurement.o src/target/main_measurement.c

all_multistep: mkdir_bin cpu dfunctions distances distributions cupt_parser cupt_load jsonl_parser jsonl_load sanitize graph measurement logging stats sorted_array obj_udpipe_interface
	set -eu
	$(CC) $(C_VERSION) $(CFLAGS_OPTIM) $(CPPFLAGS) $(OPT_LEVEL_OPTIM) $(CPP_MACROS) -flto -o bin/main_measurement src/target/main_measurement.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/cupt_load.o build/jsonl_parser.o build/jsonl_load.o build/sanitize.o build/graph.o build/measurement.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS)



# ======== converting jsonl to word2vec format ========

setup_word2vec_learner:
	$(CC) --version
	git --version
	git clone https://github.com/dav/word2vec
	cd word2vec
	make build

jsonl_to_word2vec_format:
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_STANDARD) $(CPPFLAGS) $(OPT_LEVEL) $(CPP_MACROS_JSONL) -o bin/jsonl_to_word2vec_format src/main_jsonl_to_word2vec_format.c $(LINKER_FLAGS)

# ======== test ========

test_graph_relative_proportion: ensure_prepared_udpipe share_udpipe cpu dfunctions distances distributions cupt_parser jsonl_parser sanitize graph logging stats sorted_array obj_udpipe_interface
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_TEST) $(CPPFLAGS_TEST) $(OPT_LEVEL_TEST) $(LDFLAGS) -DENABLE_AVX256=0 -DTEST_GRAPH_RELATIVE_PROPORTION -o test/test_graph_relative_proportion test/main_test.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/jsonl_parser.o build/sanitize.o build/graph.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	./test/test_graph_relative_proportion

test_entropy_shannon_weaver: ensure_prepared_udpipe share_udpipe cpu dfunctions distances distributions cupt_parser jsonl_parser sanitize graph logging stats sorted_array obj_udpipe_interface
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_TEST) $(CPPFLAGS_TEST) $(OPT_LEVEL_TEST) $(LDFLAGS) -DENABLE_AVX256=0 -DTEST_ENTROPY_SHANNON_WEAVER -o test/test_entropy_shannon_weaver test/main_test.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/jsonl_parser.o build/sanitize.o build/graph.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	./test/test_entropy_shannon_weaver

test_entropy_renyi: ensure_prepared_udpipe share_udpipe cpu dfunctions distances distributions cupt_parser jsonl_parser sanitize graph logging stats sorted_array obj_udpipe_interface
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_TEST) $(CPPFLAGS_TEST) $(OPT_LEVEL_TEST) $(LDFLAGS) -DENABLE_AVX256=0 -DTEST_ENTROPY_RENYI -o test/test_entropy_renyi test/main_test.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/jsonl_parser.o build/sanitize.o build/graph.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	./test/test_entropy_renyi

test_entropy_patil_taillie: ensure_prepared_udpipe share_udpipe cpu dfunctions distances distributions cupt_parser jsonl_parser sanitize graph logging stats sorted_array obj_udpipe_interface
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_TEST) $(CPPFLAGS_TEST) $(OPT_LEVEL_TEST) $(LDFLAGS) -DENABLE_AVX256=0 -DTEST_ENTROPY_PATIL_TAILLIE -o test/test_entropy_patil_taillie test/main_test.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/jsonl_parser.o build/sanitize.o build/graph.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	./test/test_entropy_patil_taillie

test_entropy_q_logarithmic: ensure_prepared_udpipe share_udpipe cpu dfunctions distances distributions cupt_parser jsonl_parser sanitize graph logging stats sorted_array obj_udpipe_interface
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_TEST) $(CPPFLAGS_TEST) $(OPT_LEVEL_TEST) $(LDFLAGS) -DENABLE_AVX256=0 -DTEST_ENTROPY_Q_LOGARITHMIC -o test/test_entropy_q_logarithmic test/main_test.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/jsonl_parser.o build/sanitize.o build/graph.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	./test/test_entropy_q_logarithmic

test_equivalence_entropy: ensure_prepared_udpipe share_udpipe cpu dfunctions distances distributions cupt_parser jsonl_parser sanitize graph logging stats sorted_array obj_udpipe_interface
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_TEST) $(CPPFLAGS_TEST) $(OPT_LEVEL_TEST) $(LDFLAGS) -DENABLE_AVX256=0 -DTEST_EQUIVALENCE_ENTROPY -o test/test_equivalence_entropy test/main_test.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/jsonl_parser.o build/sanitize.o build/graph.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	./test/test_equivalence_entropy

test_all: ensure_prepared_udpipe share_udpipe cpu dfunctions distances distributions cupt_parser jsonl_parser sanitize graph logging stats sorted_array obj_udpipe_interface
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS_TEST) $(CPPFLAGS_TEST) $(OPT_LEVEL_TEST) $(LDFLAGS) -DENABLE_AVX256=0 -DTEST_ALL -o test/test_all test/main_test.c build/cpu.o build/dfunctions.o build/distances.o build/distributions.o build/cupt_parser.o build/jsonl_parser.o build/sanitize.o build/graph.o build/logging.o build/stats.o build/sorted_array.o build/udpipe_interface.o $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	./test/test_all

# ======== formatting ========

format_include:
	clang-format -i src/include/*.h src/include/*/*.h
	clang-format -i src/include/*.hpp

format_target:
	clang-format -i src/target/*.c src/target/*/*.c
	clang-format -i src/target/*.cpp

format_all: format_include format_target
