# ------------------------------
# DIVERSITY MEASUREMENT MAKEFILE
# ------------------------------
#
# To compile, run:
# 	make measurement_full
#

OPT_LEVEL_DEBUG = -Og
OPT_LEVEL_STANDARD = -O2
OPT_LEVEL_OPTIM = -O3
OPT_LEVEL_TEST = -O2
OPT_LEVEL_COMPACT_STANDARD = -Os
OPT_LEVEL_COMPACT_OPTIM = -Oz

ifneq ($(VERBOSE), 1)
    .SILENT:
endif

ifeq ($(COMPACT), 1)
	OPT_LEVEL = $(OPT_LEVEL_COMPACT_STANDARD)
else
	OPT_LEVEL = $(OPT_LEVEL_STANDARD)
endif

ifeq ($(PEDANTIC), 1)
	PEDANTIC_FLAG = -pedantic -Werror
endif

ifeq ($(NATIVE), 1)
	NATIVE_FLAG = -march=native
    ifeq ($(AVX), 0)
        ENABLE_AVX256 = 0
        ENABLE_AVX512 = 0
    else
    	ifeq ($(grep avx2 /proc/cpuinfo), "")
    	    ENABLE_AVX256 = 0
    	else
    	    ENABLE_AVX256 = 1
    	endif
    	
    	ifeq ($(grep avx512 /proc/cpuinfo), "")
    	    ENABLE_AVX512 = 0
    	else
    	    ENABLE_AVX512 = 1
    	endif
    endif
else
    ENABLE_AVX256 = 0
    ENABLE_AVX512 = 0
endif

ifeq ($(OPTIM), 1)
	FAST_MATH_FLAG = -ffast-math
	ifeq ($(COMPACT), 1)
		OPT_LEVEL = $(OPT_LEVEL_COMPACT_OPTIM)
	else
		OPT_LEVEL = $(OPT_LEVEL_OPTIM)
	endif
endif

ifeq ($(DEBUG), 1)
	DEBUG_FLAG = -DDEBUG -g3
	OPT_LEVEL = $(OPT_LEVEL_DEBUG)
else
	DEBUG_FLAG = -DNDEBUG
endif


ifeq ($(COMPILER), GNU)
	CC = gcc
	CCPP = g++
	DEBUGGER = gdb
else
	CC = clang
	CCPP = clang++
	DEBUGGER = lldb
endif

# JSONL only
#TOKENIZATION_METHOD = 0 # REGEX
#TOKENIZATION_METHOD = 1 # UDPipe
TOKENIZATION_METHOD = 2 # UDPipe embedded

ifeq ($(PYTHON_BUILD), 1)
	TOKENIZATION_METHOD = 0 # REGEX
endif

# -----------------
# COMPILING OPTIONS
# -----------------

C_VERSION = -std=c99
CXX_VERSION = -std=c++11
CFLAGS = -Wall -Wextra -Wformat -Wformat-security $(PEDANTIC_FLAG) $(NATIVE_FLAG) $(FAST_MATH_FLAG) $(DEBUG_FLAG)
CPPFLAGS = -Isrc/include -Iudpipe/src_lib_only -I/usr/include
CPPFLAGS_CXX = -I/usr/include/x86_64-linux-gnu/c++/11 -I/usr/include/c++/11
CPPFLAGS_TEST = $(CPPFLAGS) -Itest/include
PATH_LIB_DIVERSUTILS = $(HOME)/.local/lib/diversutils
LDFLAGS_MINIMUM = -L$(PATH_LIB_DIVERSUTILS) -L/usr/lib/x84_64-linux-gnu
ifeq ($(origin CONDA_PREFIX), undefined)
	LDFLAGS = $(LDFLAGS_MINIMUM)
else
	LDFLAGS = -L$(CONDA_PREFIX)/lib $(LDFLAGS_MINIMUM)
endif
PYTHON_INCLUDE_FLAG = -I$(HOME)/.local/conda/include/python3.11
LINKER_FLAGS = -lc -lm -lrt -pthread
LINKER_FLAG_CXX = -lstdc++
LINKER_FLAG_UDPIPE = -ludpipe
LINKER_FLAG_DIVERSUTILS = -ldiversutils
ifeq ($(TOKENIZATION_METHOD), 0)
	LINKER_FLAGS_EXTRA = $(LINKER_FLAG_CXX)
else
	LINKER_FLAGS_EXTRA = $(LINKER_FLAG_CXX) $(LINKER_FLAG_UDPIPE)
endif

# ----------------------
# UPDATE LD_LIBRARY_PATH
# ----------------------

export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):$(PATH_LIB_DIVERSUTILS)

update_bashrc:
	set -eu ; if [ "$$(grep "$# >>> diversutils initialise >>>" ~/.bashrc)" = "" ] ; then echo "# >>> diversutils initialise >>>" >> $(HOME)/.bashrc ; echo "export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:\$$HOME/.local/lib/diversutils" >> $(HOME)/.bashrc ; echo "# <<< diversutils initialise <<<" >> $(HOME)/.bashrc ; echo "INFO: \"\033[1m\033[32m$(HOME)/.bashrc\"\033[0m has been updated." ; fi
	set -eu ; echo "INFO: \033[1mConsider running \"\033[33msource ~/.bashrc\033[0m\033[1m\" if you cannot link against *.so files.\033[0m"
	
# ----------------
# MACRO PARAMETERS
# ----------------

NUM_ROW_THREADS = -1
NUM_MATRIX_THREADS = -1
ROW_GENERATION_BATCH_SIZE = -1

NUM_FILE_READING_THREADS = 4

ENABLE_NON_DISPARITY_MULTITHREADING = 1

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

ENABLE_SW_E_PRIME_CAMARGO1993_MULTITHREADING = 1

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

ENABLE_TOKEN_UTF8_NORMALISATION = 1

CPP_MACRO_MULTITHREADING = -DNUM_ROW_THREADS=$(NUM_ROW_THREADS) -DNUM_MATRIX_THREADS=$(NUM_MATRIX_THREADS) -DROW_GENERATION_BATCH_SIZE=$(ROW_GENERATION_BATCH_SIZE) -DNUM_FILE_READING_THREADS=$(NUM_FILE_READING_THREADS) -DENABLE_ITERATIVE_DISTANCE_COMPUTATION=$(ENABLE_ITERATIVE_DISTANCE_COMPUTATION) -DENABLE_MULTITHREADED_ROW_GENERATION=$(ENABLE_MULTITHREADED_ROW_GENERATION) -DENABLE_MULTITHREADED_MATRIX_GENERATION=$(ENABLE_MULTITHREADED_MATRIX_GENERATION) -DENABLE_NON_DISPARITY_MULTITHREADING=$(ENABLE_NON_DISPARITY_MULTITHREADING) -DENABLE_SW_E_PRIME_CAMARGO1993_MULTITHREADING=$(ENABLE_SW_E_PRIME_CAMARGO1993_MULTITHREADING)

CPP_MACRO_AVX = -DENABLE_AVX256=$(ENABLE_AVX256) -DENABLE_AVX512=$(ENABLE_AVX512)

CPP_MACRO_DISPARITY = -DENABLE_DISPARITY_FUNCTIONS=$(ENABLE_DISPARITY_FUNCTIONS) -DENABLE_STIRLING=$(ENABLE_STIRLING) -DENABLE_RICOTTA_SZEIDL=$(ENABLE_RICOTTA_SZEIDL) -DENABLE_PAIRWISE=$(ENABLE_PAIRWISE) -DENABLE_LEXICOGRAPHIC=$(ENABLE_LEXICOGRAPHIC) -DENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY=$(ENABLE_CHAO_ET_AL_FUNCTIONAL_DIVERSITY) -DENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY=$(ENABLE_SCHEINER_SPECIES_PHYLOGENETIC_FUNCTIONAL_DIVERSITY) -DENABLE_LEINSTER_COBBOLD_DIVERSITY=$(ENABLE_LEINSTER_COBBOLD_DIVERSITY) -DENABLE_FUNCTIONAL_EVENNESS=$(ENABLE_FUNCTIONAL_EVENNESS) -DENABLE_FUNCTIONAL_DISPERSION=$(ENABLE_FUNCTIONAL_DISPERSION) -DENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED=$(ENABLE_FUNCTIONAL_DIVERGENCE_MODIFIED)

CPP_MACRO_NON_DISPARITY = -DENABLE_NON_DISPARITY_FUNCTIONS=$(ENABLE_NON_DISPARITY_FUNCTIONS) -DENABLE_SHANNON_WEAVER_ENTROPY=$(ENABLE_SHANNON_WEAVER_ENTROPY) -DENABLE_GOOD_ENTROPY=$(ENABLE_GOOD_ENTROPY) -DENABLE_RENYI_ENTROPY=$(ENABLE_RENYI_ENTROPY) -DENABLE_PATIL_TAILLIE_ENTROPY=$(ENABLE_PATIL_TAILLIE_ENTROPY) -DENABLE_Q_LOGARITHMIC_ENTROPY=$(ENABLE_Q_LOGARITHMIC_ENTROPY) -DENABLE_SIMPSON_INDEX=$(ENABLE_SIMPSON_INDEX) -DENABLE_SIMPSON_DOMINANCE_INDEX=$(ENABLE_SIMPSON_DOMINANCE_INDEX) -DENABLE_HILL_NUMBER_STANDARD=$(ENABLE_HILL_NUMBER_STANDARD) -DENABLE_HILL_EVENNESS=$(ENABLE_HILL_EVENNESS) -DENABLE_BERGER_PARKER_INDEX=$(ENABLE_BERGER_PARKER_INDEX) -DENABLE_JUNGE1994_PAGE20=$(ENABLE_JUNGE1994_PAGE20) -DENABLE_JUNGE1994_PAGE22=$(ENABLE_JUNGE1994_PAGE22) -DENABLE_BRILLOUIN_DIVERSITY=$(ENABLE_BRILLOUIN_DIVERSITY) -DENABLE_MCINTOSH_INDEX=$(ENABLE_MCINTOSH_INDEX) -DENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975=$(ENABLE_SW_ENTROPY_OVER_LOG_N_SPECIES_PIELOU1975) -DENABLE_SW_E_HEIP=$(ENABLE_SW_E_HEIP) -DENABLE_SW_E_ONE_MINUS_D=$(ENABLE_SW_E_ONE_MINUS_D) -DENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964=$(ENABLE_SW_E_ONE_OVER_LN_D_WILLIAMS1964) -DENABLE_SW_E_MINUS_LN_D_PIELOU1977=$(ENABLE_SW_E_MINUS_LN_D_PIELOU1977) -DENABLE_SW_F_2_1_ALATALO1981=$(ENABLE_SW_F_2_1_ALATALO1981) -DENABLE_SW_G_2_1_MOLINARI1989=$(ENABLE_SW_G_2_1_MOLINARI1989) -DENABLE_SW_E_BULLA1994=$(ENABLE_SW_E_BULLA1994) -DENABLE_SW_O_BULLA1994=$(ENABLE_SW_O_BULLA1994) -DENABLE_SW_E_MCI_PIELOU1969=$(ENABLE_SW_E_MCI_PIELOU1969) -DENABLE_SW_E_PRIME_CAMARGO1993=$(ENABLE_SW_E_PRIME_CAMARGO1993) -DENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL=$(ENABLE_SW_E_VAR_SMITH_AND_WILSON1996_ORIGINAL)

CPP_MACRO_TIMING = -DENABLE_TIMINGS=$(ENABLE_TIMINGS)

CPP_MACRO_RECOMPUTE = -DENABLE_SENTENCE_COUNT_RECOMPUTE_STEP=$(ENABLE_SENTENCE_COUNT_RECOMPUTE_STEP) -DSENTENCE_COUNT_RECOMPUTE_STEP=$(SENTENCE_COUNT_RECOMPUTE_STEP) -DSENTENCE_RECOMPUTE_STEP_USE_LOG10=$(SENTENCE_RECOMPUTE_STEP_USE_LOG10) -DSENTENCE_COUNT_RECOMPUTE_STEP_LOG10=$(SENTENCE_COUNT_RECOMPUTE_STEP_LOG10) -DENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP=$(ENABLE_DOCUMENT_COUNT_RECOMPUTE_STEP) -DDOCUMENT_COUNT_RECOMPUTE_STEP=$(DOCUMENT_COUNT_RECOMPUTE_STEP) -DDOCUMENT_RECOMPUTE_STEP_USE_LOG10=$(DOCUMENT_RECOMPUTE_STEP_USE_LOG10) -DDOCUMENT_COUNT_RECOMPUTE_STEP_LOG10=$(DOCUMENT_COUNT_RECOMPUTE_STEP_LOG10)

CPP_MACRO_IO = -DINPUT_PATH=$(INPUT_PATH) -DOUTPUT_PATH=$(OUTPUT_PATH) -DOUTPUT_PATH_TIMING=$(OUTPUT_PATH_TIMING) -DOUTPUT_PATH_MEMORY=$(OUTPUT_PATH_MEMORY) -DW2V_PATH=$(W2V_PATH)

CPP_MACRO_OTHER = -DTARGET_COLUMN=$(TARGET_COLUMN) -DTOKENIZATION_METHOD=$(TOKENIZATION_METHOD) -DENABLE_TOKEN_UTF8_NORMALISATION=$(ENABLE_TOKEN_UTF8_NORMALISATION) -DCOMPILATION_DIR=\"$(PWD)\"

CPP_MACROS = $(CPP_MACRO_MULTITHREADING) $(CPP_MACRO_AVX) $(CPP_MACRO_DISPARITY) $(CPP_MACRO_NON_DISPARITY) $(CPP_MACRO_TIMING) $(CPP_MACRO_RECOMPUTE) $(CPP_MACRO_IO) $(CPP_MACRO_OTHER)

# JSONL to W2V format variables

TARGET_NUMBER_DOCUMENTS = 1000000

CPP_MACROS_JSONL = -DTARGET_NUMBER_DOCUMENTS=$(TARGET_NUMBER_DOCUMENTS)

################

udpipe:
	echo "INFO: Ensuring UDPipe is updated"
	git submodule update --init --recursive udpipe

# ======== multi-step build ========

TGT = src/target
INC = src/include
BLD = build
BIN = bin
TST = test

udpipe/src_lib_only/udpipe.cpp: udpipe
	echo "INFO: Ensuring UDPipe is prepared"
	set -eu ; make -s --directory=udpipe/src
	set -eu ; make -s --directory=udpipe/src_lib_only

$(PATH_LIB_DIVERSUTILS)/libudpipe.so: udpipe/src_lib_only/udpipe.cpp $(TGT)/udpipe/interface.cpp
	echo "INFO: Ensuring directory \"\033[1m\033[32m$(PATH_LIB_DIVERSUTILS)\033[0m\" exists"
	mkdir -p $(PATH_LIB_DIVERSUTILS)
	echo "INFO: Generating shared library \"\033[1m\033[32m$@\033[0m\""
	$(CCPP) $(CPP_VERSION) $(CFLAGS) $(CPPFLAGS) $(CPPFLAGS_CXX) $(OPT_LEVEL) $(CPP_MACROS) $(LDFLAGS) -shared -fPIC -o $@ $^ $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	export LD_LIBRARY_PATH=$(PATH_LIB_DIVERSUTILS):$$LD_LIBRARY_PATH

BLD_DIRECTORIES = $(BLD) $(BLD)/cupt $(BLD)/jsonl $(BLD)/sorted_array $(BLD)/unicode $(BLD)/cfgparser $(BLD)/udpipe/interface

$(BLD_DIRECTORIES): %:
	set -eu ; echo "INFO: Creating \"\033[1m\033[32m$@\033[0m\""
	set -eu ; mkdir -p $@
    
$(BLD)/.placeholder: $(BLD_DIRECTORIES)
	set -eu ; if ! [ -f $@ ] ; then touch $@ ; fi

$(BIN)/.placeholder:
	set -eu ; mkdir -p $(BIN) ; touch $@

$(INC)/dfunctions.h: $(INC)/graph.h $(INC)/general_constants.h
$(INC)/distances.h: $(INC)/graph.h
$(INC)/distributions.h: $(INC)/graph.h
$(INC)/measurement.h: $(INC)/graph.h $(INC)/sorted_array/array.h
$(INC)/cupt/parser.h: $(INC)/cupt/constants.h
$(INC)/cupt/load.h: $(INC)/graph.h $(INC)/sorted_array/array.h $(INC)/measurement.h
$(INC)/jsonl/parser.h: $(INC)/jsonl/constants.h
$(INC)/jsonl/load.h: $(INC)/graph.h $(INC)/sorted_array/array.h $(INC)/measurement.h
$(INC)/sorted_array/array.h: $(INC)/sorted_array/constants.h
$(INC)/unicode/utf8.h: $(INC)/unicode/unicode.h


$(TGT)/cpu.c: $(INC)/cpu.h
$(TGT)/graph.c: $(INC)/general_constants.h $(INC)/distances.h $(INC)/cupt/parser.h $(INC)/stats.h $(INC)/logging.h
$(TGT)/dfunctions.c: $(INC)/dfunctions.h $(INC)/graph.h $(INC)/general_constants.h
$(TGT)/distances.c: $(INC)/distances.h
$(TGT)/distributions.c: $(INC)/distributions.h
$(TGT)/stats.c: $(INC)/stats.h
$(TGT)/logging.c.c: $(INC)/logging.h
$(TGT)/measurement.c.c: $(INC)/measurement.h $(INC)/dfunctions.h $(INC)/distributions.h $(INC)/graph.h $(INC)/cpu.h $(INC)/sorted_array/array.h $(INC)/logging.h $(INC)/stats.h
$(TGT)/sanitize.c: $(INC)/sanitize.h
$(TGT)/cupt/parser.c: $(INC)/cupt/parser.h $(INC)/cupt/constants.h
$(TGT)/cupt/load.c: $(INC)/cupt/parser.h $(INC)/cupt/constants.h $(INC)/cupt/load.h
$(TGT)/jsonl/parser.c: $(INC)/jsonl/parser.h $(INC)/jsonl/constants.h
$(TGT)/jsonl/load.c: $(INC)/jsonl/parser.h $(INC)/jsonl/constants.h $(INC)/jsonl/load.h
$(TGT)/cfgparser/parser.c: $(INC)/cfgparser/parser.h
$(TGT)/unicode/utf8.c: $(INC)/unicode/unicode.h $(INC)/unicode/utf8.h
$(TGT)/sorted_array/array.c: $(INC)/sorted_array/array.h $(INC)/sorted_array/constants.h
$(TGT)/udpipe/interface.cpp: $(INC)/udpipe_interface.hpp
$(TGT)/udpipe/interface/conversion.c: $(INC)/udpipe/interface/cinterface.h $(INC)/udpipe/interface/conversion.h
$(TGT)/main_measurement.c: $(INC)/cpu.h $(INC)/graph.h $(INC)/distributions.h $(INC)/stats.h $(INC)/dfunctions.h $(INC)/sorted_array/array.h $(INC)/logging.h $(INC)/measurement.h $(INC)/jsonl/parser.h $(INC)/jsonl/load.h $(INC)/cupt/parser.h $(INC)/cupt/load.h $(INC)/udpipe/interface/cinterface.h $(INC)/macroconfig.h

DIVERSUTILS_CFILES = $(TGT)/cpu.c $(TGT)/graph.c $(TGT)/dfunctions.c $(TGT)/distances.c $(TGT)/distributions.c $(TGT)/stats.c $(TGT)/logging.c $(TGT)/measurement.c $(TGT)/sanitize.c $(TGT)/cupt/parser.c $(TGT)/cupt/load.c $(TGT)/jsonl/parser.c $(TGT)/jsonl/load.c $(TGT)/sorted_array/array.c $(TGT)/unicode/utf8.c $(TGT)/cfgparser/parser.c $(TGT)/udpipe/interface/conversion.c
DIVERSUTILS_C_OBJECTS = $(BLD)/cpu.o $(BLD)/graph.o $(BLD)/dfunctions.o $(BLD)/distances.o $(BLD)/distributions.o $(BLD)/stats.o $(BLD)/logging.o $(BLD)/measurement.o $(BLD)/sanitize.o $(BLD)/cupt/parser.o $(BLD)/cupt/load.o $(BLD)/jsonl/parser.o $(BLD)/jsonl/load.o $(BLD)/sorted_array/array.o $(BLD)/unicode/utf8.o $(BLD)/cfgparser/parser.o $(BLD)/udpipe/interface/conversion.o
DIVERSUTILS_CXX_OBJECTS = $(BLD)/udpipe/interface.o

$(DIVERSUTILS_C_OBJECTS) $(BLD)/main_measurement.o: $(BLD)/%.o: $(TGT)/%.c $(BLD)/.placeholder
	echo "INFO: Building \"\033[1m\033[32m$@\033[0m\" from \"\033[1m\033[34m$<\033[0m\""
	$(CC) $< $(C_VERSION) $(CFLAGS) $(CPPFLAGS) $(OPT_LEVEL) $(CPP_MACROS) -c -o $@
$(DIVERSUTILS_CXX_OBJECTS): $(BLD)/%.o: $(TGT)/%.cpp $(BLD)/.placeholder
	echo "INFO: Building \"\033[1m\033[32m$@\033[0m\" from \"\033[1m\033[34m$<\033[0m\""
	$(CCPP) $< $(CXX_VERSION) $(CFLAGS) $(CPPFLAGS) $(CPPFLAGS_CXX) $(OPT_LEVEL) $(CPP_MACROS) -c -o $@

$(BIN)/main_measurement: $(BIN)/.placeholder $(DIVERSUTILS_C_OBJECTS) $(DIVERSUTILS_CXX_OBJECTS) $(BLD)/main_measurement.o update_bashrc
	echo "INFO: Linking \"\033[1m\033[31m$@\033[0m\" from \"\033[1m\033[32m$(DIVERSUTILS_C_OBJECTS) $(DIVERSUTILS_CXX_OBJECTS)\033[0m\""
	$(CCPP) $(DIVERSUTILS_C_OBJECTS) $(DIVERSUTILS_CXX_OBJECTS) $(BLD)/main_measurement.o $(CXX_VERSION) $(CFLAGS) $(CPPFLAGS) $(CPPFLAGS_CXX) $(OPT_LEVEL) $(CPP_MACROS) $(LDFLAGS) -flto -o $@ $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	echo "INFO: You may now run \"\033[1m\033[33mcd bin ; ./main_measurement\033[0m\""

measurement: $(BLD)/.placeholder $(BIN)/.placeholder $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(BIN)/main_measurement

# for backwards compatibility

measurement_standard: measurement
measurement_link: measurement
measurement_full: measurement

# for Python interface

$(PATH_LIB_DIVERSUTILS)/libdiversutils.so: $(DIVERSUTILS_CFILES)
	echo "INFO: Ensuring directory \"\033[1m\033[32m$(PATH_LIB_DIVERSUTILS)\033[0m\" exists"
	mkdir -p $(PATH_LIB_DIVERSUTILS)
	echo "INFO: Generating shared library \"\033[1m\033[32m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS) $(OPT_LEVEL) $(LDFLAGS) -shared -fPIC -o $@ $^ $(LINKER_FLAGS) $(LINKER_FLAG_UDPIPE)

shared_libraries: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(PATH_LIB_DIVERSUTILS)/libdiversutils.so

all: shared_libraries measurement

# ======== converting jsonl to word2vec format ========

setup_word2vec_learner:
	$(CC) --version
	git --version
	git clone https://github.com/dav/word2vec
	cd word2vec
	make build

jsonl_to_word2vec_format:
	$(CC) --version
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS) $(OPT_LEVEL) $(CPP_MACROS_JSONL) -o bin/jsonl_to_word2vec_format src/main_jsonl_to_word2vec_format.c $(LINKER_FLAGS)

# ======== test ========

$(TST)/include/test_general.h: $(INC)/logging.h
$(TST)/include/test_entropy.h: $(TST)/include/test_general.h $(INC)/graph.h $(INC)/dfunctions.h $(INC)/distances.h
$(TST)/include/test_equivalence.h: $(TST)/include/test_general.h $(INC)/graph.h $(INC)/dfunctions.h $(INC)/distances.h
$(TST)/include/test_graph.h: $(TST)/include/test_general.h $(INC)/graph.h $(INC)/dfunctions.h $(INC)/distances.h

$(TST)/main_test.c: $(TST)/include/test_general.h $(TST)/include/test_graph.h $(TST)/include/test_entropy.h $(TST)/include/test_equivalence.h

$(TST)/test_all: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(DIVERSUTILS_C_OBJECTS) $(TST)/include/test_general.h $(TST)/include/test_graph.h $(TST)/include/test_entropy.h $(TST)/include/test_equivalence.h $(TST)/main_test.c
	echo "INFO: Generating \"\033[1m\033[31m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS_TEST) $(OPT_LEVEL) $(LDFLAGS) -DENABLE_AVX256=0 -DENABLE_AVX512=0 -DTEST_ALL -o test/test_all test/main_test.c $(DIVERSUTILS_C_OBJECTS) $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
$(TST)/test_graph_relative_proportion: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(DIVERSUTILS_C_OBJECTS) $(TST)/include/test_general.h $(TST)/include/test_graph.h $(TST)/main_test.c
	echo "INFO: Generating \"\033[1m\033[31m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS_TEST) $(OPT_LEVEL) $(LDFLAGS) -DENABLE_AVX256=0 -DENABLE_AVX512=0 -DTEST_GRAPH_RELATIVE_PROPORTION -o test/test_graph_relative_proportion test/main_test.c $(DIVERSUTILS_C_OBJECTS) $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
$(TST)/test_entropy_shannon_weaver: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(DIVERSUTILS_C_OBJECTS) $(TST)/include/test_general.h $(TST)/include/test_entropy.h $(TST)/main_test.c
	echo "INFO: Generating \"\033[1m\033[31m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS_TEST) $(OPT_LEVEL) $(LDFLAGS) -DENABLE_AVX256=0 -DENABLE_AVX512=0 -DTEST_ENTROPY_SHANNON_WEAVER -o test/test_entropy_shannon_weaver test/main_test.c $(DIVERSUTILS_C_OBJECTS) $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
$(TST)/test_entropy_patil_taillie: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(DIVERSUTILS_C_OBJECTS) $(TST)/include/test_general.h $(TST)/include/test_entropy.h $(TST)/main_test.c
	echo "INFO: Generating \"\033[1m\033[31m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS_TEST) $(OPT_LEVEL) $(LDFLAGS) -DENABLE_AVX256=0 -DENABLE_AVX512=0 -DTEST_ENTROPY_PATIL_TAILLIE -o test/test_entropy_patil_taillie test/main_test.c $(DIVERSUTILS_C_OBJECTS) $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
$(TST)/test_entropy_renyi: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(DIVERSUTILS_C_OBJECTS) $(TST)/include/test_general.h $(TST)/include/test_entropy.h $(TST)/main_test.c
	echo "INFO: Generating \"\033[1m\033[31m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS_TEST) $(OPT_LEVEL) $(LDFLAGS) -DENABLE_AVX256=0 -DENABLE_AVX512=0 -DTEST_ENTROPY_RENYI -o test/test_entropy_renyi test/main_test.c $(DIVERSUTILS_C_OBJECTS) $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
$(TST)/test_entropy_q_logarithmic: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(DIVERSUTILS_C_OBJECTS) $(TST)/include/test_general.h $(TST)/include/test_entropy.h $(TST)/main_test.c
	echo "INFO: Generating \"\033[1m\033[31m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS_TEST) $(OPT_LEVEL) $(LDFLAGS) -DENABLE_AVX256=0 -DENABLE_AVX512=0 -DTEST_ENTROPY_Q_LOGARITHMIC -o test/test_entropy_q_logarithmic test/main_test.c $(DIVERSUTILS_C_OBJECTS) $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
$(TST)/test_equivalence_entropy: $(PATH_LIB_DIVERSUTILS)/libudpipe.so $(DIVERSUTILS_C_OBJECTS) $(TST)/include/test_general.h $(TST)/include/test_equivalence.h $(TST)/main_test.c
	echo "INFO: Generating \"\033[1m\033[31m$@\033[0m\""
	$(CC) $(C_VERSION) $(CFLAGS) $(CPPFLAGS_TEST) $(OPT_LEVEL) $(LDFLAGS) -DENABLE_AVX256=0 -DENABLE_AVX512=0 -DTEST_EQUIVALENCE_ENTROPY -o test/test_equivalence_entropy test/main_test.c $(DIVERSUTILS_C_OBJECTS) $(LINKER_FLAGS) $(LINKER_FLAGS_EXTRA)
	

# ======== format ========

format_include:
	clang-format -i src/include/*.h src/include/*/*.h
	clang-format -i src/include/*.hpp

format_target:
	clang-format -i src/target/*.c src/target/*/*.c
	clang-format -i src/target/*/*.cpp

format_all: format_include format_target

# ======== docs ========

sphinx:
	sphinx-build -b html docs build/html

