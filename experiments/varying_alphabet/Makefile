# --- Base directory (make this flexible) ---
BASEDIR ?= ../..

# --- Binaries ---
BINARY1 = federated_local
BINARY2 = federated_central
BINARY3 = krr_clones_feldman2020
# BINARY4 = generic_clones_feldman2020
# BINARY5 = generic_clones_feldman2023
# BINARY6 = generic_clones_ours2022
BINARIES = $(BINARY1) $(BINARY2) $(BINARY3) $(BINARY4) $(BINARY5) $(BINARY6)

# --- Program-specific sources ---
PROG1_MAIN = $(BASEDIR)/srcprog/stat/multiple_mechanisms/local/federated_local.cpp
PROG2_MAIN = $(BASEDIR)/srcprog/stat/multiple_mechanisms/central/federated_central.cpp
PROG3_MAIN = $(BASEDIR)/srcprog/stat/multiple_mechanisms/shuffling/krr_clones_feldman2020.cpp
# PROG4_MAIN = $(BASEDIR)/srcprog/stat/multiple_mechanisms/shuffling/generic_clones_feldman2020.cpp
# PROG5_MAIN = $(BASEDIR)/srcprog/stat/multiple_mechanisms/shuffling/generic_clones_feldman2023.cpp
# PROG6_MAIN = $(BASEDIR)/srcprog/stat/multiple_mechanisms/shuffling/generic_clones_ours2022.cpp

# --- Shared source directories --- a list
CODEDIRS= \
$(BASEDIR)/srclib/datasets/gdset \
$(BASEDIR)/srclib/dist \
$(BASEDIR)/srclib/loss \
$(BASEDIR)/srclib/mat \
$(BASEDIR)/srclib/mechanisms \
$(BASEDIR)/srclib/mechanisms/central \
$(BASEDIR)/srclib/mechanisms/line \
$(BASEDIR)/srclib/mechanisms/plan \
$(BASEDIR)/srclib/priors \
$(BASEDIR)/srclib/privMetrics \
$(BASEDIR)/srclib/spaces \
$(BASEDIR)/srclib/stat/tools \
$(BASEDIR)/srclib/stat/distributions \
$(BASEDIR)/srclib/stat/empirical \
$(BASEDIR)/srclib/stat/estimators

# Include directories (a list). 
# Note: $(shell brew --frefix glpk) translates to "/opt/homebrew/opt/glpk"
INCDIRS= \
$(BASEDIR)/include \
$(BASEDIR)/srclib/emd \
$(shell brew --prefix glpk)/include \
$(shell brew --prefix armadillo)/include \
$(shell brew --prefix libomp)/include 

# Library directories (may be a list if many item)
LIBDIRS = \
$(shell brew --prefix glpk)/lib \
$(shell brew --prefix armadillo)/lib \
$(shell brew --prefix libomp)/lib


CC = g++
# optimation options from O0 to O3
OPT = -O3 -mcpu=native -Xpreprocessor -fopenmp
# generate files that encode make rules for the .h dependencies
DEPFLAGS = -MP -MD

# C/C++ Compiler flags: automatically add the -I onto each include directory
CXXFLAGS = -std=c++17 -Wall -Wextra -g $(foreach D,$(INCDIRS),-I$(D)) $(OPT) $(DEPFLAGS)

# Linker flags: automatically add the -L onto each library directory
LDFLAGS= $(foreach D,$(LIBDIRS),-L$(D)) -lglpk -larmadillo -lomp -framework Accelerate 

# --- Shared sources (everything except the main programs) ---
# for-style iteration (foreach) and regular expression completions (wildcard)
COMMON_SRCS=$(foreach D,$(CODEDIRS),$(wildcard $(D)/*.cpp))
# a list of required object file names: obtained from COMMON_SRCS 
# pattern replacement '%'.
COMMON_OBJS=$(patsubst %.cpp,%.o,$(COMMON_SRCS))

# --- Program-specific objects (sources are placed at the top) ---
PROG1_OBJS = $(COMMON_OBJS) $(patsubst %.cpp,%.o,$(PROG1_MAIN))
PROG2_OBJS = $(COMMON_OBJS) $(patsubst %.cpp,%.o,$(PROG2_MAIN))
PROG3_OBJS = $(COMMON_OBJS) $(patsubst %.cpp,%.o,$(PROG3_MAIN))
PROG4_OBJS = $(COMMON_OBJS) $(patsubst %.cpp,%.o,$(PROG4_MAIN))
PROG5_OBJS = $(COMMON_OBJS) $(patsubst %.cpp,%.o,$(PROG5_MAIN))
PROG6_OBJS = $(COMMON_OBJS) $(patsubst %.cpp,%.o,$(PROG6_MAIN))

DEPFILES=$(patsubst %.cpp,%.d,$(COMMON_SRCS)) \
			$(patsubst %.cpp,%.d,$(PROG1_MAIN)) \
			$(patsubst %.cpp,%.d,$(PROG2_MAIN)) \
			$(patsubst %.cpp,%.d,$(PROG3_MAIN)) \
			$(patsubst %.cpp,%.d,$(PROG4_MAIN)) \
			$(patsubst %.cpp,%.d,$(PROG5_MAIN)) \
			$(patsubst %.cpp,%.d,$(PROG6_MAIN)) 

# --- Default target ---
all: $(BINARIES)

# linking object files:
# Important: link Flags must come after specifying output files 
# (resolved by the target $@) and the object files (resolved by the prerequisit variable $^)
$(BINARY1): $(PROG1_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(BINARY2): $(PROG2_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(BINARY3): $(PROG3_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(BINARY4): $(PROG4_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(BINARY5): $(PROG5_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(BINARY6): $(PROG6_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# Compiling source files: 
# Note that every cpp file (%.cpp) is translated 
# into one object file (%.o). $@ resolves the target and $< resolves the 
# first prerequisite, which is the only one (we don't need $^)
%.o:%.cpp
	$(CC) $(CXXFLAGS) -c -o $@ $<

# --- Clean ---
clean:
	rm -rf $(BINARIES) $(COMMON_OBJS) $(PROG1_OBJS) $(PROG2_OBJS) $(DEPFILES)
	rm -rf $(BINARIES) $(COMMON_OBJS) $(PROG3_OBJS) $(PROG4_OBJS) $(DEPFILES)
	rm -rf $(BINARIES) $(COMMON_OBJS) $(PROG5_OBJS) $(PROG6_OBJS) $(DEPFILES)

# shell commands are a set of keystrokes away
distribute: clean
	tar zcvf dist.tgz *

# @ silences the printing of the command
# $(info ...) prints output
diff:
	$(info The status of the repository, and the volume of per-file changes:)
	@git status
	@git diff --stat

# include the dependencies
-include $(DEPFILES)

# add .PHONY so that the non-targetfile - rules work even if a file with the same name exists.
.PHONY: all clean distribute diff
