#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

LOCAL_LIB_DIR=$(MESA_DIR)/lib

INCLUDE = -I$(MESA_DIR)/include
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -linterp_2d
LOAD_OTHER = -L$(MESA_LIB_DIR) -lnum -linterp_1d -lutils -lconst $(LOAD_MATRIX)

#################################################################

# STEP 3: build

FWARN=-pedantic -fcheck=all -fbacktrace -g
FFLAGS=-O2 $(FWARN)

LIB = libkapCN.a
LIB_OBJS = mod_kapCN.o

TEST = ../test_kapCN 
TEST_OBJS = test_kapCN.o $(LIB)

all: $(LIB) $(TEST)

$(LIB): $(LIB_OBJS)
	ar crs $(LIB) $(LIB_OBJS)

$(TEST): $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST) $^ $(INCLUDE) $(LOAD_LOCAL) $(LOAD_OTHER)

%.o: %.f90
	$(COMPILE_TO_TEST) $(FCopenmp) $(FCfree) $(INCLUDE) -c $<

check: $(TEST)
	./$<

clean:
	rm -f *.mod *.o *~ $(LIB) $(TEST)
