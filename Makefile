#====================================================================================================
# Solver Options
#====================================================================================================
# TM:    Taub-Mathews
# GAMMA: Constant-$\Gamma$

# Type of Equation of state: [TM/GAMMA]
SOL_OPTIONS += -DEOS=TM


# Adiabatic index of ideal gas
SOL_OPTIONS += -DGAmma=1.66666666666
SOL_OPTIONS += -DGAmma_1=0.66666666666

#====================================================================================================
# Macro Definitions
#====================================================================================================

# Compiler
CC   := g++

CC_FLAGS      = -lm --std=c++1z -O3 -lgsl -lgslcblas
CC_DEBUG_FLAG = -g -Wall


# GSL paths
GSL_DIR := /software/gsl/default
INC := -I${GSL_DIR}/include
LIB := -L${GSL_DIR}/lib


EXECUTABLE := RelativisticRiemannSolver



## Source files
CC_FILE = fluid.c functions.c main.c plot.c rarefaction_waves.c root_finder.c \
          shock_waves.c utilities.c load_parameter.c
SRC_PATH = src

## Object files
CC_OBJ_FILE = $(CC_FILE:.c=.o)
OBJ_PATH    = objective
CC_OBJ      = $(patsubst %,$(OBJ_PATH)/%,$(CC_OBJ_FILE))


# Header files
HEADER_FILE = macro.h prototypes.h global.h struct.h
INC_PATH    = include
HEADER      = $(patsubst %,$(INC_PATH)/%,$(HEADER_FILE))

CC_ALL_FLAG   = $(CC_DEBUG_FLAG)  $(INC) $(CC_FLAGS) $(SOL_OPTIONS)

# Compiling
$(OBJ_PATH)/%.o:$(SRC_PATH)/%.c $(HEADER)
	$(CC) $(CC_ALL_FLAG) -o $@ -c $<


# Linking
$(EXECUTABLE): $(CC_OBJ)
	$(CC) -o $@ $^  $(LIB)  -lm -lgsl -lgslcblas
	mv $(EXECUTABLE) bin/

.PHONY:
clean:
	@rm -rf ${OBJ_PATH}/*.o rm bin/${EXECUTABLE}
