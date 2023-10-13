#====================================================================================================
# Solver options
#====================================================================================================
# Equation of state type: [TM/GAMMA]
SOL_OPTIONS += -DEOS=TM

# Adiabatic index of ideal gas
SOL_OPTIONS += -DGAmma=1.66666666666
SOL_OPTIONS += -DGAmma_1=0.66666666666



#====================================================================================================
# Macro definitions
#====================================================================================================
# compiler
CC     = gcc
CFLAGS =  -g

EXECUTABLE := SR_Riemann_solver

# paths
SRC_DIR := ./
OBJ_DIR := ./obj
INC_DIR := ./include
GSL_DIR := /software/gsl/default



#====================================================================================================
# Preliminary definitions
#====================================================================================================
SRC_FILES := $(shell find $(SRC_DIR) -name "*.c")
OBJ_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(SRC_FILES:.c=.o)))



#====================================================================================================
# Targets
#====================================================================================================
.PHONY: clean all compile

all: compile

compile: ${EXECUTABLE}

# link objects into executable
${EXECUTABLE}: ${OBJ_FILES}
	@printf "========================================\n"
	@printf "Linking to %s ... " $@
	@${CC}  -o ${EXECUTABLE} ${SOL_OPTIONS} ${OBJ_FILES} -I${INC_DIR} -lgsl -lgslcblas -lm -L${GSL_DIR}/lib -I${GSL_DIR}/include
	@printf "Successful!\n"

# create objects
${OBJ_DIR}/%.o:%.c
	@${CC} ${CFLAGS} ${SOL_OPTIONS} -I${INC_DIR} -I${GSL_DIR}/include -c $< -o $@
	@printf "Compiling %s\n" $<

# clean up
clean:
	rm -f ${OBJ_DIR}/*
	rm -f ${EXECUTABLE}
