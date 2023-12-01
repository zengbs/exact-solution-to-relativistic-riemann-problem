MACRO :=  -DGAmma=1.66666666666 -DGAmma_1=0.66666666666 -DEOS=TM


CC   := g++

CC_FLAGS      = -lm --std=c++1z -O3 -lgsl -lgslcblas
CC_DEBUG_FLAG = -g -Wall


INC := -I/home/rocky/softwares/gsl/include
LIB := -L/home/rocky/softwares/gsl/lib64


EXECUTABLE := a.out


## Source files
CC_FILE = fluid.c functions.c main.c plot.c rarefaction_waves.c root_finder.c shock_waves.c utilities.c
SRC_PATH = src
CC_SRC = $(patsubst %,$(SRC_PATH)/%,$(CC_FILE))

## Object files
CC_OBJ_FILE = $(CC_FILE:.c=.o)
OBJ_PATH    = objective
CC_OBJ      = $(patsubst %,$(OBJ_PATH)/%,$(CC_OBJ_FILE))


# Header files
HEADER_FILE = macro.h prototypes.h global.h struct.h
INC_PATH    = includes
HEADER      = $(patsubst %,$(INC_PATH)/%,$(HEADER_FILE))

CC_ALL_FLAG   = $(CC_DEBUG_FLAG)  $(INC) $(CC_FLAGS) $(MACRO)

# Compiling
$(OBJ_PATH)/%.o : $(SRC_PATH)/%.c $(HEADER)
	$(CC) -o $@ -c $<  $(CC_ALL_FLAG)


# Linking
$(EXECUTABLE): $(CC_OBJ)
	g++ -o $@ $^  $(LIB)  -lm -lgsl -lgslcblas
	mv $(EXECUTABLE) bin/

.PHONY:
clean:
	@rm -rf ${OBJ_PATH}/*.o rm ${EXECUTABLE}

