# solver options
SOL_OPTIONS += -DEOS=TM
SOL_OPTIONS += -DGAmma=1.66666666666
SOL_OPTIONS += -DGAmma_1=0.66666666666

# object files-------------------------------------------------------

CORE_OBJ = Fluid.o Functions.o Main.o RootFinder.o Utility.o ShockWaves.o RarefactionWaves.o Plot.o

ALL_OBJ = ${CORE_OBJ}

# macro ddefinitions------------------------------------------------
CC = gcc

CFLAGS =  -g


#DLIB = /opt/math/gsl/default/lib/
#DINC = /opt/math/gsl/default/include/
DLIB = /software/gsl/default/lib/
DINC = /software/gsl/default/include/

EXECUTABLE = a.out

SRC = $(ALL_OBJ:.o=.c)


# implicit rules-----------------------------------------------------
%.o:%.c
	${CC}  ${CFLAGS} ${SOL_OPTIONS} -I${DINC} -c $<

# targets------------------------------------------------------------

.PHONY: clean all compile

all: compile

compile: ${EXECUTABLE}

${EXECUTABLE}: ${ALL_OBJ}
#	${CC}  -o ${EXEDIR}a.out ${ALL_OBJ} -lm ${DLIB}libgslcblas.so.0.0.0 ${DLIB}libgsl.so.23.1.0    # shared library
	${CC}  -o ${EXECUTABLE} ${SOL_OPTION} ${ALL_OBJ} -lgsl -lgslcblas -lm -L/software/gsl/default/lib -I/software/gsl/default/include
clean:
	rm -f *.o Makedepend ${EXECUTABLE}
