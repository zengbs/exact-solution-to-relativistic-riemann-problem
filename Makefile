# object files-------------------------------------------------------

CORE_OBJ = Fluid.o Functions.o Main.o RootFinder.o Utility.o ShockWaves.o

ALL_OBJ = ${CORE_OBJ}

# macro ddefinitions------------------------------------------------
CC = gcc

CFLAGS =  -g


DLIB = /opt/math/gsl/default/lib/
DINC = /opt/math/gsl/default/include/

BIN = a.out

SRC = $(ALL_OBJ:.o=.c)


# implicit rules-----------------------------------------------------
%.o:%.c
	${CC}  ${CFLAGS} -I${DINC} -c $<

# targets------------------------------------------------------------

.PHONY: clean all compile

all: compile

compile: ${BIN}

${BIN}: ${ALL_OBJ}
	${CC}  -o ${EXEDIR}a.out ${ALL_OBJ} -lm ${DLIB}libgslcblas.so.0.0.0 ${DLIB}libgsl.so.23.1.0    # shared library
clean:
	rm -f *.o Makedepend $(BIN)
