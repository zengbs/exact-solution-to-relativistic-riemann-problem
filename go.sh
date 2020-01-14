gcc Main.c Fluid.c Functions.c RootFinder.c ContactWaves.c RarefactionWaves.c ShockWaves.c Utility.c Plot.c -I/usr/local/include/gsl -c
gcc Main.o Fluid.o Functions.o RootFinder.o ContactWaves.c RarefactionWaves.c ShockWaves.c Utility.c Plot.c -lgsl -lgslcblas -lm
rm *.o
