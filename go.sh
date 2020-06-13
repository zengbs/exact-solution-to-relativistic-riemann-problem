gcc Main.c Fluid.c Functions.c RootFinder.c RarefactionWaves.c ShockWaves.c Utility.c Plot.c -I/software/gsl/default/include -c -DEOS=GAMMA
gcc Main.o Fluid.o Functions.o RootFinder.o RarefactionWaves.c ShockWaves.c Utility.c Plot.c -lgsl -lgslcblas -lm -L/software/gsl/default/lib
rm *.o
