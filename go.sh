rm *.dat

EOS=GAMMA
GAmma=1.66666666666
GAmma_1=0.66666666666
gcc -Wall Main.c Fluid.c Functions.c RootFinder.c RarefactionWaves.c ShockWaves.c Utility.c Plot.c -I/software/gsl/default/include -c -DEOS=$EOS -DGAmma=$GAmma -DGAmma_1=$GAmma_1
gcc Main.o Fluid.o Functions.o RootFinder.o RarefactionWaves.c ShockWaves.c Utility.c Plot.c -lgsl -lgslcblas -lm -L/software/gsl/default/lib -I/software/gsl/default/include -DEOS=$EOS -DGAmma=$GAmma -DGAmma_1=$GAmma_1
./a.out

EOS=GAMMA
GAmma=1.3333333333
GAmma_1=0.3333333333
gcc -Wall Main.c Fluid.c Functions.c RootFinder.c RarefactionWaves.c ShockWaves.c Utility.c Plot.c -I/software/gsl/default/include -c -DEOS=$EOS -DGAmma=$GAmma -DGAmma_1=$GAmma_1
gcc Main.o Fluid.o Functions.o RootFinder.o RarefactionWaves.c ShockWaves.c Utility.c Plot.c -lgsl -lgslcblas -lm -L/software/gsl/default/lib -I/software/gsl/default/include -DEOS=$EOS -DGAmma=$GAmma -DGAmma_1=$GAmma_1
./a.out

EOS=TM
gcc -Wall Main.c Fluid.c Functions.c RootFinder.c RarefactionWaves.c ShockWaves.c Utility.c Plot.c -I/software/gsl/default/include -c -DEOS=$EOS -DGAmma=$GAmma -DGAmma_1=$GAmma_1
gcc Main.o Fluid.o Functions.o RootFinder.o RarefactionWaves.c ShockWaves.c Utility.c Plot.c -lgsl -lgslcblas -lm -L/software/gsl/default/lib -I/software/gsl/default/include -DEOS=$EOS -DGAmma=$GAmma -DGAmma_1=$GAmma_1
./a.out


rm *.o
