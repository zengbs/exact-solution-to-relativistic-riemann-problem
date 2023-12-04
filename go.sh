rm -f *.dat

make clean
make
echo "Executing SR_Riemann_solver"
./SR_Riemann_solver
