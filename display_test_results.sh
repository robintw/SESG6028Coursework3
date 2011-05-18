S_RES=`cat s_output`
O_RES=`cat omp_output`
M_RES=`cat mpi_output`

echo "============================="
echo "======= Test Results ========"
echo "============================="
echo -e "Serial result\t" $S_RES
echo -e "OpenMP result\t" $O_RES
echo -e "MPI result\t" $M_RES
echo "============================="
