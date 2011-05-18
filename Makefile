define INPUT
100
100
100
0.01
100
endef

define INPUTMPI
100
100
100
0.01
100
2
2
2
endef

export INPUT
export INPUTMPI

all:
	$(MAKE) -C mpi
	$(MAKE) -C openmp
	$(MAKE) -C serial
	
test: clean all
	@# Write a title
	@echo
	@echo "Running tests..."
	@# Run the serial code and print the output
	@echo "$$INPUT" | ./serial/laplace | tail -n 1 | cut -d " " -f 3 > s_output
	
	@THREADS=8
	@export OMP_NUM_THREADS=THREADS
	@# Run the OpenMP code and print the output
	@echo "$$INPUT" | ./openmp/laplace | tail -n 1 | cut -d " " -f 3 > omp_output
	
	@# Run the MPI code and print the output
	@echo "$$INPUTMPI" | mpirun -np 8 ./mpi/laplace | tail -n 1 | cut -d " " -f 3 > mpi_output
	
	@./display_test_results.sh
	@rm mpi_output omp_output s_output
	
clean:
	@echo "Removing old binary files"
	@rm -f ./mpi/*.o
	@rm -f ./mpi/laplace
	@rm -f ./openmp/*.o
	@rm -f ./openmp/laplace
	@rm -f ./serial/*.o
	@rm -f ./serial/laplace