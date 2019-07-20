spmv_mpi:
	python3 init_vector.py
	mpicc -o spmv_mpi main.c util.c -lm
clean:
	rm spmv_mpi
