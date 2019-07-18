//
//  main.c
//  SpMV-MPI
//
//  Created by Muhammed Cavusoglu on 15.07.2019.
//  Copyright © 2019 Muhammed Cavusoglu. All rights reserved.
//

#include "util.h"
#include "mpi.h"

int main(int argc, char * argv[]) {
    int num_procs, rank, num_iters;
    int *row_ptr, *col_ind, num_rows, num_cols, num_vals;
    float *values;
    
    const char *filename = "cavity02.mtx";
    num_iters = 5;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Read the matrix, initialize the vector
    if(rank == MASTER) {
        read_matrix(&row_ptr, &col_ind, &values, filename, &num_rows, &num_cols, &num_vals);
    }
    
    for (int i = 0; i < num_iters; i++) {
        if(rank == MASTER) {
            // Scatter matrix data (row_ptr, col_ind, values)
            
            // Scatter the vector
            
            // Determine remote vectors needed? Expand local vector array to receive those? Reindex local col_ind to point to the right vector entry
            
            
            
        } else { // Workers
            
        }
        
    }
    
    printf("Helloww: %d\n", rank);
    
    return 0;
}
