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
    int num_procs, rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == MASTER) {
        printf("Hello, World!\n");
    } else { // Workers
        
    }
    
    return 0;
}
