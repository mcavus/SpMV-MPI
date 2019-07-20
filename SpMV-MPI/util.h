//
//  util.h
//  SpMV-MPI
//
//  Created by Muhammed Cavusoglu on 15.07.2019.
//  Copyright © 2019 Muhammed Cavusoglu. All rights reserved.
//

#ifndef util_h
#define util_h

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define MASTER 0

void read_matrix(int **row_ptr, int **col_ind, float **values, const char *filename, int *num_rows, int *num_cols, int *num_vals);

#endif /* util_h */
