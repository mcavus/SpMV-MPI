//
//  main.c
//  SpMV-MPI
//
//  Created by Muhammed Cavusoglu on 15.07.2019.
//  Copyright © 2019 Muhammed Cavusoglu. All rights reserved.
//

#include "util.h"
#include "mpi.h"

void multiply(int num_rows, int *row_ptr, int *col_ind, float *values, float *vector, float *result);
int reindex(int *local_col_ind, int *remote_vector_item_counts, int *remote_vector_ptr, int **remote_vector_item_indices, float **my_vector_data, int local_matrix_item_count, int rank, int num_procs, int num_rows, int local_vector_size);

int main(int argc, char * argv[]) {
    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (argc != 2) {
        if (rank == MASTER) {
            printf("Invalid command, enter:\n1. number of repetitions\n");
        }
    } else {
        int num_iters = atoi(argv[1]);
        for (int iter = 0; iter < num_iters; iter++) {
            int num_rows, local_matrix_item_count, *local_col_ind, *local_row_ptr;
            float *local_matrix_values, *local_vector_data, *all_vector_data, *local_result, *total_result;
            
            int remote_vector_item_counts[num_procs], remote_vector_ptr[num_procs], to_send_vector_item_counts[num_procs], to_send_vector_ptr[num_procs];
            int to_send_total = 0, remote_vector_item_count = 0, *remote_vector_item_indices, *to_send_vector_item_indices;
            float *to_send_vector_data;
            
            if(rank == MASTER) {
                // Read the matrix into CSR format
                // const char *filename = argv[3];
                const char *filename = "cavity02.mtx";
                int *row_ptr, *col_ind, num_cols, num_vals;
                float *values;
                read_matrix(&row_ptr, &col_ind, &values, filename, &num_rows, &num_cols, &num_vals);
                
                // Read the vector
                float *vector_values = malloc(num_rows * sizeof(float));
                /*
                 for (int i = 0; i < num_rows; i++) {
                 vector_values[i] = 1.0;
                 }
                 */
                char *vec=NULL;
                size_t vec_size;
                FILE *vector_file = fopen("x317.vec","r");
                for (int i = 0; i < num_rows; i++) {
                    getline(&vec, &vec_size, vector_file);
                    vector_values[i] = strtod(vec, NULL);
                }
                
                // Vector items and starting indices of each processor
                int vector_item_count = ceil((float) num_rows / num_procs);
                int vector_item_count_last = num_rows - (num_procs - 1) * vector_item_count; // for the last processor
                
                int vector_item_counts[num_procs], vector_item_start_indices[num_procs];
                for (int p = 0; p < num_procs-1; p++) {
                    vector_item_counts[p] = vector_item_count;
                    vector_item_start_indices[p] = p * vector_item_count;
                }
                vector_item_counts[num_procs - 1] = vector_item_count_last;
                vector_item_start_indices[num_procs - 1] = (num_procs - 1) * vector_item_count;
                
                // Matrix items and starting indices of each processor
                int matrix_item_counts[num_procs], matrix_item_start_indices[num_procs];
                for (int p = 0; p < num_procs; p++) {
                    matrix_item_counts[p] = row_ptr[p * vector_item_count + vector_item_counts[p]] - row_ptr[p * vector_item_count];
                    matrix_item_start_indices[p] = row_ptr[p * vector_item_count];
                }
                
                // Local vector items
                local_vector_data = malloc(vector_item_count * sizeof(float));
                
                // Scatter vector items
                MPI_Bcast(&num_rows, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                MPI_Scatterv(vector_values, vector_item_counts, vector_item_start_indices, MPI_FLOAT, local_vector_data, vector_item_count, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
                
                // Scatter the matrix
                MPI_Scatter(matrix_item_counts, 1, MPI_INT, &local_matrix_item_count, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                local_row_ptr = malloc((vector_item_count + 1) * sizeof(int));
                MPI_Scatterv(row_ptr, vector_item_counts, vector_item_start_indices, MPI_INT, local_row_ptr, vector_item_count, MPI_INT, MASTER, MPI_COMM_WORLD);
                local_row_ptr[vector_item_count] = local_matrix_item_count;
                local_col_ind = malloc(local_matrix_item_count* sizeof(int));
                MPI_Scatterv(col_ind, matrix_item_counts, matrix_item_start_indices, MPI_INT, local_col_ind, local_matrix_item_count, MPI_INT, MASTER, MPI_COMM_WORLD);
                local_matrix_values = malloc(local_matrix_item_count * sizeof(float));
                MPI_Scatterv(values, matrix_item_counts, matrix_item_start_indices, MPI_FLOAT, local_matrix_values, local_matrix_item_count, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
                free(row_ptr); free(col_ind); free(values);
                
                // Append remote vector items needed to the local vector,
                // Re-index local col_ind
                // Return total number of remote vector items needed
                remote_vector_item_count = reindex(local_col_ind, remote_vector_item_counts, remote_vector_ptr, &remote_vector_item_indices, &local_vector_data, local_matrix_item_count, rank, num_procs, num_rows, vector_item_count);
                
                // All to all comm. for the number of remote vector items needed
                MPI_Alltoall(remote_vector_item_counts,1,MPI_INT,to_send_vector_item_counts,1,MPI_INT,MPI_COMM_WORLD);
                
                for (int p = 0; p < num_procs; p++)
                    to_send_total = to_send_total + to_send_vector_item_counts[p];
                to_send_vector_item_indices = malloc(to_send_total * sizeof(int));
                to_send_vector_ptr[0] = 0;
                for (int p = 1; p < num_procs; p++)
                    to_send_vector_ptr[p] = to_send_vector_ptr[p-1] + to_send_vector_item_counts[p-1];
                
                // All to all comm. for the index of remote vector items needed
                MPI_Alltoallv(remote_vector_item_indices, remote_vector_item_counts, remote_vector_ptr, MPI_INT, to_send_vector_item_indices, to_send_vector_item_counts, to_send_vector_ptr, MPI_INT, MPI_COMM_WORLD);
                
                // Vector items to be sent
                to_send_vector_data = malloc(to_send_total * sizeof(float));
                for (int i = 0; i < to_send_total; i++)
                    to_send_vector_data[i] = local_vector_data[to_send_vector_item_indices[i]];
                
                // Remote vector items are appended to local vector items
                float *recvVecData = local_vector_data + vector_item_count;
                
                // All to all comm. for the vector items
                MPI_Alltoallv(to_send_vector_data, to_send_vector_item_counts, to_send_vector_ptr, MPI_FLOAT, recvVecData, remote_vector_item_counts, remote_vector_ptr, MPI_FLOAT, MPI_COMM_WORLD);
                
                // Compute local result
                local_result = malloc(vector_item_count * sizeof(float));
                multiply(vector_item_count, local_row_ptr, local_col_ind, local_matrix_values, local_vector_data, local_result);
                
                // Gather local results of other processors
                total_result = malloc(num_rows * sizeof(float));
                MPI_Gatherv(local_result, vector_item_count, MPI_FLOAT, total_result, vector_item_counts, vector_item_start_indices, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
                
                // Write results to output file
                FILE *output_file = fopen("x317.vec","w");
                for (int i = 0; i < num_rows; i++)
                    fprintf(output_file, "%.5lg\n", total_result[i]);
                fclose(output_file);
                
                free(vector_values);
                free(total_result);
                
            } else { // Workers
                int *col_ind, *matrix_item_count, *matrix_item_indices, *row_ptr, *vector_item_counts, *vector_item_start_indices;
                float *matrix_values;
                
                MPI_Bcast(&num_rows, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                int vector_item_count = ceil((float) num_rows / num_procs);
                int size = vector_item_count;
                if (rank == num_procs-1)
                    vector_item_count = num_rows - (rank * vector_item_count);
                
                // Local vector items
                local_vector_data = malloc(vector_item_count * sizeof(float));
                
                // Receive vector items from master
                MPI_Scatterv(all_vector_data, vector_item_counts, vector_item_start_indices, MPI_FLOAT, local_vector_data, vector_item_count, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
                
                MPI_Scatter(matrix_item_count, 1, MPI_INT, &local_matrix_item_count, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                
                // Receive matrix items from master
                local_row_ptr = malloc((vector_item_count + 1) * sizeof(int));
                MPI_Scatterv(row_ptr, vector_item_counts, vector_item_start_indices, MPI_INT, local_row_ptr, vector_item_count, MPI_INT, MASTER, MPI_COMM_WORLD);
                int row_ptr_offset = local_row_ptr[0];
                for (int i = 0; i < vector_item_count; i++)
                    local_row_ptr[i] = local_row_ptr[i] - row_ptr_offset;
                local_row_ptr[vector_item_count] = local_matrix_item_count;
                
                local_col_ind = malloc(local_matrix_item_count * sizeof(int));
                MPI_Scatterv(col_ind, matrix_item_count, matrix_item_indices, MPI_INT, local_col_ind, local_matrix_item_count, MPI_INT, MASTER, MPI_COMM_WORLD);
                
                local_matrix_values = malloc(local_matrix_item_count * sizeof(float));
                MPI_Scatterv(matrix_values, matrix_item_count, matrix_item_indices, MPI_FLOAT, local_matrix_values, local_matrix_item_count, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
                
                // Append remote vector items needed to the local vector,
                // Re-index local col_ind
                // Return total number of remote vector items needed
                remote_vector_item_count = reindex(local_col_ind, remote_vector_item_counts, remote_vector_ptr, &remote_vector_item_indices, &local_vector_data, local_matrix_item_count, rank, num_procs, num_rows, vector_item_count);
                
                // All to all comm. for the number of remote vector items needed
                MPI_Alltoall(remote_vector_item_counts, 1, MPI_INT, to_send_vector_item_counts, 1, MPI_INT, MPI_COMM_WORLD);
                
                for (int p = 0; p < num_procs; p++)
                    to_send_total = to_send_total + to_send_vector_item_counts[p];
                to_send_vector_item_indices = malloc(to_send_total * sizeof(int));
                to_send_vector_ptr[0] = 0;
                for (int p = 1; p < num_procs; p++)
                    to_send_vector_ptr[p] = to_send_vector_ptr[p-1] + to_send_vector_item_counts[p-1];
                
                // All to all comm. for the index of remote vector items needed
                MPI_Alltoallv(remote_vector_item_indices, remote_vector_item_counts, remote_vector_ptr, MPI_INT, to_send_vector_item_indices, to_send_vector_item_counts, to_send_vector_ptr, MPI_INT, MPI_COMM_WORLD);
                
                // Vector items to be sent
                to_send_vector_data = (float *)malloc(sizeof(float) * to_send_total);
                for (int i = 0; i < to_send_total; i++)
                    to_send_vector_data[i] = local_vector_data[to_send_vector_item_indices[i]-rank*size];
                
                // Remote vector items are appended to local vector items
                float *recv_vector_data = local_vector_data + vector_item_count;
                
                // All to all comm. for the vector items
                MPI_Alltoallv(to_send_vector_data, to_send_vector_item_counts, to_send_vector_ptr, MPI_FLOAT, recv_vector_data, remote_vector_item_counts, remote_vector_ptr, MPI_FLOAT, MPI_COMM_WORLD);
                
                // Compute local result
                local_result = malloc(vector_item_count * sizeof(float));
                multiply(vector_item_count, local_row_ptr, local_col_ind, local_matrix_values, local_vector_data, local_result);
                
                // Send local result to master
                MPI_Gatherv(local_result, vector_item_count, MPI_FLOAT, total_result, vector_item_counts, vector_item_start_indices, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
            }
            
            MPI_Barrier(MPI_COMM_WORLD); //Check all the mallocs to free
            
            // Clean up
            free(local_row_ptr);
            free(local_col_ind);
            free(local_matrix_values);
            free(local_vector_data);
            free(to_send_vector_item_indices);
            free(to_send_vector_data);
            free(local_result);
            free(remote_vector_item_indices);
        }
    }
    
    MPI_Finalize();
    return 0;
}

// CSR Matrix-Vector Multiplication
void multiply(int num_rows, int *row_ptr, int *col_ind, float *values, float *vector, float *result) {
    int count = 0;
    for (int i = 0; i < num_rows; i++) {
        result[i] = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            result[i] += values[count] * vector[col_ind[j]];
            count++;
        }
    }
}

// Append remote vector items needed to the local vector,
// Re-index local col_ind
// Return total number of remote vector items needed
int reindex(int *local_col_ind, int *remote_vector_item_counts, int *remote_vector_ptr, int **remote_vector_item_indices, float **my_vector_data, int local_matrix_item_count, int rank, int num_procs, int num_rows, int local_vector_size) {
    int size = ceil((float) num_rows / num_procs);
    int local_vector_start_index = rank * size;
    int local_vector_end_index = (rank + 1) * size;
    
    int temp_remote_vector_item_counts[num_procs];
    for (int i = 0; i < num_procs; i++)
        remote_vector_item_counts[i] = 0;
    *remote_vector_item_indices = malloc(local_matrix_item_count * sizeof(int));
    
    // Select remote vector items
    int remote_vec_item_count = 0;
    for (int i = 0; i < local_matrix_item_count; i++) {
        if ((local_col_ind[i] < local_vector_start_index) || (local_col_ind[i] >= local_vector_end_index)) {
            remote_vector_item_counts[local_col_ind[i]/size]++;
            remote_vec_item_count++;
        }
    }
    remote_vector_ptr[0] = 0;
    temp_remote_vector_item_counts[0] = 0;
    for (int i = 1; i < num_procs; i++) {
        remote_vector_ptr[i] = remote_vector_ptr[i-1] + remote_vector_item_counts[i-1];
        temp_remote_vector_item_counts[i]=0;
    }
    
    // Expand local vector array to store remote vector items
    *my_vector_data = realloc(*my_vector_data, (local_vector_size+remote_vec_item_count) * sizeof(float));
    
    // Re-index local col_ind to point to local vector items and the remote items appended to the end
    int proc = 0;
    for (int i = 0; i < local_matrix_item_count; i++) {
        if ((local_col_ind[i] < local_vector_start_index) || (local_col_ind[i] >= local_vector_end_index)) {
            proc = local_col_ind[i] / size;
            (*remote_vector_item_indices)[remote_vector_ptr[proc] + temp_remote_vector_item_counts[proc]] = local_col_ind[i];
            
            local_col_ind[i] = local_vector_size + remote_vector_ptr[proc] + temp_remote_vector_item_counts[proc];
            temp_remote_vector_item_counts[proc]++;
        } else {
            local_col_ind[i] -= (rank * size);
        }
    }
    
    // Return total number of remote vector items needed
    return remote_vec_item_count;
}
