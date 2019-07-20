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

// FIXME: CHANGE
int reindex(int *colInd, int numE, int myid, int numprocs, int nrows, int myVecSize,
                     int *remoteVecCount, int *remoteVecPtr, int **remoteVecIndex, float **myVecData);

int main(int argc, char * argv[]) {
    int num_procs, rank, num_rows;
    int local_matrix_item_count, *local_col_ind, *local_row_ptr;
    
    // FIXME: CHANGE
    float *local_matrix_values, *local_vector_data, *all_vector_data, *local_result, *total_result, *sendVecData;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // FIXME: CHANGE
    int remoteVecCount[num_procs], remoteVecPtr[num_procs];
    int tosendVecCount[num_procs], tosendVecPtr[num_procs];
    int tosendTotal = 0, numRemoteVec = 0;
    int *remoteVecIndex, *tosendVecIndex;

    if(rank == MASTER) {
        // Read the matrix into CSR format
        const char *filename = "cavity02.mtx";
        int *row_ptr, *col_ind, num_cols, num_vals;
        float *values;
        read_matrix(&row_ptr, &col_ind, &values, filename, &num_rows, &num_cols, &num_vals);
        
        // Read the vector
        float *vector_values = (float *) malloc(num_rows * sizeof(float));
        /*
         for (int i = 0; i < nrows; i++) {
            vector_values[i] = 1.0;
         }
         */
        char *vec=NULL;
        size_t vec_size;
        FILE *vector_file = fopen("x317.vec","r");
        for (int i = 0; i < num_rows; i++) {
            getline(&vec, &vec_size, vector_file);
            vector_values[i] = strtod(vec,NULL);
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
        local_vector_data = (float *)malloc(sizeof(float) * vector_item_count);
        
        // Scatter vector items
        MPI_Bcast(&num_rows, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        MPI_Scatterv(vector_values, vector_item_counts, vector_item_start_indices, MPI_FLOAT, local_vector_data, vector_item_count, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
        
        // FIXME: Scatter matrix entries to each processor
        // by sending partial Row pointers, Column Index and Values
        MPI_Scatter(matrix_item_counts,1,MPI_INT,&local_matrix_item_count,1,MPI_INT,0,MPI_COMM_WORLD);
        
        local_row_ptr = (int *)malloc(sizeof(int) * (vector_item_count+1));
        MPI_Scatterv(row_ptr,vector_item_counts,vector_item_start_indices,MPI_INT,
                     local_row_ptr,vector_item_count,MPI_INT,0,MPI_COMM_WORLD);
        local_row_ptr[vector_item_count] = local_matrix_item_count;
        
        local_col_ind = (int *)malloc(sizeof(int) * local_matrix_item_count);
        MPI_Scatterv(col_ind,matrix_item_counts,matrix_item_start_indices,MPI_INT,
                     local_col_ind,local_matrix_item_count,MPI_INT,0,MPI_COMM_WORLD);
        
        local_matrix_values = (float *)malloc(sizeof(float) * local_matrix_item_count);
        MPI_Scatterv(values,matrix_item_counts,matrix_item_start_indices,MPI_FLOAT,
                     local_matrix_values,local_matrix_item_count,MPI_FLOAT,0,MPI_COMM_WORLD);
        free(row_ptr); free(col_ind); free(values);
        
        //*** Step 2 - Analyze elements ***//
        // Determine the remote vectors needed, expand local Vector array to receive these.
        // Re-index local Column Index array to point to the right Vector Data entry
        numRemoteVec = reindex(local_col_ind,local_matrix_item_count,rank,num_procs,num_rows,vector_item_count,remoteVecCount,remoteVecPtr,&remoteVecIndex,&local_vector_data);
        
        //*** Step 3 & 4 - Determine which processors to send what ***//
        // Each process send to all others number of remote vector entries needed
        MPI_Alltoall(remoteVecCount,1,MPI_INT,tosendVecCount,1,MPI_INT,MPI_COMM_WORLD);
        
        for (int p = 0; p < num_procs; p++)
            tosendTotal = tosendTotal + tosendVecCount[p];
        tosendVecIndex = (int *)malloc(sizeof(int) * tosendTotal);
        tosendVecPtr[0] = 0;
        for (int p = 1; p < num_procs; p++)
            tosendVecPtr[p] = tosendVecPtr[p-1] + tosendVecCount[p-1];
        
        // Each process send to all others index of remote vector entries needed
        MPI_Alltoallv(remoteVecIndex,remoteVecCount,remoteVecPtr,MPI_INT,tosendVecIndex,
                      tosendVecCount,tosendVecPtr,MPI_INT,MPI_COMM_WORLD);
        
        // Get the vector entries other processors needed, ready to send
        sendVecData = (float *)malloc(sizeof(float) * tosendTotal);
        for (int i = 0; i < tosendTotal; i++)
            sendVecData[i] = local_vector_data[tosendVecIndex[i]];
        
        // Local Vector entries at the front, remote entries are concatenated starting from:
        float *recvVecData = local_vector_data + vector_item_count;
        
        //*** Step 5 - All processors send the actual vector entries to each other ***//
        MPI_Alltoallv(sendVecData,tosendVecCount,tosendVecPtr,MPI_FLOAT,recvVecData,
                      remoteVecCount,remoteVecPtr,MPI_FLOAT,MPI_COMM_WORLD);
        
        //*** Step 6 - Compute the results ***//
        local_result = (float *)malloc(sizeof(float) * vector_item_count);
        multiply(vector_item_count,local_row_ptr,local_col_ind,local_matrix_values,local_vector_data,local_result);
        
        //*** Step 7 - Gather results from all processors ***//
        total_result = (float *)malloc(sizeof(float) * num_rows);
        MPI_Gatherv(local_result,vector_item_count,MPI_FLOAT,total_result,vector_item_counts,vector_item_start_indices,
                    MPI_FLOAT,0,MPI_COMM_WORLD);
        
        // Creating output filename
        int c = 0, temp = num_rows;
        while (temp != 0) {
            temp = temp / 10;
            c++;
        }
        char fname[c+5];
        sprintf(fname,"o%d.vec",num_rows);
        
        // Writing results to output file
        FILE *fpout;
        fpout = fopen(fname,"w");
        for (int i = 0; i < num_rows; i++)
            fprintf(fpout,"%.15lg\n",total_result[i]);
        
        fclose(fpout);
        free(total_result);
        
    } else { // Workers
        
        int *colInd, *eCount, *eDispls, *rowptr, *vecDataSize, *vecDataDispls;
        float *matVal;
        int size;
        
        MPI_Bcast(&num_rows,1,MPI_INT,0,MPI_COMM_WORLD);
        int vector_item_count = ceil((float)num_rows / num_procs);
        size = vector_item_count;
        if (rank == num_procs-1)
            vector_item_count = num_rows - rank*vector_item_count;
        // Prepare to store local vector entries
        local_vector_data = (float *)malloc(sizeof(float) * vector_item_count);
        
        // Receive vector entries from root processor
        MPI_Scatterv(all_vector_data,vecDataSize,vecDataDispls,MPI_FLOAT,
                     local_vector_data,vector_item_count,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        MPI_Scatter(eCount,1,MPI_INT,&local_matrix_item_count,1,MPI_INT,0,MPI_COMM_WORLD);
        
        // Receive matrix entries from root processor
        local_row_ptr = (int *)malloc(sizeof(int) * (vector_item_count+1));
        MPI_Scatterv(rowptr,vecDataSize,vecDataDispls,MPI_INT,
                     local_row_ptr,vector_item_count,MPI_INT,0,MPI_COMM_WORLD);
        int rowptrOffset = local_row_ptr[0];
        for (int i = 0; i < vector_item_count; i++)
            local_row_ptr[i] = local_row_ptr[i] - rowptrOffset; //offset
        local_row_ptr[vector_item_count] = local_matrix_item_count;
        
        local_col_ind = (int *)malloc(sizeof(int) * local_matrix_item_count);
        MPI_Scatterv(colInd,eCount,eDispls,MPI_INT,
                     local_col_ind,local_matrix_item_count,MPI_INT,0,MPI_COMM_WORLD);
        
        local_matrix_values = (float *)malloc(sizeof(float) * local_matrix_item_count);
        MPI_Scatterv(matVal,eCount,eDispls,MPI_FLOAT,
                     local_matrix_values,local_matrix_item_count,MPI_FLOAT,0,MPI_COMM_WORLD);
        
        // Determine the remote vectors needed, expand local Vector array to receive these.
        // Re-index local Column Index array to point to the right Vector Data entry
        numRemoteVec = reindex(local_col_ind,local_matrix_item_count,rank,num_procs,num_rows,vector_item_count,remoteVecCount,remoteVecPtr,&remoteVecIndex,&local_vector_data);
        
        // Each process send to all others number of remote vector entries needed
        MPI_Alltoall(remoteVecCount,1,MPI_INT,tosendVecCount,1,MPI_INT,MPI_COMM_WORLD);
        
        for (int p = 0; p < num_procs; p++)
            tosendTotal = tosendTotal + tosendVecCount[p];
        tosendVecIndex = (int *)malloc(sizeof(int) * tosendTotal);
        tosendVecPtr[0] = 0;
        for (int p = 1; p < num_procs; p++)
            tosendVecPtr[p] = tosendVecPtr[p-1] + tosendVecCount[p-1];
        
        // Each process send to all others index of remote vector entries needed
        MPI_Alltoallv(remoteVecIndex,remoteVecCount,remoteVecPtr,MPI_INT,
                      tosendVecIndex,tosendVecCount,tosendVecPtr,MPI_INT,MPI_COMM_WORLD);
        // Get the vector entries other processors needed, ready to send
        sendVecData = (float *)malloc(sizeof(float) * tosendTotal);
        for (int i = 0; i < tosendTotal; i++)
            sendVecData[i] = local_vector_data[tosendVecIndex[i]-rank*size];
        
        // Local Vector entries at the front, remote entries are concatenated starting from:
        float *recvVecData = local_vector_data + vector_item_count;
        
        // All processors send the actual vector entries to each other
        MPI_Alltoallv(sendVecData,tosendVecCount,tosendVecPtr,MPI_FLOAT,recvVecData,
                      remoteVecCount,remoteVecPtr,MPI_FLOAT,MPI_COMM_WORLD);
        
        // Compute the local results
        local_result = (float *)malloc(sizeof(float) * vector_item_count);
        multiply(vector_item_count,local_row_ptr,local_col_ind,local_matrix_values,local_vector_data,local_result);
        
        // Send the results to the root processor
        MPI_Gatherv(local_result,vector_item_count,MPI_FLOAT,total_result,vecDataSize,vecDataDispls,
                    MPI_FLOAT,0,MPI_COMM_WORLD);
    }
    
    free(local_row_ptr); free(local_col_ind); free(local_matrix_values); free(tosendVecIndex);
    free(sendVecData); free(local_result); free(local_vector_data); free(remoteVecIndex);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
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

// Append remote vector items needed to the local vector, and re-index local col_ind
// FIXME: CHANGE PARAMS
int reindex(int *local_col_ind, int myNumE, int rank, int num_procs, int num_rows, int myVecSize, int *remoteVecCount, int *remoteVecPtr, int **remoteVecIndex, float **myVecData) {
    int size = ceil((float) num_rows / num_procs);
    int local_vector_start_index = rank * size;
    int local_vector_end_index = (rank + 1) * size;
    
    // FIXME: CHANGE
    int tempRemoteVecCount[num_procs];
    
    // FIXME: CHANGE
    for (int i = 0; i < num_procs; i++)
        remoteVecCount[i] = 0;
    *remoteVecIndex = (int *)malloc(sizeof(int) * myNumE);
    
    // Select remote vectors (counts, indices, displacements from each remote processor)
    int numRemoteVec = 0;
    for (int i = 0; i < myNumE; i++) {
        if ((local_col_ind[i] < local_vector_start_index) || (local_col_ind[i] >= local_vector_end_index)) {
            remoteVecCount[local_col_ind[i]/size]++;
            numRemoteVec++;
        }
    }
    remoteVecPtr[0] = 0;
    tempRemoteVecCount[0] = 0;
    for (int i = 1; i < num_procs; i++) {
        remoteVecPtr[i] = remoteVecPtr[i-1] + remoteVecCount[i-1];
        tempRemoteVecCount[i]=0;
    }
    
    // Realloc local vector array to store received remote vector entries
    *myVecData = (float *)realloc(*myVecData,sizeof(float) * (myVecSize+numRemoteVec));
    
    // Re-index the local Column Index array to point to local Vector Data entry at the beginning and remote Vector Data entries at the end of Vector Data array.
    int proc = 0;
    for (int i = 0; i < myNumE; i++) {
        if ((local_col_ind[i] < local_vector_start_index) || (local_col_ind[i] >= local_vector_end_index)) {
            proc = local_col_ind[i]/size;
            (*remoteVecIndex)[remoteVecPtr[proc]+tempRemoteVecCount[proc]] = local_col_ind[i];
            
            local_col_ind[i] = myVecSize + remoteVecPtr[proc] + tempRemoteVecCount[proc];
            tempRemoteVecCount[proc]++;
        }
        else {
            local_col_ind[i] -= (rank * size);
        }
    }
    
    //return total number of remote vector entries needed
    return numRemoteVec;
}
