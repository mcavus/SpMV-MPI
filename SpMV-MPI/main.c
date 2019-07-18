//
//  main.c
//  SpMV-MPI
//
//  Created by Muhammed Cavusoglu on 15.07.2019.
//  Copyright © 2019 Muhammed Cavusoglu. All rights reserved.
//

#include "util.h"
#include "mpi.h"

int prepareRemoteVec(int *colInd, int numE, int myid, int numprocs, int nrows, int myVecSize,
                     int *remoteVecCount, int *remoteVecPtr, int **remoteVecIndex, double **myVecData);
void computeMatVec(int proc, int n, int *rowptr, int *colind, double *matval, double *vec, double *res);

int main(int argc, char * argv[]) {
    int num_procs, myid, name_len, p=0, i=0, nrows, myNumE, vSize=0, vSize2=0;
    double *myVecData, *myMatVal, *result, *totalResult, *sendVecData;
    int *myColInd, *myRowptr;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    int remoteVecCount[num_procs], remoteVecPtr[num_procs];
    int tosendVecCount[num_procs], tosendVecPtr[num_procs];
    int tosendTotal = 0, numRemoteVec = 0;
    int *remoteVecIndex, *tosendVecIndex;
    
    int *rowptr, *colInd, num_cols, num_vals;
    double *values;
    double *allVecData;
    int vecDataSize[num_procs], vecDataDispls[num_procs];
    int eCount[num_procs]; int eDispls[num_procs];
    
    // First Process - Input & Output / Distribute Data & Receive Results / Compute
    if ( myid == 0 ) {
        
        const char *filename = "cavity02.mtx";
        read_matrix(&rowptr, &colInd, &values, filename, &nrows, &num_cols, &num_vals);
        
        allVecData = (double *) malloc(nrows * sizeof(double));
        for (i = 0; i < nrows; i++) {
            allVecData[i] = 1.0;
            //y[i] = 0.0;
        }
        
        vSize = ceil((double)nrows / num_procs); //partial vector size for each
        vSize2 = nrows - (num_procs-1)*vSize;    //vector size for last process
        
        for (p = 0; p < num_procs-1; p++) {
            vecDataSize[p] = vSize;
            vecDataDispls[p] = p*vSize;
        }
        vecDataSize[num_procs-1] = vSize2;
        vecDataDispls[num_procs-1] = (num_procs-1)*vSize;
        
        // Matrix entries count and displacement for each processor
        for (p = 0; p < num_procs; p++) {
            eCount[p] = rowptr[p*vSize+vecDataSize[p]] - rowptr[p*vSize];
            eDispls[p] = rowptr[p*vSize];
        }
    }
    
    for (int iter = 0; iter < 1; iter++) {
        if(myid == MASTER) {
            
            // Prepare to store local vector entries
            myVecData = (double *)malloc(sizeof(double) * vSize);
            
            // Scatter vector entries to each processor
            MPI_Bcast(&nrows,1,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Scatterv(allVecData,vecDataSize,vecDataDispls,MPI_DOUBLE,
                         myVecData,vSize,MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            // Scatter matrix entries to each processor
            // by sending partial Row pointers, Column Index and Values
            MPI_Scatter(eCount,1,MPI_INT,&myNumE,1,MPI_INT,0,MPI_COMM_WORLD);
            
            myRowptr = (int *)malloc(sizeof(int) * (vSize+1));
            MPI_Scatterv(rowptr,vecDataSize,vecDataDispls,MPI_INT,
                         myRowptr,vSize,MPI_INT,0,MPI_COMM_WORLD);
            myRowptr[vSize] = myNumE;
            
            myColInd = (int *)malloc(sizeof(int) * myNumE);
            MPI_Scatterv(colInd,eCount,eDispls,MPI_INT,
                         myColInd,myNumE,MPI_INT,0,MPI_COMM_WORLD);
            
            myMatVal = (double *)malloc(sizeof(double) * myNumE);
            MPI_Scatterv(values,eCount,eDispls,MPI_DOUBLE,
                         myMatVal,myNumE,MPI_DOUBLE,0,MPI_COMM_WORLD);
            free(rowptr); free(colInd); free(values);
            
            //*** Step 2 - Analyze elements ***//
            // Determine the remote vectors needed, expand local Vector array to receive these.
            // Re-index local Column Index array to point to the right Vector Data entry
            numRemoteVec = prepareRemoteVec(myColInd,myNumE,myid,num_procs,nrows,vSize,remoteVecCount,remoteVecPtr,&remoteVecIndex,&myVecData);
            
            //*** Step 3 & 4 - Determine which processors to send what ***//
            // Each process send to all others number of remote vector entries needed
            MPI_Alltoall(remoteVecCount,1,MPI_INT,tosendVecCount,1,MPI_INT,MPI_COMM_WORLD);
            
            for (p = 0; p < num_procs; p++)
                tosendTotal = tosendTotal + tosendVecCount[p];
            tosendVecIndex = (int *)malloc(sizeof(int) * tosendTotal);
            tosendVecPtr[0] = 0;
            for (p = 1; p < num_procs; p++)
                tosendVecPtr[p] = tosendVecPtr[p-1] + tosendVecCount[p-1];
            
            // Each process send to all others index of remote vector entries needed
            MPI_Alltoallv(remoteVecIndex,remoteVecCount,remoteVecPtr,MPI_INT,tosendVecIndex,
                          tosendVecCount,tosendVecPtr,MPI_INT,MPI_COMM_WORLD);
            
            // Get the vector entries other processors needed, ready to send
            sendVecData = (double *)malloc(sizeof(double) * tosendTotal);
            for (i = 0; i < tosendTotal; i++)
                sendVecData[i] = myVecData[tosendVecIndex[i]];
            
            // Local Vector entries at the front, remote entries are concatenated starting from:
            double *recvVecData = myVecData + vSize;
            
            //*** Step 5 - All processors send the actual vector entries to each other ***//
            MPI_Alltoallv(sendVecData,tosendVecCount,tosendVecPtr,MPI_DOUBLE,recvVecData,
                          remoteVecCount,remoteVecPtr,MPI_DOUBLE,MPI_COMM_WORLD);
            
            //*** Step 6 - Compute the results ***//
            result = (double *)malloc(sizeof(double) * vSize);
            computeMatVec(myid, vSize,myRowptr,myColInd,myMatVal,myVecData,result);
            
            //*** Step 7 - Gather results from all processors ***//
            totalResult = (double *)malloc(sizeof(double) * nrows);
            MPI_Gatherv(result,vSize,MPI_DOUBLE,totalResult,vecDataSize,vecDataDispls,
                        MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            

            // Write the result back to vector
            for (i = 0; i < nrows; i++) {
                allVecData[i] = totalResult[i];
                printf("NEW VEXX %d: %lf\n", i, allVecData[i]);
            }
            
            
            // Creating output filename
            int c = 0, temp = nrows;
            while (temp != 0) {
                temp = temp / 10;
                c++;
            }
            char fname[c+5];
            sprintf(fname,"o%d.vec",nrows);
            
            // Writing results to output file
            FILE *fpout;
            fpout = fopen(fname,"w");
            for (i = 0; i < nrows; i++)
                fprintf(fpout,"%.15lg\n",totalResult[i]);
            
            fclose(fpout);
            free(totalResult);
            
            // Other processes: Receive Data / Compute / Send Back
        } else {
            int *colInd, *eCount, *eDispls, *rowptr, *vecDataSize, *vecDataDispls;
            double *matVal;
            int size;
            
            MPI_Bcast(&nrows,1,MPI_INT,0,MPI_COMM_WORLD);
            vSize = ceil((double)nrows / num_procs);
            size = vSize;
            if (myid == num_procs-1)
                vSize = nrows - myid*vSize;
            // Prepare to store local vector entries
            myVecData = (double *)malloc(sizeof(double) * vSize);
            
            // Receive vector entries from root processor
            MPI_Scatterv(allVecData,vecDataSize,vecDataDispls,MPI_DOUBLE,
                         myVecData,vSize,MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            MPI_Scatter(eCount,1,MPI_INT,&myNumE,1,MPI_INT,0,MPI_COMM_WORLD);
            
            // Receive matrix entries from root processor
            myRowptr = (int *)malloc(sizeof(int) * (vSize+1));
            MPI_Scatterv(rowptr,vecDataSize,vecDataDispls,MPI_INT,
                         myRowptr,vSize,MPI_INT,0,MPI_COMM_WORLD);
            int rowptrOffset = myRowptr[0];
            for (i = 0; i < vSize; i++)
                myRowptr[i] = myRowptr[i] - rowptrOffset; //offset
            myRowptr[vSize] = myNumE;
            
            myColInd = (int *)malloc(sizeof(int) * myNumE);
            MPI_Scatterv(colInd,eCount,eDispls,MPI_INT,
                         myColInd,myNumE,MPI_INT,0,MPI_COMM_WORLD);
            
            myMatVal = (double *)malloc(sizeof(double) * myNumE);
            MPI_Scatterv(matVal,eCount,eDispls,MPI_DOUBLE,
                         myMatVal,myNumE,MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            // Determine the remote vectors needed, expand local Vector array to receive these.
            // Re-index local Column Index array to point to the right Vector Data entry
            numRemoteVec = prepareRemoteVec(myColInd,myNumE,myid,num_procs,nrows,vSize,remoteVecCount,remoteVecPtr,&remoteVecIndex,&myVecData);
            
            // Each process send to all others number of remote vector entries needed
            MPI_Alltoall(remoteVecCount,1,MPI_INT,tosendVecCount,1,MPI_INT,MPI_COMM_WORLD);
            
            for (p = 0; p < num_procs; p++)
                tosendTotal = tosendTotal + tosendVecCount[p];
            tosendVecIndex = (int *)malloc(sizeof(int) * tosendTotal);
            tosendVecPtr[0] = 0;
            for (p = 1; p < num_procs; p++)
                tosendVecPtr[p] = tosendVecPtr[p-1] + tosendVecCount[p-1];
            
            // Each process send to all others index of remote vector entries needed
            MPI_Alltoallv(remoteVecIndex,remoteVecCount,remoteVecPtr,MPI_INT,
                          tosendVecIndex,tosendVecCount,tosendVecPtr,MPI_INT,MPI_COMM_WORLD);
            // Get the vector entries other processors needed, ready to send
            sendVecData = (double *)malloc(sizeof(double) * tosendTotal);
            for (i = 0; i < tosendTotal; i++)
                sendVecData[i] = myVecData[tosendVecIndex[i]-myid*size];
            
            // Local Vector entries at the front, remote entries are concatenated starting from:
            double *recvVecData = myVecData + vSize;
            
            // All processors send the actual vector entries to each other
            MPI_Alltoallv(sendVecData,tosendVecCount,tosendVecPtr,MPI_DOUBLE,recvVecData,
                          remoteVecCount,remoteVecPtr,MPI_DOUBLE,MPI_COMM_WORLD);
            
            // Compute the local results
            result = (double *)malloc(sizeof(double) * vSize);
            computeMatVec(myid, vSize,myRowptr,myColInd,myMatVal,myVecData,result);
            
            // Send the results to the root processor
            MPI_Gatherv(result,vSize,MPI_DOUBLE,totalResult,vecDataSize,vecDataDispls,
                        MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    free(myRowptr); free(myColInd); free(myMatVal); free(tosendVecIndex);
    free(sendVecData); free(result); free(myVecData); free(remoteVecIndex);
    
    MPI_Finalize();
    return 0;
}

//CRS Matrix-Vector multiplication
void computeMatVec(int myid, int nrows, int *rowptr, int *colind, double *matval, double *vec, double *res) {
    int i, j, count;
    count = 0;
    for (i = 0; i < nrows; i++) {
        res[i] = 0.0;
        for (j = rowptr[i]; j < rowptr[i+1]; j++) {
            res[i] += matval[count] * vec[colind[j]];
            count++;
        }
    }
}

// Determine the remote vectors needed, expand local Vector array to receive these.
// Re-index local Column Index array to point to the right Vector Data entry
int prepareRemoteVec(int *myColInd, int myNumE, int myid, int numprocs, int nrows, int myVecSize,
                     int *remoteVecCount, int *remoteVecPtr, int **remoteVecIndex, double **myVecData) {
    int i;
    int size = ceil((double)nrows / numprocs);
    int localVecIndStart = myid * size;
    int localVecIndEnd = (myid + 1) * size;
    int tempRemoteVecCount[numprocs];
    
    for (i = 0; i < numprocs; i++)
        remoteVecCount[i] = 0;
    *remoteVecIndex = (int *)malloc(sizeof(int) * myNumE);
    
    // Select remote vectors (counts, indices, displacements from each remote processor)
    int numRemoteVec = 0;
    for (i = 0; i < myNumE; i++) {
        if ((myColInd[i] < localVecIndStart) || (myColInd[i] >= localVecIndEnd)) {
            remoteVecCount[myColInd[i]/size]++;
            numRemoteVec++;
        }
    }
    remoteVecPtr[0] = 0;
    tempRemoteVecCount[0] = 0;
    for (i = 1; i < numprocs; i++) {
        remoteVecPtr[i] = remoteVecPtr[i-1] + remoteVecCount[i-1];
        tempRemoteVecCount[i]=0;
    }
    
    // Realloc local vector array to store received remote vector entries
    *myVecData = (double *)realloc(*myVecData,sizeof(double) * (myVecSize+numRemoteVec));
    
    // Re-index the local Column Index array to point to local Vector Data entry at the beginning and remote Vector Data entries at the end of Vector Data array.
    int proc = 0;
    for (i = 0; i < myNumE; i++) {
        if ((myColInd[i] < localVecIndStart) || (myColInd[i] >= localVecIndEnd)) {
            proc = myColInd[i]/size;
            (*remoteVecIndex)[remoteVecPtr[proc]+tempRemoteVecCount[proc]] = myColInd[i];
            
            myColInd[i] = myVecSize + remoteVecPtr[proc] + tempRemoteVecCount[proc];
            tempRemoteVecCount[proc]++;
        }
        else {
            myColInd[i] -= (myid * size);
        }
    }
    
    //return total number of remote vector entries needed
    return numRemoteVec;
}
