#ifndef CLUSTER_H
#define CLUSTER_H

#include <stdio.h>
#include <stdbool.h>
#include "matrix.h"
#include "ising.h"

void initializeClusterVars(Parameters *para)
{
    if (para->dim == 3)
    {
        fprintf(stderr, "3-Dimensional use of Cluster-Algorithm isn't supported");
        exit(EXIT_FAILURE);
    }   
  
    int i;
    bool  *data;   
    
    data = malloc(para->N*para->N*sizeof(*data));
    para->cv.cluster = malloc(para->N*sizeof(*para->cv.cluster));
    
    for (i=0; i<para->N; ++i)
    {
        para->cv.cluster[i] = data + i*para->N;
    } 
    
    para->cv.addProbability = 0; // muss später mit 1 - exp(-2*para->J/T) initialisiert werden
}

// vor jedem Schritt muss die Clustermatrix wieder zurückgesetzt werden
void resetClusterMatrix(Parameters *para)
{
    int i,j;
    for (i = 0; i < para->N; i++)
        for (j = 0; j < para->N; j++)
            para->cv.cluster[i][j] = false;
}


void growCluster(Parameters *para, int row, int col, int clusterSpin);


void tryAdd(Parameters *para, int row, int col, int clusterSpin) {
    int **spins2D = (int **)para->spins;
    
    if (spins2D[row][col] == clusterSpin)
        if (mt_random()/ (double) MT_MAX < para->cv.addProbability)
        {
           growCluster(para, row, col, clusterSpin);
        }
}

//wird aufgerufen, um den Aufbau des Clusters zu starten und ruft sich dann über "tryAdd" rekursiv auf
void growCluster(Parameters *para, int row, int col, int clusterSpin)
{
    printf("%f \n", para->cv.addProbability);
    int **spins2D = (int **)para->spins;

    int rowPrev = (row == 0 ? para->N - 1 : row - 1);
    int rowNext = (row == para->N - 1 ? 0 : row + 1);
    int colPrev = (col == 0 ? para->N - 1 : col - 1);
    int colNext = (col == para->N - 1 ? 0 : col + 1);

    para->cv.cluster[row][col] = true;
    spins2D[row][col] = -spins2D[row][col];
    
    if (!para->cv.cluster[rowPrev][col])
        tryAdd(para, rowPrev, col, clusterSpin);
    if (!para->cv.cluster[rowNext][col])
        tryAdd(para, rowNext, col, clusterSpin);
    if (!para->cv.cluster[row][colPrev])
        tryAdd(para, row, colPrev, clusterSpin);
    if (!para->cv.cluster[row][colNext])
        tryAdd(para, row, colNext, clusterSpin);
} 

void clusterDelete(Parameters *para)
{
    free(*para->cv.cluster);
    free(para->cv.cluster);
}

#endif










