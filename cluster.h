#ifndef CLUSTER_H
#define CLUSTER_H

#include <stdio.h>
#include <stdbool.h>
#include "matrix.h"
#include "ising.h"

void initializeClusterVars(Parameters *para)
{
    if (para->dim == 2)
    {
        para->cv.cluster = (void *)matrixMalloc2DBool(para->N, para->N);
    }
    else if (para->dim == 3)
    {
        para->cv.cluster = (void *)matrixMalloc3DBool(para->N, para->N, para->N);
    }
    
    para->cv.addProbability = 0; // muss später mit 1 - exp(-2*para->J/T) initialisiert werden
}

// vor jedem Schritt muss die Clustermatrix wieder zurückgesetzt werden
void resetClusterMatrix(Parameters *para)
{
    int i,j,k;
    if (para->dim == 2)
    {
        bool **mat = (bool **)para->cv.cluster;
        for (i = 0; i < para->N; i++)
            for (j = 0; j < para->N; j++)
                mat[i][j] = false;
    }
    else if (para->dim == 3)
    {
        bool ***mat = (bool ***)para->cv.cluster;
        for (i = 0; i < para->N; i++)
            for (j = 0; j < para->N; j++)
                for (k = 0; k < para->N; k++)
                    mat[i][j][k] = false;
    }
}


void growCluster(Parameters *para, int row, int col, int dep, int clusterSpin);


void tryAdd(Parameters *para, int row, int col, int dep, int clusterSpin)
{
    if (para->dim == 2)
    {
        int **spins2D = (int **)para->spins;        
        if (spins2D[row][col] == clusterSpin)
            if (mt_random()/ (double) MT_MAX < para->cv.addProbability)
            {
               growCluster(para, row, col, dep, clusterSpin);
            }
    }
    else if (para->dim == 3)
    {
        int ***spins3D = (int ***)para->spins;
        if (spins3D[row][col][dep] == clusterSpin)
            if (mt_random()/ (double) MT_MAX < para->cv.addProbability)
            {
                growCluster(para, row, col, dep, clusterSpin);
            }
    }
}

//wird aufgerufen, um den Aufbau des Clusters zu starten und ruft sich dann über "tryAdd" rekursiv auf
void growCluster(Parameters *para, int row, int col, int dep, int clusterSpin)
{
    if (para->dim == 2)
    {
        int **spins2D = (int **)para->spins;
        bool ** cluster = (bool **)para->cv.cluster;

        int rowPrev = (row == 0 ? para->N - 1 : row - 1);
        int rowNext = (row == para->N - 1 ? 0 : row + 1);
        int colPrev = (col == 0 ? para->N - 1 : col - 1);
        int colNext = (col == para->N - 1 ? 0 : col + 1);

        cluster[row][col] = true;
        spins2D[row][col] = -spins2D[row][col];
        
        if (!cluster[rowPrev][col])
            tryAdd(para, rowPrev, col, 0, clusterSpin);
        if (!cluster[rowNext][col])
            tryAdd(para, rowNext, col, 0, clusterSpin);
        if (!cluster[row][colPrev])
            tryAdd(para, row, colPrev, 0, clusterSpin);
        if (!cluster[row][colNext])
            tryAdd(para, row, colNext, 0, clusterSpin);
    }
    else if (para->dim == 3)
    {
        int ***spins3D = (int ***)para->spins;
        bool ***cluster = (bool ***)para->cv.cluster;
    
        int rowPrev = (row == 0 ? para->N - 1 : row - 1);
        int rowNext = (row == para->N - 1 ? 0 : row + 1);
        int colPrev = (col == 0 ? para->N - 1 : col - 1);
        int colNext = (col == para->N - 1 ? 0 : col + 1);
        int depPrev = (dep == 0 ? para->N - 1 : dep - 1);
        int depNext = (dep == para->N - 1 ? 0 : dep + 1);
     
        cluster[row][col][dep] = true;
        spins3D[row][col][dep] = -spins3D[row][col][dep];
   
        if (!cluster[rowPrev][col][dep])
            tryAdd(para, rowPrev, col, dep, clusterSpin);
        if (!cluster[rowNext][col][dep])
            tryAdd(para, rowNext, col, dep, clusterSpin);
        if (!cluster[row][colPrev][dep])
            tryAdd(para, row, colPrev, dep, clusterSpin);
        if (!cluster[row][colNext][dep])
            tryAdd(para, row, colNext, dep, clusterSpin);
        if (!cluster[row][col][depPrev])
            tryAdd(para, row, col, depPrev, clusterSpin);
        if (!cluster[row][col][depNext])
            tryAdd(para, row, col, depNext, clusterSpin);
    }
} 

void clusterDelete(Parameters *para)
{
    if(para->dim == 2)
    {
        bool **mat = (bool **)para->cv.cluster;
        free(*mat);
        free(mat);
    }
    else if (para->dim == 3)
    {
        bool ***mat = (bool ***)para->cv.cluster;
        free(**mat);
        free(*mat);
        free(mat);
    }
}

#endif










